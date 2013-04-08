////////////////////////////////////////////////////////////////////////
//
//  CosmicTracker
//
//  Tracker to reconstruct cosmic ray muons
// 
//  tjyang@fnal.gov
// 
//  This algorithm is based on orignal idea in Track3Dreco
// 
////////////////////////////////////////////////////////////////////////

// C++ includes
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

// ROOT includes
#include "TVectorD.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TH1D.h"
#include "TVirtualFitter.h"

//2-D weighted fit of time vs wire
std::vector<double> vwire;
std::vector<double> vtime;
std::vector<double> vph;

void myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
  //minimisation function computing the sum of squares of residuals
  f = 0;
  for (size_t i = 0; i<vwire.size(); ++i){
    double x = vwire[i];
    double y = par[0]+par[1]*x+par[2]*x*x;
    //using ph^2 to suppress low ph hits
    f += vph[i]*vph[i]*(vtime[i]-y)*(vtime[i]-y);
  }
}

namespace trkf {
   
  struct SortByWire {
    bool operator() (art::Ptr<recob::Hit> const& h1, art::Ptr<recob::Hit> const& h2) const { 
      return 
	h1->Wire()->RawDigit()->Channel() < 
	h2->Wire()->RawDigit()->Channel() ;
    }
  };

  class CosmicTracker : public art::EDProducer {
    
  public:
    
    explicit CosmicTracker(fhicl::ParameterSet const& pset);
    ~CosmicTracker();
    
    //////////////////////////////////////////////////////////
    void reconfigure(fhicl::ParameterSet const& p);
    void produce(art::Event& evt); 
    void beginJob();
    void endJob();

  private:

    double          fKScut;              ///< tolerance for cluster matching based on KS test.

    double          ftmatch;             ///< tolerance for time matching (in time samples) 
    
    double          fsmatch;             ///< tolerance for distance matching (in cm)

    std::string     fClusterModuleLabel; ///< label for input cluster collection

    //testing histograms
    TH1D *dt[3];
    TH1D *dtime[3];
    TH1D *testsig[3];
    TH1D *hsig[3];
  
  }; // class CosmicTracker

}

namespace trkf {

//-------------------------------------------------
CosmicTracker::CosmicTracker(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
  produces< std::vector<recob::Track>                        >();
  produces< std::vector<recob::SpacePoint>                   >();
  produces< art::Assns<recob::Track,      recob::Cluster>    >();
  produces< art::Assns<recob::Track,      recob::SpacePoint> >();
  produces< art::Assns<recob::SpacePoint, recob::Hit>        >();
  produces< art::Assns<recob::Track,      recob::Hit>        >();
}

//-------------------------------------------------
CosmicTracker::~CosmicTracker()
{
}

void CosmicTracker::reconfigure(fhicl::ParameterSet const& pset)
{
  fClusterModuleLabel     = pset.get< std::string >("ClusterModuleLabel");
  fKScut                  = pset.get< double >("KScut");
  ftmatch                 = pset.get< double >("TMatch");
  fsmatch                 = pset.get< double >("SMatch");
}

//-------------------------------------------------
void CosmicTracker::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  dt[0] = tfs->make<TH1D>("dt0","dt0",100,0,100);
  dt[1] = tfs->make<TH1D>("dt1","dt1",100,0,100);
  dt[2] = tfs->make<TH1D>("dt2","dt2",100,0,100);

  dtime[0] = tfs->make<TH1D>("dtime0","dtime0",100,-50,50);
  dtime[1] = tfs->make<TH1D>("dtime1","dtime1",100,-50,50);
  dtime[2] = tfs->make<TH1D>("dtime2","dtime2",100,-50,50);

  testsig[0] = tfs->make<TH1D>("testsig0","testsig0",4096,0,4096);
  testsig[1] = tfs->make<TH1D>("testsig1","testsig1",4096,0,4096);
  testsig[2] = tfs->make<TH1D>("testsig2","testsig2",4096,0,4096);

  hsig[0] = tfs->make<TH1D>("hsig0","hsig0",4096,0,4096);
  hsig[1] = tfs->make<TH1D>("hsig1","hsig1",4096,0,4096);
  hsig[2] = tfs->make<TH1D>("hsig2","hsig2",4096,0,4096);

  for (int i = 0; i<3; ++i) dtime[i]->Sumw2();


}

void CosmicTracker::endJob()
{
}

//------------------------------------------------------------------------------------//
void CosmicTracker::produce(art::Event& evt){
  
  // get services
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;
  art::ServiceHandle<util::DetectorProperties> detprop;

  std::unique_ptr<std::vector<recob::Track>      >              tcol (new std::vector<recob::Track>);	   
  std::unique_ptr<std::vector<recob::SpacePoint> > 	        spcol(new std::vector<recob::SpacePoint>);
  std::unique_ptr<art::Assns<recob::Track, recob::SpacePoint> > tspassn(new art::Assns<recob::Track, recob::SpacePoint>);
  std::unique_ptr<art::Assns<recob::Track, recob::Cluster> >    tcassn(new art::Assns<recob::Track, recob::Cluster>);
  std::unique_ptr<art::Assns<recob::Track, recob::Hit> >        thassn(new art::Assns<recob::Track, recob::Hit>);
  std::unique_ptr<art::Assns<recob::SpacePoint, recob::Hit> >   shassn(new art::Assns<recob::SpacePoint, recob::Hit>);

  /*
  //print wire ends
  double xyzStart[3];
  double xyzEnd[3];
  for (size_t ip = 0; ip<geom->Nplanes(); ++ip){
    for (size_t iw = 0; iw<geom->Nwires(ip); ++iw){
      geom->WireEndPoints(0,0,ip,iw,xyzStart,xyzEnd);
      std::cout<<ip<<" "<<iw<<" "<<xyzStart[0]<<" "<<xyzStart[1]<<" "<<xyzStart[2]<<
	" "<<xyzEnd[0]<<" "<<xyzEnd[1]<<" "<<xyzEnd[2]<<std::endl;
    }
  }
  */

  double timetick = detprop->SamplingRate()*1e-3;    //time sample in us
  double presamplings = detprop->TriggerOffset(); // presamplings in ticks  
  double plane_pitch = geom->PlanePitch(0,1);   //wire plane pitch in cm 
  double wire_pitch = geom->WirePitch(0,1,0);    //wire pitch in cm
  double Efield_drift = larprop->Efield(0);  // Electric Field in the drift region in kV/cm
  double Temperature = larprop->Temperature();  // LAr Temperature in K

  double driftvelocity = larprop->DriftVelocity(Efield_drift,Temperature);    //drift velocity in the drift region (cm/us)
  double timepitch = driftvelocity*timetick;                         //time sample (cm) 

  int nts = detprop->NumberTimeSamples();
  int nplanes = geom->Nplanes();

  std::vector<TH1D*> signals[nplanes];

   // get input Cluster object(s).
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
      art::fill_ptr_vector(clusterlist, clusterListHandle);

  art::FindManyP<recob::Hit> fm(clusterListHandle, evt, fClusterModuleLabel);

  std::vector<int> Cls[nplanes];

  for (size_t iclu = 0; iclu<clusterlist.size(); ++iclu){

    double t0 = clusterlist[iclu]->StartPos()[1];
    double t1 = clusterlist[iclu]->EndPos()[1];
    t0 -= detprop->GetXTicksOffset(clusterlist[iclu]->View(),0,0);
    t1 -= detprop->GetXTicksOffset(clusterlist[iclu]->View(),0,0);
    
    switch(clusterlist[iclu]->View()){
    case geo::kU :
      Cls[0].push_back(iclu);
      break;
    case geo::kV :
      Cls[1].push_back(iclu);
      break;
    case geo::kW :
      Cls[2].push_back(iclu);
      break;
    default :
      break;
    }

  }

  //calibrate drift times between wire planes using single muons
  std::vector<double> meantime[nplanes];
  for (int i = 0; i<nplanes; ++i){
    for (size_t ic = 0; ic<Cls[i].size(); ++ic){
      TH1D *sig = new TH1D(Form("sig_%d_%d",i,int(ic)),Form("sig_%d_%d",i,int(ic)),nts,0,nts);
      TH1D *sigint = new TH1D(Form("sigint_%d_%d",i,int(ic)),Form("sigint_%d_%d",i,int(ic)),nts,0,nts);      
      std::vector< art::Ptr<recob::Hit> > hitlist = fm.at(Cls[i][ic]);
      std::sort(hitlist.begin(), hitlist.end(), trkf::SortByWire());
      for(auto theHit = hitlist.begin(); theHit != hitlist.end();  theHit++){
	
	double time = (*theHit)->PeakTime();
	time -= detprop->GetXTicksOffset((*theHit)->WireID().Plane,
					 (*theHit)->WireID().TPC,
					 (*theHit)->WireID().Cryostat);

	double charge = (*theHit)->Charge();
	int bin = sig->FindBin(time);
	sig->SetBinContent(bin,sig->GetBinContent(bin)+charge);
	for (int j = bin; j<=sig->GetNbinsX(); ++j){
	  sigint->SetBinContent(j,sigint->GetBinContent(j)+charge);
	}
      }
      if (sigint->Integral()) sigint->Scale(1./sigint->GetBinContent(sigint->GetNbinsX()));
      signals[i].push_back(sigint);
      if (hitlist.size()>10){
	meantime[i].push_back(sig->GetMean());
      }
      delete sig;
    }
  }

  bool singletrack = true;
  for (int i = 0; i<nplanes; ++i){
    singletrack = singletrack&&meantime[i].size()==1;
  }
  if (singletrack){
    for (int i = 0; i<nplanes; ++i){
      for (int j = i+1; j<nplanes; ++j){
	dtime[i+j-1]->Fill(meantime[j][0]-meantime[i][0]);
      }
    }
    for (int i = 0; i<nplanes; ++i){
      for (size_t k = 0; k<signals[i].size(); ++k){
	if (fm.at(Cls[i][k]).size()<10) continue;
	for (int j = 0; j<signals[i][k]->GetNbinsX(); ++j){
	  double binc = signals[i][k]->GetBinContent(j+1);
	  testsig[i]->SetBinContent(j+1,binc);
	}
      }
    }
  }


  //matching clusters between different views
  std::vector<int> matched(clusterlist.size());
  for (size_t i = 0; i<clusterlist.size(); ++i) matched[i] = 0;

  std::vector< std::vector<int> > matchedclusters;

//  for (int i = 0; i<nplanes-1; ++i){
//    for (int j = i+1; j<nplanes; ++j){
  for (int i = 0; i<nplanes; ++i){
    for (int j = 0; j<nplanes; ++j){
      for (size_t c1 = 0; c1<Cls[i].size(); ++c1){
	for (size_t c2 = 0; c2<Cls[j].size(); ++c2){
	  // check if both are the same view
	  if (clusterlist[Cls[i][c1]]->View()==
	      clusterlist[Cls[j][c2]]->View()) continue;
	  // check if both are already in the matched list
	  if (matched[Cls[i][c1]]==1&&matched[Cls[j][c2]]==1) continue;
	  // KS test between two views in time
	  double ks = signals[i][c1]->KolmogorovTest(signals[j][c2]);
	  int imatch = -1; //track candidate index
	  int iadd = -1; //cluster index to be inserted
	  if (ks>fKScut){//pass KS test
	    // check both clusters with all matched clusters
	    // if one is already matched, 
	    // check if need to add the other to the same track candidate
	    for (size_t l = 0; l<matchedclusters.size(); ++l){
	      for (size_t m = 0; m<matchedclusters[l].size(); ++m){
		if (matchedclusters[l][m]==Cls[i][c1]){
		  imatch = l; //track candidate
		  iadd = j; //consider the other cluster
		}
		else if (matchedclusters[l][m]==Cls[j][c2]){
		  imatch = l; //track candidate
		  iadd = i; //consider the other cluster
		}
	      }
	    }
	    if (imatch>=0){
	      if (iadd == i){
		bool matchview = false;
		// check if one matched cluster has the same view
		for (size_t ii = 0; ii<matchedclusters[imatch].size(); ++ii){
		  if (clusterlist[matchedclusters[imatch][ii]]->View()==
		      clusterlist[Cls[i][c1]]->View()){
		    matchview = true;
		    //replace if the new cluster has more hits
		    if (fm.at(Cls[i][c1])>fm.at(matchedclusters[imatch][ii])){
		      matched[matchedclusters[imatch][ii]] = 0;
		      matchedclusters[imatch][ii] = Cls[i][c1];
		      matched[Cls[i][c1]] = 1;
		    }
		  }
		}
		if (!matchview){//not matched view found, just add
		  matchedclusters[imatch].push_back(Cls[i][c1]);
		  matched[Cls[i][c1]] = 1;
		}
	      }
	      else {
		bool matchview = false;
		for (size_t jj = 0; jj<matchedclusters[imatch].size(); ++jj){
		  if (clusterlist[matchedclusters[imatch][jj]]->View()==
		      clusterlist[Cls[j][c2]]->View()){
		    matchview = true;
		    //replace if it has more hits
		    if (fm.at(Cls[j][c2])>fm.at(matchedclusters[imatch][jj])){
		      matched[matchedclusters[imatch][jj]] = 0;
		      matchedclusters[imatch][jj] = Cls[j][c2];
		      matched[Cls[j][c2]] = 1;
		    }
		  }
		}
		if (!matchview){
		  matchedclusters[imatch].push_back(Cls[j][c2]);
		  matched[Cls[j][c2]] = 1;
		}		
	      }
	    }
	    else{
	      std::vector<int> tmp;
	      tmp.push_back(Cls[i][c1]);
	      tmp.push_back(Cls[j][c2]);
	      matchedclusters.push_back(tmp);
	    }
	    matched[Cls[i][c1]]=1;
	    matched[Cls[j][c2]]=1;
	  }//pass KS test
	}//c2
      }//c1
    }//j
  }//i

  for (size_t i = 0; i<matchedclusters.size(); ++i){
    if (matchedclusters[i].size()) mf::LogVerbatim("CosmicTracker")<<"Track candidate "<<i<<":";
    for (size_t j = 0; j<matchedclusters[i].size(); ++j){
      mf::LogVerbatim("CosmicTracker")<<matchedclusters[i][j];
    }
  } 

  for (int i = 0; i<nplanes; ++i){
    for (size_t j = 0; j<signals[i].size(); ++j){
      delete signals[i][j];
    }
  }

  /////////////////////////////////////////////////////
  /////// 2D Track Matching and 3D Track Reconstruction
  /////////////////////////////////////////////////////
  
  ///Prepare fitter
  TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
  TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
  fitter->SetFCN(myfcn);
  double arglist[10];
  arglist[0] = -1;
  fitter->ExecuteCommand("SET PRIN",arglist,1);

  //fit each cluster in 2D using pol2, iterate once to remove outliers
  for (size_t itrk = 0; itrk<matchedclusters.size(); ++itrk){//loop over tracks

    //all the clusters associated with the current track
    art::PtrVector<recob::Cluster> clustersPerTrack;
    for (size_t iclu = 0; iclu<matchedclusters[itrk].size(); ++iclu){
      art::Ptr <recob::Cluster> cluster(clusterListHandle,matchedclusters[itrk][iclu]);
      clustersPerTrack.push_back(cluster);
    }

    //save time/hit information along track trajectory
    std::vector<std::map<int,double> > vtimemap;
    std::vector<std::map<int,art::Ptr<recob::Hit> > > vhitmap;

    for (size_t iclu = 0; iclu<matchedclusters[itrk].size(); ++iclu){//loop over clusters

      vwire.clear();
      vtime.clear();
      vph.clear();
      //fit hits time vs wire with pol2
      std::vector< art::Ptr<recob::Hit> > hits = fm.at(matchedclusters[itrk][iclu]);
      std::sort(hits.begin(), hits.end(), trkf::SortByWire());
      double dtdw = 0;
      if (clusterlist[matchedclusters[itrk][iclu]]->StartPos()[1]-
	  clusterlist[matchedclusters[itrk][iclu]]->EndPos()[1]){
	dtdw = (clusterlist[matchedclusters[itrk][iclu]]->EndPos()[0]-
		clusterlist[matchedclusters[itrk][iclu]]->StartPos()[0])/
	  (clusterlist[matchedclusters[itrk][iclu]]->EndPos()[1]-
	   clusterlist[matchedclusters[itrk][iclu]]->StartPos()[1]);
      }
      fitter->SetParameter(0,"p0",clusterlist[matchedclusters[itrk][iclu]]->StartPos()[1]-dtdw-detprop->GetXTicksOffset(clusterlist[matchedclusters[itrk][iclu]]->View(),0,0),0.1,0,0);
      fitter->SetParameter(1,"p1",clusterlist[matchedclusters[itrk][iclu]]->dTdW(),0.1,0,0);
      fitter->SetParameter(2,"p2",0,0.1,0,0);
      
      for (size_t ihit = 0; ihit<hits.size(); ++ihit){//loop over hits
	geo::WireID hitWireID = hits[ihit]->WireID();
	unsigned int w = hitWireID.Wire;
	vwire.push_back(w);
	double time = hits[ihit]->PeakTime();
	time -= detprop->GetXTicksOffset(hits[ihit]->WireID().Plane,
					 hits[ihit]->WireID().TPC,
					 hits[ihit]->WireID().Cryostat);
	vtime.push_back(time);
	vph.push_back(hits[ihit]->Charge());
      }
      arglist[0] = 0;
      if (vwire.size()>2) fitter->ExecuteCommand("MIGRAD", arglist, 0);
      else{
	fitter->SetParameter(0,"p0",vtime[0],0.1,0,0);
	fitter->SetParameter(1,"p1",0,0.1,0,0);
	fitter->SetParameter(2,"p2",0,0.1,0,0);
      }
      //remove outliers
      for (auto iw = vwire.begin(), it = vtime.begin(), iph = vph.begin(); iw!=vwire.end(); ){
	double y = fitter->GetParameter(0)+
	  fitter->GetParameter(1)*(*iw)+
	  fitter->GetParameter(2)*(*iw)*(*iw);
	if (std::abs(*it-y)>30){
	  iw = vwire.erase(iw);
	  it = vtime.erase(it);
	  iph = vph.erase(iph);
	}
	else{
	  ++iw;
	  ++it;
	  ++iph;
	}
      }
      
      //refit
      if (vwire.size()>2) fitter->ExecuteCommand("MIGRAD", arglist, 0);

      std::map<int,double> timemap;
      std::map<int,double> phmap;
      std::map<int,art::Ptr<recob::Hit> > hitmap;

      //find hit on each wire along the fitted line
      for (size_t ihit = 0; ihit<hits.size(); ++ihit){//loop over hits
	geo::WireID hitWireID = hits[ihit]->WireID();
	unsigned int w = hitWireID.Wire;
	vwire.push_back(w);
	double time = hits[ihit]->PeakTime();
	time -= detprop->GetXTicksOffset(hits[ihit]->WireID().Plane,
					 hits[ihit]->WireID().TPC,
					 hits[ihit]->WireID().Cryostat);
	double ph = hits[ihit]->Charge();
	
	if (ph>(phmap[w])){
	  double y = fitter->GetParameter(0)+
	    fitter->GetParameter(1)*w+
	    fitter->GetParameter(2)*w*w;
	  if (std::abs(time-y)<20){
	    phmap[w] = ph;
	    timemap[w] = time;
	    hitmap[w] = hits[ihit];
	  }
	}
      }//ihit
      vtimemap.push_back(timemap);
      vhitmap.push_back(hitmap);
    }//iclu
    
     /// Find two clusters with the most numbers of hits, and time ranges
    int iclu1 = -1;
    int iclu2 = -1;
    int iclu3 = -1;
    unsigned maxnumhits0 = 0;
    unsigned maxnumhits1 = 0;
    
    double tmin[vtimemap.size()];
    double tmax[vtimemap.size()];
    for (size_t iclu = 0; iclu<vtimemap.size(); ++iclu){
      tmin[iclu] = 1e9;
      tmax[iclu] = -1e9;
    }
    
    for (size_t iclu = 0; iclu<vtimemap.size(); ++iclu){
      for (auto itime = vtimemap[iclu].begin(); itime!=vtimemap[iclu].end(); ++itime){
	if (itime->second>tmax[iclu]){
	  tmax[iclu] = itime->second;
	}
	if (itime->second<tmin[iclu]){
	  tmin[iclu] = itime->second;
	}
      }
      if (vtimemap[iclu].size()>maxnumhits0){
	if (iclu1!=-1){
	  iclu2 = iclu1;
	  maxnumhits1 = maxnumhits0;
	}
	iclu1 = iclu;
	maxnumhits0 = vtimemap[iclu].size();
      }
      else if (vtimemap[iclu].size()>maxnumhits1){
	iclu2 = iclu;
	maxnumhits1 = vtimemap[iclu].size();
      }
    }
    
    std::swap(iclu1,iclu2); //now iclu1 has fewer hits than iclu2

    for (int iclu = 0; iclu<(int)vtimemap.size(); ++iclu){
      if (iclu!=iclu1&&iclu!=iclu2) iclu3 = iclu;
    }
    
    if (iclu1!=-1&&iclu2!=-1){
      //select hits in a common time range
      auto ihit = vhitmap[iclu1].begin();
      auto itime = vtimemap[iclu1].begin();
      while (itime!=vtimemap[iclu1].end()){
	if (itime->second<std::max(tmin[iclu1],tmin[iclu2])-ftmatch||
	    itime->second>std::min(tmax[iclu1],tmax[iclu2])+ftmatch){
	  vtimemap[iclu1].erase(itime++);
	  vhitmap[iclu1].erase(ihit++);
	}
	else{
	  ++itime;
	  ++ihit;
	}
      }

      ihit = vhitmap[iclu2].begin();
      itime = vtimemap[iclu2].begin();
      while (itime!=vtimemap[iclu2].end()){
	if (itime->second<std::max(tmin[iclu1],tmin[iclu2])-ftmatch||
	    itime->second>std::min(tmax[iclu1],tmax[iclu2])+ftmatch){
	  vtimemap[iclu2].erase(itime++);
	  vhitmap[iclu2].erase(ihit++);
	}
	else{
	  ++itime;
	  ++ihit;
	}
      }
      
      //if one cluster is empty, replace it with iclu3
      if (!vtimemap[iclu1].size()){
	if (iclu3!=-1){
	  std::swap(iclu3,iclu1);
	}
      }
      if (!vtimemap[iclu2].size()){
	if (iclu3!=-1){
	  std::swap(iclu3,iclu2);
	  std::swap(iclu1,iclu2);
	}
      }
      if ((!vtimemap[iclu1].size())||(!vtimemap[iclu2].size())) continue;

      size_t spStart = spcol->size();
      std::vector<recob::SpacePoint> spacepoints;
      TVector3 startpointVec,endpointVec, DirCos;

      bool rev = false;
      auto times1 = vtimemap[iclu1].begin();
      auto timee1 = vtimemap[iclu1].end();
      --timee1;
      auto times2 = vtimemap[iclu2].begin();
      auto timee2 = vtimemap[iclu2].end();
      --timee2;

      double ts1 = times1->second;
      double te1 = timee1->second;
      double ts2 = times2->second;
      double te2 = timee2->second;
      
      auto hit1s = vhitmap[iclu1].begin();
      auto hit1e = vhitmap[iclu1].end();
      --hit1e;
      
      auto hit2s = vhitmap[iclu2].begin();
      auto hit2e = vhitmap[iclu2].end();
      --hit2e;
      
      //find out if we need to flip ends
      if (std::abs(ts1-ts2)+std::abs(te1-te2)>std::abs(ts1-te2)+std::abs(te1-ts2)){
	rev = true;
	hit2s = vhitmap[iclu2].end();
	--hit2s;
	times2 = vtimemap[iclu2].end();
	--times2;
	hit2e = vhitmap[iclu2].begin();
	timee2 = vtimemap[iclu2].begin();
      }
      double y,z;
      geom->ChannelsIntersect((hit1s->second)->Wire()->RawDigit()->Channel(),
			      (hit2s->second)->Wire()->RawDigit()->Channel(),
			      y,z);	     
      startpointVec.SetXYZ((ts1-presamplings)*timepitch+2*plane_pitch,y,z);
      geom->ChannelsIntersect((hit1e->second)->Wire()->RawDigit()->Channel(),
			      (hit2e->second)->Wire()->RawDigit()->Channel(),
			      y,z);	     
      endpointVec.SetXYZ((te1-presamplings)*timepitch+2*plane_pitch,y,z);
      
      DirCos = endpointVec - startpointVec;
      //SetMag casues a crash if the magnitude of the vector is zero
      try
	{
	  DirCos.SetMag(1.0);//normalize vector
	}
      catch(...){std::cout<<"The Spacepoint is infinitely small"<<std::endl;
	continue;
      }
      std::vector<double> vtracklength;
      
      for (size_t iclu = 0; iclu<vtimemap.size(); ++iclu){
	
	double tracklength = 0;
	if (vtimemap[iclu].size()==1){
	  tracklength = wire_pitch;
	}
	else{
	  double t0, w0;
	  for (auto iw = vtimemap[iclu].begin(); iw!=vtimemap[iclu].end(); ++iw){
	    if (iw==vtimemap[iclu].begin()){
	      w0 = iw->first;
	      t0 = iw->second;
	    }
	    else{
	      tracklength += std::sqrt(std::pow((iw->first-w0)*wire_pitch,2)+std::pow((iw->second-t0)*timepitch,2));
	      w0 = iw->first;
	      t0 = iw->second;	     
	    }
	    hsig[iclu]->Fill(iw->second);
	  }
	}
	vtracklength.push_back(tracklength);
      }
      
      std::map<int,int> maxhitsMatch;

      auto ihit1 = vhitmap[iclu1].begin();
      for (auto itime1 = vtimemap[iclu1].begin(); itime1!=vtimemap[iclu1].end(); ++itime1, ++ihit1){//loop over min-hits
	art::PtrVector<recob::Hit> sp_hits;
	sp_hits.push_back(ihit1->second);
	double hitcoord[3];
	double length1 = 0;
	if (vtimemap[iclu1].size()==1){
	  length1 = wire_pitch;
	}
	else{
	  for (auto iw1 = vtimemap[iclu1].begin(); iw1!=itime1; ++iw1){
	    auto iw2 = iw1;
	    ++iw2;
	    length1 += std::sqrt(std::pow((iw1->first-iw2->first)*wire_pitch,2)+std::pow((iw1->second-iw2->second)*timepitch,2));
	  }
	}
	double difference = 1e10; //distance between two matched hits
	auto matchedtime = vtimemap[iclu2].end();
	auto matchedhit  = vhitmap[iclu2].end();
	
	auto ihit2 = vhitmap[iclu2].begin();
	for (auto itime2 = vtimemap[iclu2].begin(); itime2!=vtimemap[iclu2].end(); ++itime2, ++ihit2){//loop over max-hits
	  if (maxhitsMatch[itime2->first]) continue;
	  double length2 = 0;
	  if (vtimemap[iclu2].size()==1){
	    length2 = wire_pitch;
	  }
	  else{
	    for (auto iw1 = vtimemap[iclu2].begin(); iw1!=itime2; ++iw1){
	      auto iw2 = iw1;
	      ++iw2;
	      length2 += std::sqrt(std::pow((iw1->first-iw2->first)*wire_pitch,2)+std::pow((iw1->second-iw2->second)*timepitch,2));
	    }
	  }
	  if (rev) length2 = vtracklength[iclu2] - length2;
	  length2 = vtracklength[iclu1]/vtracklength[iclu2]*length2;
	  bool timematch = std::abs(itime1->second-itime2->second)<ftmatch;
	  if (timematch &&std::abs(length2-length1)<difference){
	    difference = std::abs(length2-length1);
	    matchedtime = itime2;
	    matchedhit = ihit2;
	  }
	}//loop over hits2
	if (difference<1){
	  hitcoord[0] = matchedtime->second*detprop->GetXTicksCoefficient();
	  hitcoord[1] = -1e10;
	  hitcoord[2] = -1e10;
	  geom->ChannelsIntersect((ihit1->second)->Wire()->RawDigit()->Channel(),
				  (matchedhit->second)->Wire()->RawDigit()->Channel(),
				  hitcoord[1],hitcoord[2]);
	  if (hitcoord[1]>-1e9&&hitcoord[2]>-1e9){
	    maxhitsMatch[matchedtime->first] = 1;
	    sp_hits.push_back(matchedhit->second);
	  }
	}
	if (sp_hits.size()>1){
	  double err[6] = {util::kBogusD};
	  recob::SpacePoint mysp(hitcoord, err, util::kBogusD, spStart + spacepoints.size());//3d point at end of track
	  spacepoints.push_back(mysp);
	  spcol->push_back(mysp);	
	  util::CreateAssn(*this, evt, *spcol, sp_hits, *shassn);
	}
      }//loop over hits1
      
      size_t spEnd = spcol->size();
      
      // Add the 3D track to the vector of the reconstructed tracks
      if(spacepoints.size()>0){
	
	// make a vector of the trajectory points along the track
	std::vector<TVector3> xyz(spacepoints.size());
	for(size_t s = 0; s < spacepoints.size(); ++s){
	  xyz[s] = TVector3(spacepoints[s].XYZ());
	}
	
	///\todo really should fill the direction cosines with unique values 
	std::vector<TVector3> dircos(spacepoints.size(), DirCos);
	
	std::vector< std::vector<double> > dQdx;
	std::vector<double> mom(2, util::kBogusD);
	tcol->push_back(recob::Track(xyz, dircos, dQdx, mom, tcol->size()));
	
	// make associations between the track and space points
	util::CreateAssn(*this, evt, *tcol, *spcol, *tspassn, spStart, spEnd);
	
	// now the track and clusters
	util::CreateAssn(*this, evt, *tcol, clustersPerTrack, *tcassn);
	
	// and the hits and track
	art::FindManyP<recob::Hit> fmh(clustersPerTrack, evt, fClusterModuleLabel);
	for(size_t cpt = 0; cpt < clustersPerTrack.size(); ++cpt)
	  util::CreateAssn(*this, evt, *tcol, fmh.at(cpt), *thassn);
	
	/*
	std::vector<art::Ptr<recob::Hit> > trkhits;
	for(auto ihit = vhitmap[iclu1].begin(); ihit!=vhitmap[iclu1].end();++ihit){
	  trkhits.push_back(ihit->second);
	}
	for(auto ihit = vhitmap[iclu2].begin(); ihit!=vhitmap[iclu2].end();++ihit){
	  trkhits.push_back(ihit->second);
	}
	util::CreateAssn(*this, evt, *tcol, trkhits, *thassn);
	*/
      }
    }//if iclu1&&iclu2
  }//itrk		 

  mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
  mf::LogVerbatim("Summary") << "CosmicTracker Summary:";
  for(unsigned int i = 0; i<tcol->size(); ++i) mf::LogVerbatim("Summary") << tcol->at(i) ;
  
  evt.put(std::move(tcol));
  evt.put(std::move(spcol));
  evt.put(std::move(tspassn));
  evt.put(std::move(tcassn));
  evt.put(std::move(thassn));
  evt.put(std::move(shassn));

  return;
}

  DEFINE_ART_MODULE(CosmicTracker);

} // namespace
