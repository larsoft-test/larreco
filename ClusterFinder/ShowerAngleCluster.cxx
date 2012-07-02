////////////////////////////////////////////////////////////////////////
//
// \file ShowerAngleCluster.cxx
//
// biagio.rossi@lhep.unibe.ch   (FWMK : argoneut specific)
// thomas.strauss@lhep.unibe.ch (ART  : general detector)
//
// andrzej.szelc@yale.edu (port to detector agnostic version)
//
// This algorithm is designed to reconstruct showers
// 
///////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

// Framework includes
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

// LArSoft includes
#include "Simulation/sim.h"
#include "ClusterFinder/ShowerAngleCluster.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"
#include "Utilities/AssociationUtil.h"

#include "SimulationBase/simbase.h"
#include "RawData/RawDigit.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "SummaryData/summary.h"


// ***************** //

//------------------------------------------------------------------------------
cluster::ShowerAngleCluster::ShowerAngleCluster(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
  produces< std::vector<recob::Cluster> >();
  produces< art::Assns<recob::Cluster, recob::Hit>  >();
   
}


void cluster::ShowerAngleCluster::reconfigure(fhicl::ParameterSet const& pset) 
{
  fClusterModuleLabel 		=pset.get< std::string >("ClusterModuleLabel");
  fVertexCLusterModuleLabel	=pset.get<std::string > ("VertexClusterModuleLabel");
  fMCGeneratorLabel		=pset.get<std::string > ("MCGeneratorLabel");
  fLarGeantlabel		=pset.get<std::string > ("LarGeantlabel");     
  fUseMCVertex			=pset.get<int > ("UseMCVertex");
}

// ***************** //
cluster::ShowerAngleCluster::~ShowerAngleCluster()
{
}

namespace cluster {
struct SortByWire 
{
  bool operator() (recob::Hit const& h1, recob::Hit const& h2) const 
  { return 
      h1.Wire()->RawDigit()->Channel() < h2.Wire()->RawDigit()->Channel() ;
  }
};
}

// ***************** //
void cluster::ShowerAngleCluster::beginJob()
{

  //temporary:
  unsigned int tpc=0;
  

  /** Get Geometry and detector properties*/
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<util::DetectorProperties> detp;
  art::ServiceHandle<util::LArProperties> larp;
  
  fNPlanes = geo->Nplanes();
  fMean_wire_pitch = geo->WirePitch(0,1,0);    //wire pitch in cm
  ftimetick=detp->SamplingRate()/1000.; 
  std::cout << " In SHowANgle "<<  larp->Efield() << std::endl;
  fdriftvelocity=larp->DriftVelocity(larp->Efield(),larp->Temperature());
  
  /**Get TFileService and define output Histograms*/
  art::ServiceHandle<art::TFileService> tfs;

  
    
  for(unsigned int i=0;i<fNPlanes;++i){
   
    int nwires=geo->Plane(i,tpc).Nwires();
    int ntimes=geo->DetHalfWidth(tpc)*2/(ftimetick*fdriftvelocity);
  

    /**Histos for the angular distribution theta of the shower*/

    fh_theta[i] = tfs->make<TH1F>(Form("fh_theta_%i",i),"Theta distribution",720,-180., 180.);

    tgx[i]=tfs->make<TH2F>(Form("charge distrib_%i",i),"charge distribution per wires",
		 nwires/8.,0,nwires*fMean_wire_pitch,ntimes/8.,0,ntimes*ftimetick*fdriftvelocity);
    
    tgx2[i]=tfs->make<TH2F>(Form("hit distrib_%i",i),"Hit distribution per wires",
			    nwires/8.,0, nwires*fMean_wire_pitch,ntimes/8.,0,ntimes*ftimetick*fdriftvelocity);  

    linefit[i]=tfs->make<TF1>(Form("linefit_%d",i),"pol1",0,4000);
    linefit2[i]=tfs->make<TF1>(Form("linefit_2_%d",i),"pol1",0,4000);	
    
    /**Histos for the angular distribution theta of the shower*/
    fh_omega_evt.push_back( tfs->make<TH1F>(Form("fh_omega_evt_%i",i),
					"Theta distribution per event",720,-180., 180.) );
       
    fh_omega_evt_reb.push_back( tfs->make<TH1F>(Form("fh_omega_evt_reb_%i",i),
					  "Theta distribution per event, rebinned",180,-180., 180.) );
   

    /**Histos for the angular distribution theta wire of the shower*/
   fh_theta_wt[i] = tfs->make<TH1F>(Form("ftheta_wire_%i",i),
				     "Theta wire distribution",720,-180., 180.);
  }  // end loop on planes
  
  ftree_cluster =tfs->make<TTree>("ShowerAngleCluster","Results");/**All-knowing tree with reconstruction information*/
    
    ftree_cluster->Branch("run",&fRun,"run/I");
    ftree_cluster->Branch("subrun",&fSubRun,"subrun/I");
    ftree_cluster->Branch("event",&fEvent,"event/I");
    ftree_cluster->Branch("nplanes",&fNPlanes,"nplanes/I");
  
    ftree_cluster->Branch("mcpdg",&mcpdg,"mcpdg/I");
    ftree_cluster->Branch("mcenergy",&mcenergy,"mcenergy/D");
   
    ftree_cluster->Branch("mcphi",&mcphi,"mcphi/D");
    ftree_cluster->Branch("mctheta",&mctheta,"mctheta/D");
       
    ftree_cluster->Branch("wire_vertex","std::vector<unsigned int>", &fWire_vertex);
    ftree_cluster->Branch("time_vertex","std::vector<double>", &fTime_vertex);

    ftree_cluster->Branch("mcwirevertex","std::vector<unsigned int>", &mcwirevertex);
    ftree_cluster->Branch("mctimevertex","std::vector<double>", &mctimevertex);

      
    ftree_cluster->Branch("wire_last","std::vector<unsigned int>", &fWire_last);
    ftree_cluster->Branch("time_last","std::vector<double>", &fTime_last);

    ftree_cluster->Branch("test_wire_start","std::vector<double>", &test_wire_start);
    ftree_cluster->Branch("test_time_start","std::vector<double>", &test_time_start);

  //  ftree_cluster->Branch("fitw_last","std::vector<double>", &wire_end);
  //  ftree_cluster->Branch("fitt_last","std::vector<double>", &time_end); 

    ftree_cluster->Branch("xyz_vertex","std::vector<double>", &xyz_vertex);
    ftree_cluster->Branch("xyz_vertex_fit","std::vector<double>", &xyz_vertex_fit);

    ftree_cluster->Branch("omega_2d","std::vector<double>", &fOmega_Mean);
    ftree_cluster->Branch("omega_2d_RMS","std::vector<double>", &fOmega_RMS);

    ftree_cluster->Branch("omega_2d_line","std::vector<double>", &fOmega_Mean_line);
    ftree_cluster->Branch("omega_2d_RMS_line","std::vector<double>", &fOmega_RMS_line);

    ftree_cluster->Branch("omega_2d_reb","std::vector<double>", &fOmega_Mean_reb);
    ftree_cluster->Branch("omega_2d_reb_RMS","std::vector<double>", &fOmega_RMS_reb);
    ftree_cluster->Branch("omega_2d_mean","std::vector<double>", &fOmega_Mean_Mean);

    ftree_cluster->Branch("slope","std::vector<double>", &slope);		
    ftree_cluster->Branch("lineslope","std::vector<double>", &lineslope);
    ftree_cluster->Branch("calcslope","std::vector<double>", &calcslope);

    ftree_cluster->Branch("RMS_wire","std::vector<double>", &fRMS_wire);
    ftree_cluster->Branch("RMS_time","std::vector<double>", &fRMS_time);

    ftree_cluster->Branch("Chisq","std::vector<double>", &fChisq);
    ftree_cluster->Branch("minwir","std::vector<double>", &fminwir);
    ftree_cluster->Branch("maxwir","std::vector<double>", &fmaxwir);
    ftree_cluster->Branch("mintime","std::vector<double>", &fmintime);

    ftree_cluster->Branch("maxtime","std::vector<double>", &fmaxtime);
    ftree_cluster->Branch("correlation","std::vector<double>", &fcorrelation);
    ftree_cluster->Branch("covariance","std::vector<double>", &fcovariance);



    ftree_cluster->Branch("Eventangleposition","std::vector<std::vector<double>>", 				&fSingleEvtAngle);
    ftree_cluster->Branch("Eventanglepositionval","std::vector<std::vector<double>>", 				&fSingleEvtAngleVal);

  // ftree_cluster->Branch("fslope_2d"," std::vector<double>", &fSlope_2d);
  // ftree_cluster->Branch("fintercept_2d","std::vector<double>", &fIntercept_2d);
//   
    ftree_cluster->Branch("ShowerPosition2D","std::vector<std::vector<double>>", 				&fShowerPosition2D);
    ftree_cluster->Branch("ShowerWidthProfile2D","std::vector<std::vector<double>>", 				&fShowerWidthProfile2D);
    ftree_cluster->Branch("ShowerChargeProfile2D","std::vector<std::vector<double>>",				&fShowerChargeProfile2D);
  }

// ***************** //
void cluster::ShowerAngleCluster::produce(art::Event& evt)
{ 
  art::ServiceHandle<util::LArProperties> larp;
  mf::LogInfo("ShowerAngleCluster") << " In SHowANgle produce "
				    <<  larp->Efield() << " " 
				    << larp->Efield(1) << " " 
				    << larp->Efield(2);

  /* Get Geometry */
  art::ServiceHandle<geo::Geometry> geo;
  fNPlanes = geo->Nplanes();
  
//     	fWire_vertex.resize(0);  // wire coordinate of vertex for each plane
//     	fTime_vertex.resize(0);  // time coordinate of vertex for each plane
// 	fWire_last.resize(0);  // wire coordinate of vertex for each plane
//     	fTime_last.resize(0);  // time coordinate of vertex for each plane
// 
//         fOmega_Mean.resize(0);    // Mean value of the 2D angular distribution (1=Ind - 0=Coll) cm,cm
//         fOmega_RMS.resize(0);;     // RMS of the 2D angular distribution  (1=Ind - 0=Coll) cm, cm
// 	
//         fOmega_Mean_reb.resize(0);    // Mean value of the 2D angular Rebinned by 4
//         fOmega_RMS_reb.resize(0);     // RMS of the 2D angular distribution  Rebinned by 4
//         fOmega_Mean_Mean.resize(0);    // Mean value of the 2D angular use mean instead of maximum
//         
// 
//         fOmega_wt_Mean.resize(0);; // Mean value of the angular distribution (1=Ind - 0=Coll) wire,time
//         fOmega_wt_RMS.resize(0);;  // RMS of the angular distribution  (1=Ind - 0=Coll) wire,time
//         fChannel_vertex.resize(0);  // wire coordinate of vertex for each plane
//          fChannel_last.resize(0);  // wire coordinate of vertex for each plane



  //fPitch.resize(0);  // Pitch calculated the old way
  fShowerWidthProfile2D.clear(); ;  // vector to show the plane shower Width distribution 
  fShowerChargeProfile2D.clear(); ;  //vector to show the plane shower Charge distribution
  fShowerPosition2D.clear(); ;  //vector to store the positions of hit values stored in the previous two vectors.
  fSingleEvtAngle.clear(); 
  fSingleEvtAngleVal.clear();

 
  fSingleEvtAngle.resize(fNPlanes); 
  fSingleEvtAngleVal.resize(fNPlanes); 
  fShowerWidthProfile2D.resize(fNPlanes); ;  // vector to show the plane shower Width distribution 
  fShowerChargeProfile2D.resize(fNPlanes); ;  //vector to show the plane shower Charge distribution
  fShowerPosition2D.resize(fNPlanes); ;  //vector to store the positions of hit values stored in the previous two vectors.


  for(unsigned int ii=0;ii<fNPlanes;ii++){   
    fSingleEvtAngle[ii].resize(180); 
    fSingleEvtAngleVal[ii].resize(180); 
    fShowerWidthProfile2D[ii].resize(0); ;  // vector to show the plane shower Width distribution 
    fShowerChargeProfile2D[ii].resize(0); ;  //vector to show the plane shower Charge distribution
    fShowerPosition2D[ii].resize(0); ;  //vector to store the positions of hit values stored in the
  }


  // fPitch.resize(fNPlanes); 
 	 
  fWire_vertex.resize(fNPlanes);
  fTime_vertex.resize(fNPlanes);
  mcwirevertex.resize(fNPlanes);  // wire coordinate of vertex for each plane 
  mctimevertex.resize(fNPlanes);  // time coordinate of vertex for each plane
  
  fWire_last.resize(fNPlanes);
  fTime_last.resize(fNPlanes);
  fChannel_vertex.resize(fNPlanes);
  fChannel_last.resize(fNPlanes);

  xyz_vertex.resize(3);
  xyz_vertex_fit.resize(3);	

  test_wire_start.resize(fNPlanes);
  test_time_start.resize(fNPlanes);

  slope.resize(fNPlanes);
  slope_wt.resize(fNPlanes);
  lineslope.resize(fNPlanes);
  calcslope.resize(fNPlanes);

  fOmega_Mean.resize(fNPlanes);    // Mean value of the 2D angular distribution (1=Ind - 0=Coll) cm,cm
  fOmega_RMS.resize(fNPlanes);     // RMS of the 2D angular distribution  (1=Ind - 0=Coll) cm, cm

  fOmega_Mean_line.resize(fNPlanes);    // Mean value of the 2D angular distribution (1=Ind - 0=Coll) cm,cm
  fOmega_RMS_line.resize(fNPlanes);     // RMS of the 2D angular distribution  (1=Ind - 0=Coll) cm, cm

  fOmega_wt_Mean.resize(fNPlanes); // Mean value of the angular distribution (1=Ind - 0=Coll) wire,time
  fOmega_wt_RMS.resize(fNPlanes);  // RMS of the angular distribution  (1=Ind - 0=Coll) wire,time      
  fOmega_Mean_reb.resize(fNPlanes);    // Mean value of the 2D angular Rebinned by 4
  fOmega_RMS_reb.resize(fNPlanes);     // RMS of the 2D angular distribution  Rebinned by 4
  fOmega_Mean_Mean.resize(fNPlanes);    // Mean value of the 2D angular use mean instead of maximum
  fRMS_wire.resize(fNPlanes);
  fRMS_time.resize(fNPlanes);
  fChisq.resize(fNPlanes);
  fminwir.resize(fNPlanes);
  fmaxwir.resize(fNPlanes);
  fmintime.resize(fNPlanes);
  fmaxtime.resize(fNPlanes);
  fcorrelation.resize(fNPlanes);
  fcovariance.resize(fNPlanes);         

  //Find run, subrun and event number:
  fRun = evt.id().run();
  fSubRun = evt.id().subRun();
  fEvent = evt.id().event();

  /**Get Clusters*/
  
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);

  art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, fClusterModuleLabel);

  std::vector< art::PtrVector < recob::Hit> > hitlist_all;
  hitlist_all.resize(fNPlanes);
 
  art::PtrVector<recob::Cluster> clusters;

  // this is temporary until the cluster coming in will represent the actual shower.
  // Currently it sums up showers which are big enough. 
  // Will cause problems with multiple showers. 
  
  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii){

    art::Ptr<recob::Cluster> cl(clusterListHandle, ii);
    std::vector< art::Ptr<recob::Hit> > hitlist = fmh.at(ii);
    unsigned int p(0),w(0), t(0),cs(0); //c=channel, p=plane, w=wire
    GetPlaneAndTPC(hitlist[0],p,cs,t,w);
    
    if(hitlist.size()>15){
      clusters.push_back(cl);

      //loop over cluster hits
      for(art::PtrVector<recob::Hit>::const_iterator a = hitlist.begin(); a != hitlist.end();  a++){ 
	GetPlaneAndTPC(*a,p,cs,t,w);
	hitlist_all[p].push_back(*a);
      }
    }
  } // end temporary loop determining big enough clusters.  
  
  // GetVertex(evt) from MC - should be cut out?;
  /// \todo Never have checks on MC in reconstruction algorithms
  if(fUseMCVertex) GetVertexN(evt);
  
  for(unsigned int i = 0; i < fNPlanes; ++i)
    AngularDistribution(hitlist_all[i]); // 2D Direction of the shower in consecutive planes

  Find2DStartPoints(hitlist_all);
  
  for(unsigned int i=0;i<fNPlanes;i++){
    hitlist_all[i].sort(cluster::SortByWire());
    fh_omega_evt[i]->Reset();
    fh_omega_evt_reb[i]->Reset();
     
    FitAngularDistributions(hitlist_all[i]);  
    Get2DVariables(hitlist_all[i],i);
  }

  

  //create Shower object section:
  
  
 
  // make an art::PtrVector of the clusters
  std::auto_ptr<std::vector<recob::Cluster> > ShowerAngleCluster(new std::vector<recob::Cluster>);
  std::auto_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);

  for(size_t iplane = 0; iplane < fNPlanes; ++iplane){
    
    // figure out the view and total charge for this cluster
    geo::View_t view = hitlist_all[iplane][0]->View();
    double totalQ    = 0.;
    for(size_t ih = 0; ih < hitlist_all[iplane].size(); ++ih)
      totalQ += hitlist_all[iplane][ih]->Charge();

    recob::Cluster temp(fWire_vertex[iplane], fWire_vertex[iplane]*0.05,
			fTime_vertex[iplane], fTime_vertex[iplane]*0.05, 
			fWire_last[iplane], fWire_last[iplane]*0.05,
			fTime_last[iplane], fTime_last[iplane]*0.05, 
			slope[iplane], slope[iplane]*0.05, lineslope[iplane],lineinterc[iplane], 
			totalQ,
			view,
			iplane);


    ShowerAngleCluster->push_back(temp);
    // associate the hits to this cluster
    util::CreateAssn(*this, evt, *(ShowerAngleCluster.get()), hitlist_all[iplane], *(assn.get()));
    mf::LogInfo("ShowerAngleCluster") << "######## in plane loop filling clusters ";
    

    
  }



  /**Fill the output tree with all information */
  ftree_cluster->Fill();

  evt.put(ShowerAngleCluster);
  evt.put(assn);
}


// ******************************* //
int cluster::ShowerAngleCluster::GetPlaneAndTPC(art::Ptr<recob::Hit> a,
						unsigned int &p,
						unsigned int &cs,
						unsigned int &t,
						unsigned int &w)
{
  art::ServiceHandle<geo::Geometry> geo;
  unsigned int c = a->Wire()->RawDigit()->Channel(); 
  geo->ChannelToWire(c,cs,t,p,w);
    
  return 0;
}



//*********************************//
// Angular distribution of the energy of the shower - Collection view
void cluster::ShowerAngleCluster::AngularDistribution(art::PtrVector < recob::Hit>  hitlist){
 
  std::cout << "------ in angular distribution, n of hits " << hitlist.size() << std::endl;
  art::ServiceHandle<geo::Geometry> geom;
  double time;
  unsigned int wire,tpc, cstat;
  unsigned int plane;

  if(hitlist.size()==0)
    return;
  
  art::Ptr<recob::Hit> theHit = (*hitlist.begin());
  time = theHit->PeakTime();  
  GetPlaneAndTPC(hitlist[0],plane,cstat,tpc,wire);
    
  unsigned int minwire=wire,maxwire=0;;
  double mintime=99999,maxtime=0.;

  tgx[plane]->Reset();
  tgx2[plane]->Reset();
  
  // this should changed on the loop on the cluster of the shower
  for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    time = (*hitIter)->PeakTime();  
    //time_C -= (presamplings+10.1);
    GetPlaneAndTPC((*hitIter),plane,cstat,tpc,wire);
    
    maxwire=wire;   

    if(time>maxtime)
	maxtime=time;

    if(time<mintime)
	mintime=time;

  }
 
 // padding of the selected TGraph in wires and time. 
  int wirepad=20;
  int timepad=wirepad*fMean_wire_pitch/(ftimetick*fdriftvelocity)+0.5;
 
  int nbinsx= (maxwire-minwire+2*wirepad)*fMean_wire_pitch;  // nbins to have 
  int nbinsy= (maxtime-mintime+2*timepad)*ftimetick*fdriftvelocity;  // nbins to have 

 
  tgx[plane]->SetBins(nbinsx,(minwire-wirepad)*fMean_wire_pitch,								(maxwire+wirepad)*fMean_wire_pitch,nbinsy,
			  (mintime-timepad)*ftimetick*fdriftvelocity,(maxtime+timepad)*ftimetick*fdriftvelocity);
  
  tgx2[plane]->SetBins(nbinsx,(minwire-wirepad)*fMean_wire_pitch,								(maxwire+wirepad)*fMean_wire_pitch,nbinsy,
			  (mintime-timepad)*ftimetick*fdriftvelocity,(maxtime+timepad)*ftimetick*fdriftvelocity);
 
  
  
  for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != 				hitlist.end();  hitIter++){
    
    time =  (*hitIter)->PeakTime();  
    GetPlaneAndTPC((*hitIter),plane,cstat,tpc,wire);
  
    tgx[plane]->Fill((double)wire*fMean_wire_pitch,
		     time*ftimetick*fdriftvelocity,(*hitIter)->Charge());
    tgx2[plane]->Fill((double)wire*fMean_wire_pitch,time*ftimetick*fdriftvelocity);		
  }
  
 
  tgx[plane]->Fit(Form("linefit_%d",plane),"QMRNCFrob=0.8");
  tgx2[plane]->Fit(Form("linefit_2_%d",plane),"QMRNCFrob=0.95");


  std::cout << "{{{-----}}}  histo stats: rms w,t " << tgx[plane]->GetRMS(1) << " " << 
  tgx[plane]->GetRMS(2) << " chisq " << linefit[plane]->GetChisquare()/linefit[plane]->GetNDF()
  << " max, min wires and times " <<  minwire << " " <<maxwire << " " <<  mintime << " " << maxtime << std::endl;


  fRMS_wire[plane]=tgx[plane]->GetRMS(1);
  fRMS_time[plane]=tgx[plane]->GetRMS(2);
  fChisq[plane]=linefit[plane]->GetChisquare()/linefit[plane]->GetNDF();
  fminwir[plane]=minwire;
  fmaxwir[plane]=maxwire;
  fmintime[plane]=mintime;
  fmaxtime[plane]= maxtime;
  fcorrelation[plane]=tgx[plane]->GetCorrelationFactor();
  fcovariance[plane]=tgx[plane]->GetCovariance();

  return (void)0;
}








// ***************** //
void cluster::ShowerAngleCluster::FitAngularDistributions(art::PtrVector < recob::Hit> hitlist){
  /** Fit function of the angular distribution (cm,cm)*/
 // art::ServiceHandle<geo::Geometry> geo;
  //unsigned int planes = geo->Nplanes();
  //TF1 *gau = new TF1("gaus","gaus",-60, 60);
 std::cout << "------ in angular distribution, n of hits " << hitlist.size() << std::endl;

  art::ServiceHandle<geo::Geometry> geom;
  double time;
  unsigned int wire;
  double BC,AC;
  double omega;
  unsigned int channel,iplane,plane,tpc,cstat;

  if(hitlist.size()==0)
    return;
art::Ptr<recob::Hit> theHit = (*hitlist.begin());
    time = theHit->PeakTime();  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, cstat, tpc, iplane, wire);



   	//tgx[plane]->Set(hitlist.size());

  // this should changed on the loop on the cluster of the shower
  for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    time = theHit->PeakTime();  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, cstat, tpc, plane, wire);
  
      BC = ((double)wire - fWire_vertex[plane])*fMean_wire_pitch; // in cm
      AC = ((double)time - fTime_vertex[plane])*ftimetick*fdriftvelocity; //in cm 
      omega = asin(  AC/sqrt(pow(AC,2)+pow(BC,2)) );
 
      if(BC<0)  // for the time being. Will check if it works for AC<0
	  { 
	  if(AC!=0)
	  omega= AC/fabs(AC)*pi-omega;  //subtract when negative, add when positive
	  else    
	  omega=pi;
	  } 

      omega = 180*omega/3.14;
      fh_theta[plane]->Fill(omega, theHit->Charge()); // Filling the histo (angle, energy of the hit)
      fh_omega_evt[plane]->Fill(omega, theHit->Charge());

  
      fh_omega_evt_reb[plane]->Fill(omega, theHit->Charge());

  }
  //for(unsigned int iplane = 0; iplane < fNPlanes; ++iplane){
  

    fOmega_Mean[iplane] =
    fh_omega_evt[iplane]->GetBinCenter(fh_omega_evt[iplane]->GetMaximumBin());// Mean value of the fit
    fOmega_RMS[iplane] = fh_omega_evt[iplane]->GetRMS(); // RMS of the fit of the angular distribution in deg

    fOmega_Mean_reb[iplane]= fh_omega_evt_reb[iplane]->GetBinCenter(fh_omega_evt_reb[iplane]->GetMaximumBin());// Mean value of the fit
    fOmega_RMS[iplane] = fh_omega_evt_reb[iplane]->GetRMS(); // RMS of the fit of the angular distribution in deg
    fOmega_Mean_Mean[iplane]= fh_omega_evt[iplane]->GetMean();// Mean value of the;    // Mean value of the 2D angular use mean instead of maximum
    
std::cout << "########## intermediate angles, plane: " << iplane << " stand, _w reb, mean " << fOmega_Mean[iplane] << " " << fOmega_Mean_reb[iplane] << " " << fOmega_Mean_Mean[iplane] << std::endl;


fOmega_Mean_line[iplane]=atan(linefit[iplane]->GetParameter(1));


for(int i=0;i<180;i++)
{fSingleEvtAngleVal[iplane][i]=fh_omega_evt_reb[iplane]->GetBinContent(i);
//fSingleEvtAngle[iplane][i]=fh_omega_evt[iplane]->GetBinContent(i);
fSingleEvtAngle[iplane][i]=(double)i*2-180;
}
  //}
//  double  Low_th  = fOmega_Mean[iplane]-(alpha*fOmega_RMS[iplane]);
//  double  High_th = fOmega_Mean[iplane]+(alpha*fOmega_RMS[iplane]);


slope[iplane] = tan((fOmega_Mean[iplane]*pi/180))*fMean_wire_pitch/(ftimetick*fdriftvelocity);
slope_wt[iplane] = tan((fOmega_Mean[iplane]*pi/180));  // ?


calcslope[iplane]=linefit2[iplane]->GetParameter(1);

std::cout << " ((------stand slope and slope from hits only ----- )) " << slope[iplane] << " " << calcslope[iplane] << "  "<< std::endl;

}















void cluster::ShowerAngleCluster::GetVertexN(art::Event& evt){

 art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  evt.getByLabel(fMCGeneratorLabel,mctruthListHandle);




 art::PtrVector<simb::MCTruth> mclist;
  for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);	
      mclist.push_back(mctparticle);
    } 


std::cout << "%%%%%%% mc size size,  "<<mclist.size() <<    std::endl;


    art::Ptr<simb::MCTruth> mc(mclist[0]);
    simb::MCParticle neut(mc->GetParticle(0));

    mcpdg=neut.PdgCode();
    mcenergy=neut.P();  
    
    if (neut.P()){
      double lep_dcosx_truth = neut.Px()/neut.P();
      double lep_dcosy_truth = neut.Py()/neut.P();
      double lep_dcosz_truth = neut.Pz()/neut.P();

     std::cout << "-----  cx,cy,cz " << lep_dcosx_truth << " " << lep_dcosy_truth << " " << lep_dcosz_truth << std::endl;


mcphi=  (lep_dcosx_truth == 0.0 && lep_dcosz_truth == 0.0) ? 0.0 : TMath::ATan2(lep_dcosx_truth,lep_dcosz_truth);
mctheta= (lep_dcosx_truth == 0.0 && lep_dcosy_truth == 0.0 && lep_dcosz_truth == 0.0) ? 0.0 : TMath::Pi()*0.5-TMath::ATan2(sqrt(lep_dcosx_truth*lep_dcosx_truth + lep_dcosz_truth*lep_dcosz_truth),lep_dcosy_truth);




mcphi=180*mcphi/TMath::Pi();
mctheta= 180*mctheta/TMath::Pi();
    std::cout << "-----  phi, theta " <<  mcphi << " " << mctheta << std::endl;

    }

    
    
    
    
// if((neut.PdgCode()==22)&& neut.StatusCode()==1) //photon - need first electron.
//      { 
//      art::Handle< std::vector<sim::Particle> > parHandle;
//      evt.getByLabel(fLarGeantlabel, parHandle);
// 
//      art::PtrVector<simb::MCParticle> pvec;
//     int fpart=0;
//     for(unsigned int i = 0; i < parHandle->size(); ++i){
//       art::Ptr<simb::MCParticle> p(parHandle, i);      
//       pvec.push_back(p);
//       if(p->PdgCode() ==11 || p->PdgCode()==-11)
// 	  {
// 		fpart=i;
// 		break;
// 	  }	
// 
//     }
// 
// 
//      std::cout << "%%%&&&&&&&&&& is PDG: " << pvec[fpart]->PdgCode() << " " << pvec[fpart]->TrackId() << std::endl;
//       
// 
//     int trackid ;
// 
// 	trackid =  mc->GetParticle(0).Daughter(0);
//         std::cout << "####### NDaughters: " << trackid << " "<<mc->GetParticle(0).NumberDaughters() <<  " "<<mc->GetParticle(1).NumberDaughters() <<  std::endl;
// 
// 	for(int xx=0;xx<mc->GetParticle(0).NumberDaughters();xx++)
// 		{
//       trackid =  mc->GetParticle(0).Daughter(xx);
//         std::cout << "####### is PDG, trackid: " << trackid << " "<<mc->GetParticle(0).NumberDaughters() << std::endl; 
//               } 
// 	unsigned int jj;
//       for(jj = 0; jj < pvec.size(); jj++) // Don't look below i.
// 	    {
// 	      if (trackid==pvec[jj]->TrackId())
// 		{
// 		   std::cout << "daughter particle "<<jj << " " << pvec[jj]->PdgCode() << std::endl; // get the pointer, 
// 			break; 
//             
//                  }            
// 
// 
//     	   }
//      neut=*pvec[fpart];
// 
//      } //end foton clause
//if((neut.PdgCode()==11 || neut.PdgCode()==-11 )&& neut.StatusCode()==1)
//NumberDaughters	(		 ) 	 const [inline]
//int simb::MCParticle::Daughter	(	const int 	i	 ) 	 const
 
 //  std::cout << "%%%%%%% particle size, partslist: "<< partslist << " "  <<  mc->NParticles() << std::endl; 
  int npart=0;
   //  while(&& npart < mc->NParticles() )
     //     {
                 std::cout << "%%%%%%%####### is PDG: "<< npart <<" " << neut.PdgCode() << std::endl; 
 	//	neut=mc->GetParticle(npart++);

       //   }       

 std::cout << "%%%%%%%####### after loop is PDG: "<< npart <<" " << neut.PdgCode() << std::endl; 
    //if((neut.PdgCode()==11 || neut.PdgCode()==-11 )&& neut.StatusCode()==1){
  
    
    xyz_vertex[0] =neut.Vx();
    xyz_vertex[1] =neut.Vy();
    xyz_vertex[2] =neut.Vz();
	
    std::cout<<"neut.Vx()= "<<neut.Vx()<<" ,y= "<<neut.Vy()<<" ,z= "<<neut.Vz()<<std::endl;
//if(((neut.PdgCode()==11 || neut.PdgCode()==-11 )&& neut.StatusCode()==1))
  //    break;

    

	
    
   art::ServiceHandle<geo::Geometry> geom;
   art::ServiceHandle<util::LArProperties> larp;
   

double drifttick=(xyz_vertex[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./ftimetick);

const double origin[3] = {0.};
for(unsigned int iplane=0;iplane<fNPlanes;iplane++)
{
double pos[3];
 unsigned int  wirevertex, t, cs;
unsigned int p;
geom->Plane(iplane).LocalToWorld(origin, pos);
	//planex[p] = pos[0];
std::cout << "plane X positionp " << iplane << " " << pos[0] << std::endl;

pos[1]=xyz_vertex[1];
pos[2]=xyz_vertex[2];
 unsigned int channel2 = geom->NearestChannel(pos,iplane);
 geom->ChannelToWire(channel2,cs,t,p,wirevertex); 
       
if(iplane!=p)
	{std::cout << " error - planes don't match " << iplane << " " << p << std::endl;
	return ;
	}

fWire_vertex[p]=wirevertex;
fTime_vertex[p]=drifttick-(pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./ftimetick)+60;
std::cout<<"wirevertex= "<<wirevertex<< " timevertex " << fTime_vertex[p] << " correction "<< (pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./ftimetick) << " " << pos[0] <<std::endl;
 mcwirevertex[p]=wirevertex;  // wire coordinate of vertex for each plane
 mctimevertex[p]=drifttick-(pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./ftimetick);  // time coordinate of vertex for each plane


}



  return (void)0;
}




//int cluster::ShowerAngleCluster::Get2Dvariables(float Wire_vertexI_wt, float Wire_vertexC_wt, float Time_I_wt, float Time_C_wt){
void cluster::ShowerAngleCluster::Get2DVariables(art::PtrVector < recob::Hit> hitlist,unsigned int iplane) {  

    
//get parameters of the slope obtained by searching for the maximum and start points
double tst=fTime_vertex[iplane];
double wst= fWire_vertex[iplane];
double slp=slope[iplane];

double wireend=fWire_last[iplane];
//double timeend=time_end[iplane];

double intercept=tst-slp*(double)wst;
art::ServiceHandle<geo::Geometry> geom;

//get slope of lines orthogonal to those found crossing the shower.
 double aprim=0;
      
	if(slp)	
	{
	aprim=-1./slp;
	}

  std::cout << "========= line params, inside 2d variables, plane: a,c " << iplane <<" " << slp << " " << intercept << std::endl;

double slp_cm=slope[iplane]/fMean_wire_pitch*(ftimetick*fdriftvelocity);
double projlength=(wireend-wst)*fMean_wire_pitch*TMath::Sqrt(1+slp_cm*slp_cm);



int nbins = (hitlist.size() > 1000) ? 100 : hitlist.size()/10;

std::cout << "---- projected length for plane: " << iplane << " " << projlength << " nbins " << nbins <<  std::endl;

TH1F * hithist=new TH1F(Form("hithist_ev_%d_pl_%d",fEvent,iplane),Form("hithist_ev_%d_pl_%d",fEvent,iplane),nbins,0,projlength);  
  
TH1F * hithist2=new TH1F(Form("hithist_ev_%d_pl_%d_w",fEvent,iplane),Form("hithist_ev_%d_pl_%d",fEvent,iplane),nbins,0,projlength);  

// start loop to calculate the profile of the shower along the shower axis.

//double extreme_intercept_end=-999999;
//double extreme_intercept_start=999999;

//double extr_wire_pos=0,extr_time_pos=0;

// int multiplier=1;   // +1 for positive angles, -1 for negative angles. to compensate that we are looking for either the highest (omega >0 ) or lowest (omega<0) intercept.
// 
// if(slp>0)
// 	{
// 	multiplier=1;
// 	}
// else if(slp<0)
// 	{
// 	multiplier=-1;
// 	}
  
  

for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    unsigned int wire,cstat, tpc,channel,plane;
    channel = theHit->Wire()->RawDigit()->Channel();
    geom->ChannelToWire(channel, cstat, tpc, plane, wire);
 
    if(iplane!=plane)
      continue;
    
    double ort_intercept=time-aprim*(double)wire;
    double wire_on_line=(ort_intercept - intercept)/(slp-aprim); 
    double time_on_line=slp*wire_on_line+intercept;     
    
   // std::cout << "plane: wire on line, time on line " << plane << " : " << wire_on_line << " " << time_on_line << " for wire and time " << wire << " " << time << " pos: " << TMath::Sqrt( pow((wire_on_line-wst)*fMean_wire_pitch,2)+pow((time_on_line-tst)*(ftimetick*fdriftvelocity),2) ) << " prof: " << TMath::Sqrt( pow((wire_on_line-(double)wire)*fMean_wire_pitch,2)+pow((time_on_line-time)*(ftimetick*fdriftvelocity),2) )  << std::endl;
  
    double linedist=TMath::Sqrt( pow((wire_on_line-wst)*fMean_wire_pitch,2)+pow((time_on_line-tst)*(ftimetick*fdriftvelocity),2));
    double ortdist=TMath::Sqrt( pow((wire_on_line-(double)wire)*fMean_wire_pitch,2)+pow((time_on_line-time)*(ftimetick*fdriftvelocity),2) );
    
    hithist->Fill(linedist);
    hithist2->Fill(linedist,ortdist);
    
    fShowerPosition2D[plane].push_back(linedist);  
    fShowerWidthProfile2D[plane].push_back(ortdist);
    fShowerChargeProfile2D[plane].push_back(theHit->Charge()); 
  }
 

 TFile * trh=new TFile("histos.root","UPDATE");
 
hithist->Write();
hithist2->Write();
 
 delete trh;
 delete hithist;
 

  return (void)0;
}






void   cluster::ShowerAngleCluster::Find2DStartPoints(std::vector< art::PtrVector < recob::Hit> > hitlist_all)
{


  //// find for which planes the correlation factor is largest (and preferably larger than 0.6)
  std::vector<double>  wire_start,wire_end;
  std::vector<double>  time_start,time_end;
 
  wire_start.resize(fNPlanes);wire_end.resize(fNPlanes);
  time_start.resize(fNPlanes);time_end.resize(fNPlanes);
  
  std::vector< int > best_planes;

  Find2DBestPlanes( best_planes);
 
//// find the hits near the beginning for the best planes.

  art::ServiceHandle<geo::Geometry> geom;
  double time;
  unsigned int wire,plane,tpc,cstat;

  double a,c;
 // double wlend,tlend;

 // double wire_bar=0,time_bar=0;
 // int nhits=0;
 
  // loop on selected planes
  for (unsigned int ii=0;ii<best_planes.size();ii++)
  {

  unsigned int iplane=best_planes[ii];

  //get first wire of plane (should be sorted by wire number)
  if(hitlist_all[iplane].size()==0) // this should never happen.
    continue;
  
  GetPlaneAndTPC((*hitlist_all[iplane].begin()),plane,cstat,tpc,wire);
  
  //error checking:
  if(iplane!=plane){
    std::cout << " error: planes mismatch  " << iplane << " "<< plane << std::endl;
    return;
    }

  //get paramters of the straight line fit. (and rescale them to cm/cm)
  a=linefit[iplane]->GetParameter(1);
  lineslope[iplane]=a*fMean_wire_pitch/(ftimetick*fdriftvelocity);
  c=linefit[iplane]->GetParameter(0);
  lineinterc[iplane]=c/(ftimetick*fdriftvelocity);

  //get slope of lines orthogonal to those found crossing the shower.
  double aprim=0;
  if(a){
    aprim=-1./a;
  }

  std::cout << "========= line params, plane: a,c " << plane << " " << a << " " << slope[iplane] << " " << c << std::endl;

  // find extreme intercepts. For the time being we don't care which one is the start point and which one is the endpoint.
  
  double extreme_intercept_high=-999999;
  double extreme_intercept_low=999999;
 Find_Extreme_Intercepts(hitlist_all[iplane],aprim,extreme_intercept_high,extreme_intercept_low);
  
  
  double extreme_intercept_start=-999999;
  double extreme_intercept_end=999999;
  
 // int multiplier=1;   // +1 for positive angles, -1 for negative angles. to compensate that we are looking for either the highest (omega >0 ) or lowest (omega<0) intercept.
  
  if(a>0) {  // for the time being assuming forward going showers
   //  multiplier=1;
     extreme_intercept_end=extreme_intercept_high;
     extreme_intercept_start=extreme_intercept_low;
  }
  else if(a<0){
   // multiplier=-1;
    extreme_intercept_start=extreme_intercept_high;
    extreme_intercept_end=extreme_intercept_low;
  }
  

  //   wlend=(extreme_intercept_end - c)/(a-aprim);   //in cm
  //   tlend=(a*wlend+c);   // in cm

 

  std::cout << "^^^^^^^^^ a^prim + max and min intercept " << aprim << " " << extreme_intercept_end << " " << extreme_intercept_start << std::endl;

  
  int wire_online_end=(extreme_intercept_end - c)/(a-aprim); 
  int wire_online_begin=(extreme_intercept_start - c)/(a-aprim); 
  double time_online_end=a*wire_online_end+c; 
  double time_online_begin=a*wire_online_begin+c;


  std::cout << " :::::::: wire_online begin point " << wire_online_begin << " " << time_online_begin << std::endl;

  //calculate the first and last cluster points on the through line:
  
  art::Ptr<recob::Hit> startHit=FindClosestHit(hitlist_all[iplane], wire_online_begin,time_online_begin);
  art::Ptr<recob::Hit> endHit=FindClosestHit(hitlist_all[iplane], wire_online_end,time_online_end);
  
  GetPlaneAndTPC(startHit,plane,cstat,tpc,wire);
  wire_start[plane]=wire;
  time_start[plane]=startHit->PeakTime();
  
  GetPlaneAndTPC(endHit,plane,cstat,tpc,wire);
  wire_end[plane]=wire;
  time_end[plane]=endHit->PeakTime();
  
  //   wlend=wire_end[plane]*fMean_wire_pitch;   // temporary - will need to get last wire coordinate for each plane and hitlist.
  //   tlend=a*wlend+c;  // in cm
  // 
  //   time_end[plane]=tlend/(ftimetick*fdriftvelocity);

  } // end of  loop on selected planes



  // temporary cross-check  
  for(unsigned int ii=0;ii<best_planes.size();ii++)
	std::cout << " ----++++ determined start wire points per planes " << best_planes[ii] << " "  << wire_start[best_planes[ii]] << " " << time_start[best_planes[ii]] <<std::endl;  

  //sort the best planes in order of correlation factor
  //first sort
  if(fabs(tgx[best_planes[0]]->GetCorrelationFactor())<fabs(tgx[best_planes[1]]->GetCorrelationFactor()))
    std::swap(best_planes[0],best_planes[1]);

  //second sort
  if(best_planes.size()>2 && 
    (fabs(tgx[best_planes[1]]->GetCorrelationFactor())<fabs(tgx[best_planes[2]]->GetCorrelationFactor())))
	std::swap(best_planes[1],best_planes[2]);

  // sanity check - see if times of the points found are close enough.
	/////////////////
  for(unsigned int ii=0;ii<fNPlanes;ii++)	
	std::cout << " ---- plane parameters " << ii << " " << tgx[ii]->GetCorrelationFactor() << " RMS ratio " << tgx[ii]->GetRMS(1) /tgx[ii]->GetRMS(2) << " Chisq " << linefit[ii]->GetChisquare()/linefit[ii]->GetNDF() << " covariance " << tgx[ii]->GetCovariance() << std::endl;
	
	
  const double origin[3] = {0.};
  std::vector <std::vector <  double > > position;

  // get starting positions for all planes
  for(unsigned int xx=0;xx<fNPlanes;xx++){
    double pos1[3];
    geom->Plane(xx).LocalToWorld(origin, pos1);
    std::vector <double > pos2;
    pos2.push_back(pos1[0]);
    pos2.push_back(pos1[1]);
    pos2.push_back(pos1[2]);
    position.push_back(pos2);
  }

  ///////////////// Clean up done up to here.
	
  //loop to check time discrepancies between points found
  for(unsigned int xx=0;xx<best_planes.size()-1;xx++){  // best_planes.size()-1 because we want to always find a next plane
    for(unsigned int yy=xx+1;yy<best_planes.size();yy++){
	std::cout << "**** difference between planes X position, planes: "<< best_planes[xx] << " "<< best_planes[yy] << " "<< position[best_planes[xx]][0] << " " << position[best_planes[yy]][0] << " " << fabs(position[best_planes[xx]][0]-position[best_planes[yy]][0]) << " " << " " << time_start[best_planes[xx]] << " " << time_start[best_planes[yy]] << " " << fabs(time_start[best_planes[xx]]-time_start[best_planes[yy]]) << std::endl;
			
	if(fabs(time_start[best_planes[xx]]-time_start[best_planes[yy]]) > 1.5*fabs(position[best_planes[xx]][0]-position[best_planes[yy]][0])){
		std::cout << " !!!! Something's wrong in the wire determination " << fabs(time_start[best_planes[xx]]-time_start[best_planes[yy]]) << " " << 1.5*fabs(position[xx][0]-position[yy][0]) << std::endl;
	}	
    } // end inner yy loop
  } // end outer xx loop

	
	

	if((fabs(time_start[best_planes[0]]-time_start[best_planes[1]]) > 1.5*fabs(position[0][0]-position[1][0])) && best_planes.size() > 2) //time discrepancy of best correlation factor pair and we have a choice.
		{std::cout << " time discrepancy of best correlation factor pair 0 and 1 " << time_start[best_planes[0]] << " " << time_start[best_planes[1]] << " "<< position[0][0] << " " << position[1][0] << std::endl;
		if(!(fabs(time_start[best_planes[0]]-time_start[best_planes[2]]) > 2.5*fabs(position[0][0]-position[2][0]))) //0,1 is bad but 0,2 is ok:
			{
			std::cout << " but 0,2 is ok " << time_start[best_planes[0]] << " " << time_start[best_planes[2]] << " "<< position[0][0] << " " << position[2][0] << std::endl;
			std::swap(best_planes[1],best_planes[2]);
			}
		else   //0,1 is not ok and 0,2 is not ok
			{
			std::cout << " 0,1 and 0,2 is not ok " << std::endl;
			if(!(fabs(time_start[best_planes[1]]-time_start[best_planes[2]]) > 2.5*fabs(position[1][0]-position[2][0]))) //0,1 and 0,2 is bad but 1,2 is ok.
				{
				std::cout << " but 1,2 is ok " << time_start[best_planes[1]] << " " << time_start[best_planes[2]] << " "<< position[1][0] << " " << position[2][0] << std::endl;
				std::swap(best_planes[0],best_planes[1]);
				std::swap(best_planes[1],best_planes[2]);  // shift zero to last position.
				}
			}

		}


	for(unsigned int ii=0;ii<best_planes.size();ii++)
	std::cout << " ----++++ determined start wire points per planes " << best_planes[ii] << " "  << wire_start[best_planes[ii]] << " " << time_start[best_planes[ii]] <<std::endl;  


	// Assuming there is no problem ( and we found the best pair that comes close in time )
	// we try to get the Y and Z coordinates for the start of the shower. 
	int chan1=geom->PlaneWireToChannel(best_planes[0],wire_start[best_planes[0]], 0);
	int chan2=geom->PlaneWireToChannel(best_planes[1],wire_start[best_planes[1]], 0);

	double y,z;
	bool wires_cross = geom->ChannelsIntersect(chan1,chan2,y,z);

	
	xyz_vertex_fit[1]=y;
	xyz_vertex_fit[2]=z;
	xyz_vertex_fit[0]=time_start[best_planes[0]]*fdriftvelocity*ftimetick+position[0][0];


	std::cout << ":::::: found x,y,z vertex " << wires_cross << " " << xyz_vertex_fit[0] << " " << y << " " << z << std::endl;


	// assume some condition - probably that correlation factor is too small. Then project the found vertex in x,y,z into wire, time coordinates in the last plane.
///////////////// Only do the following part if there are 3 planes.

	if(fNPlanes>=3){
	
	double pos[3];
	unsigned int  wirevertex, t,cstat;
	unsigned int worst_plane=2;
	if(best_planes.size()>=3)
		worst_plane=best_planes[2];
	else  //find plane that has bad correlation factor. We know, that at least two are there.
		{
		std::cout << "bplane size <3 " << best_planes.size() << std::endl;
		for(unsigned int jj=0;jj<fNPlanes;jj++)
			{
			bool exist_flag=false;		
			for(unsigned int kk=0;kk<best_planes.size();kk++)			
				{if(jj==(unsigned int)best_planes[kk])
					exist_flag=true;
				std::cout << " jj,kk, true or false " << jj << " " << kk << std::endl;	
				}

			if(!exist_flag)  // adding non_existing flag
				{
				worst_plane=jj;
				std::cout << "setting worst plane to " << jj << std::endl;
				break;
				}
			}
		}	

	geom->Plane(worst_plane).LocalToWorld(origin, pos);
	//planex[p] = pos[0];
	std::cout << "plane X positionp " << worst_plane << " " << pos[0] << std::endl;

	
	
	///////////////////////////////////////
	// geometry test:
// 	double width  = 2.*geom->TPC(0).HalfWidth();  //notice the geometry gives the 1/2 width, so multiply by 2
// 	double height = 2.*geom->TPC(0).HalfHeight(); //notice the geometry gives the 1/2 height, so multiply by 2
// 	double length =    geom->TPC(0).Length();     //notice the geometry gives the total length
// 	
// 	
// 	std::cout << "-------- height " << height << " " << length << " " <<width << std::endl;
// 	
// 	//display first and last wires:
// 	
// 	for(unsigned int iplane=0;iplane<fNPlanes;iplane++)
// 	    { std::cout << " +++ pl" << iplane << geom->Plane(iplane).Nwires()  << std::endl;
// 	      double wire1_Start[3]={0},wire1_End[3]={0};
// 	    geom->WireEndPoints(0,iplane,0,wire1_Start,wire1_End);
// 	      std::cout << "++++++++ wire positions: w nr " << 0 << " " << wire1_Start[0] << " " << wire1_Start[1] << " " << wire1_Start[2] << " "<< wire1_End[0] << " " << wire1_End[1] << " "<< wire1_End[2] << " " << std::endl;
// 	
// 	      geom->WireEndPoints(0,iplane,geom->Plane(iplane).Nwires()-1,wire1_Start,wire1_End);
// 	      std::cout << "++++++++ wire positions: w nr " << geom->Plane(iplane).Nwires()-1 << " " << wire1_Start[0] << " " << wire1_Start[1] << " " << wire1_Start[2] << " "<< wire1_End[0] << " " << wire1_End[1] << " "<< wire1_End[2] << " " << std::endl;
// 	
// 	      
// 	
// 	    }
	
	// display all wire positions
	
// 	for(unsigned int iplane=0;iplane<fNPlanes;iplane++)
// 	{
// 	  std::cout << " +++ pl" << iplane << geom->Plane(iplane).Nwires()  << std::endl;
// 	for(unsigned int wire=0;wire<geom->Plane(iplane).Nwires();wire+=10)
// 	{
// 	  double wire1_Start[3]={0},wire1_End[3]={0};
// 	 geom->WireEndPoints(0,iplane,wire,wire1_Start,wire1_End);
// 	//std::cout << "++++++++ wire positions: w nr " << wire << " " << wire1_Start[0] << " " << wire1_Start[1] << " " << wire1_Start[2] << " "<< wire1_End[0] << " " << wire1_End[1] << " "<< wire1_End[2] << " " << std::endl;
// 	 
// 	}
// 	std::cout << std::endl << std::endl << std::endl;
// 	}
// 	
// 	
// 	
// 	for(double iy=-60;iy<60;iy+=15)
// 	{
// 	  // for(double iz=length-500;iz<length-100;iz+=10)
// 	  for(double iz=1000;iz<1100;iz+=50)
// 	  {
// 	   unsigned int wirev[3]; 
// 	     std::cout << " ---- plane for y,z: " << iy << " " << iz <<" ";
// 	  for(unsigned int iplane=0;iplane<fNPlanes;iplane++)
// 	  {unsigned int p;
// 	  geom->Plane(iplane).LocalToWorld(origin, pos);
// 	  pos[1]=iy;
// 	  pos[2]=iz;
// 	    unsigned int channel2 = geom->NearestChannel(pos);
// 	  geom->ChannelToWire(channel2,t,p,wirev[iplane]); 
// 	 std::cout << " p: " << iplane << " " << wirev[iplane];
// 	  } // end of wire finding for positions
// 	   std::cout << std::endl;  
// 	  
// 	   // second loop to test intersect:
// 	   
// 	    std::cout << " ---- y,z for planes : " << std::endl; 
// 	   for(unsigned int iplane=0;iplane<fNPlanes-1;iplane++)
// 	      {
// 	      for(unsigned int iplane2=iplane+1;iplane2<fNPlanes;iplane2++)
// 		{
// 		int chan1=geom->PlaneWireToChannel(iplane,wirev[iplane], 0);
// 		int chan2=geom->PlaneWireToChannel(iplane2,wirev[iplane2], 0);
// 
// 		double y,z;
// 		bool wires_cross = geom->ChannelsIntersect(chan1,chan2,y,z);
// 	        
// 	       std::cout << iplane << " " << iplane2 << " " << y << " " << z << std::endl;
// 	      
// 	      wires_cross = geom->ChannelsIntersect(chan2,chan1,y,z);
// 	        
// 	       std::cout << " inverse " << iplane2 << " " << iplane << " " << y << " " << z << std::endl;
// 	      
// 	      
// 		}
// 	      }
// 	  
// 	 }
// 	}
	 
	//// end geometry test 
	/////////////////////////////////////////////////////////////////////////////////////////////// 
	 
	pos[1]=xyz_vertex_fit[1];
	pos[2]=xyz_vertex_fit[2];
 	unsigned int channel2 = geom->NearestChannel(pos,worst_plane);
       	geom->ChannelToWire(channel2,cstat,t,worst_plane,wirevertex); 


	art::ServiceHandle<util::LArProperties> larp;
	art::ServiceHandle<util::DetectorProperties> detp;
	double drifttick=(xyz_vertex_fit[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./ftimetick);


	double timestart=drifttick-(pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./ftimetick);//+detp->TriggerOffset();
	std::cout << " worst plane " << worst_plane <<" wirevertex= "<<wirevertex<< " timevertex " << timestart << " correction " << detp->TriggerOffset() << " " << (pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./ftimetick) << " "<<pos[0] <<std::endl;


	double min_dist=999999.;

	for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist_all[worst_plane].begin(); hitIter != hitlist_all[worst_plane].end();  hitIter++){
    		art::Ptr<recob::Hit> theHit = (*hitIter);
    		time = theHit->PeakTime() ;  
    		unsigned int plane;
    		GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
		
	
    		double dist_begin=TMath::Sqrt( pow((double)((int)wirevertex-(int)wire)*fMean_wire_pitch,2)+pow((timestart-time)*fdriftvelocity*ftimetick,2) );	

		//std::cout << "=== min_dist " << wire << " " << time <<" " << dist_begin << " " << pow((double)((int)wirevertex-(int)wire)*fMean_wire_pitch,2) << " " << ((int)wirevertex-(int)wire)*fMean_wire_pitch << " " << min_dist << std::endl; 
		
		if(dist_begin<min_dist)
			{
			min_dist=dist_begin;
			wire_start[worst_plane]=wire;
			time_start[worst_plane]=time;

			}


	} // end loop on hits.


} // end big if(fNPlanes >= 3)



std::cout << " final wire, time vertices for all planes:  " << std::endl;  

  for(unsigned int ii=0;ii<fNPlanes;ii++){
    std::cout << ii << " " << wire_start[ii] << " " << time_start[ii] << " " <<                geom->Plane(ii,0).SignalType() << " " << geo::kCollection << " " << geo::kInduction  << std::endl;
    fWire_vertex[ii]=wire_start[ii];
    fTime_vertex[ii]=time_start[ii];
    fWire_last[ii]=wire_end[ii];
    fTime_last[ii]=time_end[ii];
   }
  }

  
  
  
// ------------------ select the clusters for the best two planes. This should make sure it's looking at one TPC in the future.


void cluster::ShowerAngleCluster::Find2DBestPlanes(std::vector<int> &best_planes){
   
   
  // first loop: select planes where the shower has a large enough correlation factor
  for (unsigned int iplane=0;iplane<fNPlanes;iplane++){
    if(fabs(tgx[iplane]->GetCorrelationFactor()) > 0.6 )	
      best_planes.push_back(iplane);
    std::cout << " correlation factors in 0.6 search" << tgx[iplane]->GetCorrelationFactor() << std::endl;
  }

 
  //second loop: select planes where the shower in the plane is near horizontal (based on RMS ratio)
  for (unsigned int iplane=0;iplane<fNPlanes;iplane++){
    if(fabs(tgx[iplane]->GetCorrelationFactor()) <= 0.6 
      && fabs(tgx[iplane]->GetCorrelationFactor()) > 0.5 
      && tgx[iplane]->GetRMS(2) > 0. 
      && tgx[iplane]->GetRMS(1) /tgx[iplane]->GetRMS(2)>2. )	
	  best_planes.push_back(iplane);
    }
 
 
  //// Find which plane has the highest correlation factor, and if there is only one in best_planes, add the one  with next largest correl. factor
  unsigned int used_plane=999;
  while(best_planes.size()<2)
    {
    double mincorr=0;
    int maxplane=0;
    for (unsigned int iplane=0;iplane<fNPlanes;iplane++){
	if(fabs(tgx[iplane]->GetCorrelationFactor()) <= 0.6 
	  && fabs(tgx[iplane]->GetCorrelationFactor()) > mincorr 
	  && iplane!=used_plane ){
	    maxplane=iplane;
	    mincorr=fabs(tgx[iplane]->GetCorrelationFactor());
	    }
	} // end for loop on planes

    std::cout << "pushing back " << maxplane << std::endl;
    used_plane=maxplane;  // to cut out redundancy.
    best_planes.push_back(maxplane);

    } // end while loop

  /////test values:
  for(unsigned int ii=0;ii<best_planes.size();ii++)
    std::cout << "======+++==== " << best_planes[ii] << " " << 
    tgx[best_planes[ii]]->GetCorrelationFactor() << std::endl; 


//   for (unsigned int iplane=0;iplane<fNPlanes;iplane++){
//     std::cout << "other parameters: " <<  std::endl;
//     std::cout << "RMS:: " <<  tgx[iplane]->GetRMS(1) << " " << tgx[iplane]->GetRMS(2) << " " << tgx[iplane]->GetRMS(1) /tgx[iplane]->GetRMS(2) << std::endl;
//     
//     TH1D *PX=tgx[iplane]->ProjectionX();
//     std::cout << " nbins " << PX->GetNbinsX() << std:: endl;
//     
//     
//     // this somehow was calculated earlier? in max an minwire?
//     int maxix=-1,minix=-1;
//     for(int i=0;i<PX->GetNbinsX();i++){
//       if(PX->GetBinContent(i)>0){
// 	minix=i;
// 	break;
//       }
//     }
//     for(int i=PX->GetNbinsX();i>0;i--){
//       if(PX->GetBinContent(i)>0){
// 	maxix=i;
// 	break;
//       }
//     }
//      
//         
//     TH1D *PY=tgx[iplane]->ProjectionY();
//        
//     std::cout << " nbins " << PY->GetNbinsX() << std:: endl;
//     
//     int maxiy=-1,miniy=-1;
//     for(int i=0;i<PY->GetNbinsX();i++){
//       if(PY->GetBinContent(i)>0){
// 	miniy=i;
// 	break;
//       }
//     }
//     for(int i=PY->GetNbinsX();i>0;i--){
//       if(PY->GetBinContent(i)>0){
// 	maxiy=i;
// 	break;
//       }
//     }
//     std::cout << "max,min x,y:: " << PX->GetBinCenter(minix) << " "<< PX->GetBinCenter(maxix) << " " << PX->GetBinCenter(maxix) - PX->GetBinCenter(minix) << " y: " << PY->GetBinCenter(miniy) << " "<< PY->GetBinCenter(maxiy) << " " << PY->GetBinCenter(maxiy) - PY->GetBinCenter(miniy) << " "<<(PX->GetBinCenter(maxix) - PX->GetBinCenter(minix))/(PY->GetBinCenter(maxiy) - PY->GetBinCenter(miniy)) << std::endl;
//     
//     std::cout << " chisq " << linefit[iplane]->GetChisquare()<< " " << linefit[iplane]->GetNDF() << " "<< linefit[iplane]->GetChisquare()/linefit[iplane]->GetNDF() << std::endl;
//  
//   }

   
  return;
}




void cluster::ShowerAngleCluster::Find_Extreme_Intercepts(art::PtrVector<recob::Hit> hitlist,
			     double perpslope,
			     double &inter_high,
			     double &inter_low)  
{
  
  inter_high=-999999;
  inter_low=999999;

  unsigned int plane,tpc,wire,cstat;
  
  for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
    
    //wire_bar+=wire;
    //time_bar+=time;	
    //nhits++;
    
    double intercept=time*ftimetick*fdriftvelocity-perpslope*(double)wire*fMean_wire_pitch;
    
    if(intercept > inter_high ){
      inter_high=intercept;
    }
    if(intercept < inter_low ){
      inter_low=intercept;
    }  

    

  }   // end of first HitIter loop, at this point we should have the extreme intercepts 

}




art::Ptr<recob::Hit> cluster::ShowerAngleCluster::FindClosestHit(art::PtrVector<recob::Hit> hitlist,
			     unsigned int wire_online,
			     double time_online)
{
  
  double min_length_from_start=99999;
  art::Ptr<recob::Hit> nearHit;
   
  unsigned int plane,tpc,wire,cstat;
   
   
  for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
    
    double dist_mod=TMath::Sqrt( pow(((double)wire_online-(double)wire*fMean_wire_pitch),2)+pow((time_online-time*fdriftvelocity*ftimetick),2) );	

    if(dist_mod<min_length_from_start){
	//wire_start[plane]=wire;
	//time_start[plane]=time;
	nearHit=(*hitIter);
	min_length_from_start=dist_mod;
	}	

  } 
  
return nearHit;    
}



