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


#include "SimulationBase/simbase.h"
#include "Simulation/SimListUtils.h"
#include "RawData/RawDigit.h"
#include "Utilities/LArProperties.h"
#include "SummaryData/summary.h"


// ***************** //


cluster::ShowerAngleCluster::ShowerAngleCluster(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
  produces< std::vector<recob::Cluster> >();
}


void cluster::ShowerAngleCluster::reconfigure(fhicl::ParameterSet const& pset) 
{
  fClusterModuleLabel = pset.get< std::string >("ClusterModuleLabel");
   //fHoughLineModuleLabel=pset.get<std::string > ("HoughLineModuleLabel");
  fVertexCLusterModuleLabel=pset.get<std::string > ("VertexClusterModuleLabel");
  fMCGeneratorLabel=pset.get<std::string > ("MCGeneratorLabel");
  fLarGeantlabel=pset.get<std::string > ("LarGeantlabel");     
  fUseMCVertex=pset.get<int > ("UseMCVertex");
 

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
      h1.Wire()->RawDigit()->Channel() < 
      h2.Wire()->RawDigit()->Channel() ;
  }
};
}

// ***************** //
void cluster::ShowerAngleCluster::beginJob()
{



  /** Get Geometry*/
  art::ServiceHandle<geo::Geometry> geo;
  fNPlanes = geo->Nplanes();
  fMean_wire_pitch = geo->WirePitch(0,1,0);    //wire pitch in cm
  unsigned int tpc=0;


  /**Get TFileService and define output Histograms*/
  art::ServiceHandle<art::TFileService> tfs;

  /** Create Histos names*/
 
  char tit_h_theta[128] = {0};
  char tit_h_theta_wt[128] = {0};
  
  
 
  
  
  
  
  for(unsigned int i=0;i<fNPlanes;++i){


    //    sprintf(&tit_dedx[0],"fh_dedx_%.4i_%.4i_%i",i);
   
   int nwires=geo->Plane(i,tpc).Nwires();
   int ntimes=geo->DetHalfWidth(tpc)*2/(ftimetick*0.158);
//ftimetick      =  0.198; // time sample in us fdriftvelocity =  0.157;	
std::cout << "{{{ --- }}} init for plane" << i << " " << nwires << " " << ntimes << std::endl; 


    /**Histos for the angular distribution theta of the shower*/
    sprintf(&tit_h_theta[0],"fh_theta_%i",i);
    fh_theta[i] = tfs->make<TH1F>(tit_h_theta,"Theta distribution",720,-180., 180.);

    sprintf(&tit_h_theta[0],"charge distrib_%i",i);
    tgx[i]=tfs->make<TH2F>(tit_h_theta,"charge distribution per wires",nwires/8.,0, nwires,ntimes/80.,0,ntimes/2);
   sprintf(&tit_h_theta[0],"hit distrib_%i",i);
    tgx2[i]=tfs->make<TH2F>(tit_h_theta,"Hit distribution per wires",nwires/8.,0, nwires,ntimes/80.,0,ntimes/2);  


    linefit[i]=tfs->make<TF1>(Form("linefit_%d",i),"pol1",0,4000);
    linefit2[i]=tfs->make<TF1>(Form("linefit_2_%d",i),"pol1",0,4000);	
    
    /**Histos for the angular distribution theta of the shower*/
    sprintf(&tit_h_theta[0],"fh_omega_evt_%i",i);
    fh_omega_evt.push_back( tfs->make<TH1F>(tit_h_theta,"Theta distribution per event",720,-180., 180.) );

   
   sprintf(&tit_h_theta[0],"fh_omega_evt_reb_%i",i);
    fh_omega_evt_reb.push_back( tfs->make<TH1F>(tit_h_theta,"Theta distribution per event, rebinned",180,-180., 180.) );
   

   //fh_omega2_evt.push_back( new TH1F(tit_h_theta,"Theta distribution per event",720,-180., 180.) );    


    /**Histos for the angular distribution theta wire of the shower*/
    sprintf(&tit_h_theta[0],"ftheta_wire_%i",i);
    fh_theta_wt[i] = tfs->make<TH1F>(tit_h_theta_wt,"Theta wire distribution",720,-180., 180.);
 	}
  
  ftree_cluster =tfs->make<TTree>("ShowerAngleCluster","Results");/**All-knowing tree with reconstruction information*/
  
  
   ftree_cluster->Branch("run",&fRun,"run/I");
    ftree_cluster->Branch("subrun",&fSubRun,"subrun/I");
   ftree_cluster->Branch("event",&fEvent,"event/I");
   ftree_cluster->Branch("nplanes",&fNPlanes,"nplanes/I");
  
   ftree_cluster->Branch("wire_vertex","std::vector<unsigned int>", &fWire_vertex);
   ftree_cluster->Branch("time_vertex","std::vector<double>", &fTime_vertex);
   ftree_cluster->Branch("wire_last","std::vector<unsigned int>", &fWire_last);
   ftree_cluster->Branch("time_last","std::vector<double>", &fTime_last);


  ftree_cluster->Branch("fitw_vertex","std::vector<double>", &wire_start);
   ftree_cluster->Branch("fitt_vertex","std::vector<double>", &time_start);
   ftree_cluster->Branch("fitw_last","std::vector<double>", &wire_end);
   ftree_cluster->Branch("fitt_last","std::vector<double>", &time_end); 

 //  ftree_cluster->Branch("Pitch","std::vector<double>", &fPitch);

// this should be temporary - until the omega is sorted out.
  // ftree_cluster->Branch("fh_omega2_evt","std::vector<TH1F *>", &fh_omega2_evt);


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



   ftree_cluster->Branch("Eventangleposition","std::vector<std::vector<double>>",&fSingleEvtAngle);
ftree_cluster->Branch("Eventanglepositionval","std::vector<std::vector<double>>",&fSingleEvtAngleVal);

  // ftree_cluster->Branch("fslope_2d"," std::vector<double>", &fSlope_2d);
  // ftree_cluster->Branch("fintercept_2d","std::vector<double>", &fIntercept_2d);
//   
ftree_cluster->Branch("ShowerPosition2D","std::vector<std::vector<double>>",&fShowerPosition2D);
ftree_cluster->Branch("ShowerWidthProfile2D","std::vector<std::vector<double>>",&fShowerWidthProfile2D);
ftree_cluster->Branch("ShowerChargeProfile2D","std::vector<std::vector<double>>",&fShowerChargeProfile2D);
  


}

// ***************** //
void cluster::ShowerAngleCluster::produce(art::Event& evt)
{ 

  /**TODO THIS VALUES SHOULD BE PARAMETERS OF THE MODULE, evtl from database
  *double Efield_SI        =  0.7;     // Electric Field between Shield and Induction planes in kV/cm
  *double Efield_IC        =  0.9;     // Electric Field between Induction and Collection planes in kV/cm
  *double Temperature      = 87.6;  // LAr Temperature in K
  *check if there can be a replacement later for the product, not needed now
  *double timepitch        = fdriftvelocity*ftimetick;                  //time sample (cm) */

  /** Get Geometry */
  art::ServiceHandle<geo::Geometry> geo;
  //art::ServiceHandle<util::LArProperties> larprop;
  fNPlanes = geo->Nplanes();
  //fdriftvelocity = larprob->DriftVelocity(Efield_SI,Temperature);
  
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
    	

 
fSingleEvtAngle.resize(fNPlanes); 
fSingleEvtAngleVal.resize(fNPlanes); 

 fShowerWidthProfile2D.resize(fNPlanes); ;  // vector to show the plane shower Width distribution 
 fShowerChargeProfile2D.resize(fNPlanes); ;  //vector to show the plane shower Charge distribution
 fShowerPosition2D.resize(fNPlanes); ;  //vector to store the positions of hit values stored in the previous two vectors.


 for(unsigned int ii=0;ii<fNPlanes;ii++)
  {   
  fSingleEvtAngle[ii].resize(180); 
  fSingleEvtAngleVal[ii].resize(180); 
  }


      // fPitch.resize(fNPlanes); 
    	 
	fWire_vertex.resize(fNPlanes);
	fTime_vertex.resize(fNPlanes);
        fWire_last.resize(fNPlanes);
	fTime_last.resize(fNPlanes);
 	fChannel_vertex.resize(fNPlanes);
        fChannel_last.resize(fNPlanes);

   wire_start.resize(fNPlanes);wire_end.resize(fNPlanes);;
   time_start.resize(fNPlanes);time_end.resize(fNPlanes);;

	slope.resize(fNPlanes);
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




  /**Get Clusters*/
  std::cout << "************ What I'm getting out " << fClusterModuleLabel << " " << std::endl;
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);


   std::vector< art::PtrVector < recob::Hit> > hitlist_all;
   //art::PtrVector < recob::Hit> hitlistInd;
   hitlist_all.resize(fNPlanes);

 
    //std::auto_ptr<std::vector<recob::Shower> > Shower3DVector(new std::vector<recob::Shower>);

  for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {

      art::Ptr<recob::Cluster> cl(clusterListHandle, ii);
      

      art::PtrVector<recob::Hit> hitlist;
      hitlist = cl->Hits();
      hitlist.sort(cluster::SortByWire());
      unsigned int p(0),w(0), c(0), t(0); //c=channel, p=plane, w=wire

      art::PtrVector<recob::Hit>::const_iterator a = hitlist.begin();
      c=(*a)->Wire()->RawDigit()->Channel(); 
      geo->ChannelToWire(c,t,p,w);
      wire_start[p]=w;

//       a = hitlist.end();
//       c=(*a)->Wire()->RawDigit()->Channel(); 
//       geo->ChannelToWire(c,t,p,w);
//       wire_end[p]=w;	


      for(art::PtrVector<recob::Hit>::const_iterator a = hitlist.begin(); a != hitlist.end();  a++) //loop over cluster hits
      {
	c=(*a)->Wire()->RawDigit()->Channel(); 
	geo->ChannelToWire(c,t,p,w);

        //geo->Plane(i).View()

	//std::cout << "+++++++ planeview in hit list   " <<  geo->Plane(p,t).View() << std::endl; 

          hitlist_all[geo->Plane(p,t).View()-1].push_back(*a);

// 	if(geo->Plane(p,t).SignalType() == geo::kCollection)
// 	  {	  
// 	    hitlistCol.push_back(*a);
// 
// 	  }
// 	else if (geo->Plane(p,t).SignalType() == geo::kInduction)
// 	  {
// 	    hitlistInd.push_back(*a);
// 	  } 
      }

    //  hitlistInd.sort(cluster::SortByWire());
     // hitlistCol.sort(cluster::SortByWire());


    } // End loop on clusters.

 fRun = evt.id().run();
 fSubRun = evt.id().subRun();
 fEvent = evt.id().event();



   // GetVertex(evt);
if(fUseMCVertex)
    GetVertexN(evt);



   for(unsigned int i=0;i<fNPlanes;i++)
      {
       hitlist_all[i].sort(cluster::SortByWire());
      fh_omega_evt[i]->Reset();
      fh_omega_evt_reb[i]->Reset();
      AngularDistribution(hitlist_all[i]); // 2D Direction of the shower in consecutive planes
      FitAngularDistributions(i);  
      Get2DVariables(hitlist_all[i]);
       }
 
 
  //AngularDistribution(hitlistCol); // 2D Direction of the shower Collection
              // Fit of 2d distributions (for both planes)
  
  std::cout << "######## in main loop " << fOmega_Mean[0] << " " <<  fOmega_Mean[1] << std::endl;
  
art::PtrVector<recob::Hit> clusterHits;



// make an art::PtrVector of the clusters
std::auto_ptr<std::vector<recob::Cluster> > ShowerAngleCluster(new std::vector<recob::Cluster>);
unsigned int i=0;

for(unsigned int iplane=0;iplane<fNPlanes;iplane++)
{
 clusterHits.clear();
	int maxlength=0,maxpos=0;
 	for( i = 0; i < clusterListHandle->size(); ++i){
   		art::Ptr<recob::Cluster> prod(clusterListHandle, i);

		if((prod->View()-1)!=iplane)
			continue;

	//	std::cout << "----- looping on clusters " << i << " " << prod->View()-1 << " " <<prod->Hits().size() << std::endl;
	//	std::cout << "----- lenght and pos " << i << " " << maxlength << " "  <<maxpos << std::endl;

   		if(  (prod->Hits().size() > (unsigned int)maxlength) )
		{
        		maxpos=i;
        		maxlength=prod->Hits().size();
          //        std::cout << "++++ condition ok " << std::endl << std::endl;
		}
 	}

art::Ptr<recob::Cluster> prod(clusterListHandle, maxpos);

for(unsigned int ii=0;ii<prod->Hits().size();ii++)
  clusterHits.push_back(prod->Hits()[ii]);


int wirevert=prod->StartPos()[0];
double timevert=prod->StartPos()[1];

int wireend=wire_end[iplane];
double timeend=time_end[iplane];


if(fUseMCVertex)
{
wirevert=fWire_vertex[iplane];
timevert=fTime_vertex[iplane];
}

recob::Cluster temp(clusterHits,wirevert, prod->SigmaStartPos()[0],timevert, prod->SigmaStartPos()[1], wireend, prod->SigmaEndPos()[0],timeend, prod->SigmaEndPos()[1], slope[iplane], slope[iplane]*0.05, lineslope[iplane],lineinterc[iplane], iplane);

std::cout << "######## in plane loop filling clusters " << std::endl; 

ShowerAngleCluster->push_back(temp);
}



  /**Fill the output tree with all information */
    ftree_cluster->Fill();


  evt.put(ShowerAngleCluster);

}


// ******************************* //

// Angular distribution of the energy of the shower - Collection view
void cluster::ShowerAngleCluster::AngularDistribution(art::PtrVector < recob::Hit>  hitlist){
  std::cout << "------ in angular distribution, n of hits " << hitlist.size() << std::endl;
  int    loop = 0; // flag
  art::ServiceHandle<geo::Geometry> geom;
  double time;
  unsigned int wire;
  double BC,AC;
  double omega;
  unsigned int channel,plane;


 art::Ptr<recob::Hit> theHit = (*hitlist.begin());
    time = theHit->PeakTime();  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, tpc, plane, wire);
  
unsigned int minwire=wire,maxwire=0;;
double mintime=99999,maxtime=0.;

	tgx[plane]->Reset();
        tgx2[plane]->Reset();
   	//tgx[plane]->Set(hitlist.size());

  // this should changed on the loop on the cluster of the shower
  for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    time = theHit->PeakTime();  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, tpc, plane, wire);
    
   maxwire=wire;   

   if(time>maxtime)
	maxtime=time;

if(time<mintime)
	mintime=time;


   tgx[plane]->Fill((double)wire,time,theHit->Charge());
   tgx2[plane]->Fill((double)wire,time);		

    BC = ((double)wire - fWire_vertex[plane])*fMean_wire_pitch; // in cm
    AC = ((double)time - fTime_vertex[plane])*ftimetick*0.158; //in cm 
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
    loop++; // flag counter
  }
  std::cout << "VertexWireC= " << fWire_vertex[plane] << "   VerTimeC= " << fTime_vertex[plane] << std::endl;

double w_bar=0,t_bar=0;
int nhits=0;

for(int i=0;i<tgx2[plane]->GetNbinsX();i++)
	for(int j=0;j<tgx2[plane]->GetNbinsY();j++)
		{
		if(tgx2[plane]->GetBinContent(i,j)<=3)
			tgx2[plane]->SetBinContent(i,j,0);
		else
		      {
			w_bar+=tgx2[plane]->GetXaxis()->GetBinCenter(i);
			t_bar+=tgx2[plane]->GetYaxis()->GetBinCenter(j);
			nhits++;
			}	
		}

std::cout << "/////// w_bar, t_bar, nhits " << w_bar << " " << t_bar << " " << nhits << " " << w_bar/nhits << " " << t_bar/nhits <<std::endl;  

w_bar/=nhits;
t_bar/=nhits;

double Sxx=0,Syy=0,Sxy=0;

for(int i=0;i<tgx2[plane]->GetNbinsX();i++)
	for(int j=0;j<tgx2[plane]->GetNbinsY();j++)
		{
		if(tgx2[plane]->GetBinContent(i,j)>3)
			{Sxx+=(tgx2[plane]->GetXaxis()->GetBinCenter(i)-w_bar)*(tgx2[plane]->GetXaxis()->GetBinCenter(i)-w_bar);
			Syy+=(tgx2[plane]->GetYaxis()->GetBinCenter(j)-t_bar)*(tgx2[plane]->GetYaxis()->GetBinCenter(j)-t_bar);
			Sxy+=(tgx2[plane]->GetYaxis()->GetBinCenter(j)-t_bar)*(tgx2[plane]->GetXaxis()->GetBinCenter(i)-w_bar);;
			}
		}


std::cout << "/////// Sxx,Syy,Sxy " << Sxx << " " << Syy << " " << Sxy << "newest slope: " << Sxy/Sxx <<std::endl;  


 
  tgx[plane]->Fit(Form("linefit_%d",plane),"QMRNCFrob=0.8");
  tgx2[plane]->Fit(Form("linefit_2_%d",plane),"QMRNCFrob=0.95");


std::cout << "{{{-----}}}  histo stats: rms w,t " << tgx[plane]->GetRMS(1) << " " << tgx[plane]->GetRMS(2) << " chisq " << linefit[plane]->GetChisquare()/linefit[plane]->GetNDF() << " max, min wires and times " <<
minwire << " " <<maxwire << " " <<  mintime << " " << maxtime << std::endl;


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
void cluster::ShowerAngleCluster::FitAngularDistributions(int iplane){
  /** Fit function of the angular distribution (cm,cm)*/
  art::ServiceHandle<geo::Geometry> geo;
  //unsigned int planes = geo->Nplanes();
  //TF1 *gau = new TF1("gaus","gaus",-60, 60);



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

if((neut.PdgCode()==22)&& neut.StatusCode()==1) //photon - need first electron.
     { 
     art::Handle< std::vector<sim::Particle> > parHandle;
     evt.getByLabel(fLarGeantlabel, parHandle);

     art::PtrVector<simb::MCParticle> pvec;
    int fpart=0;
    for(unsigned int i = 0; i < parHandle->size(); ++i){
      art::Ptr<simb::MCParticle> p(parHandle, i);      
      pvec.push_back(p);
      if(p->PdgCode() ==11 || p->PdgCode()==-11)
	  {
		fpart=i;
		break;
	  }	

    }


     std::cout << "%%%&&&&&&&&&& is PDG: " << pvec[fpart]->PdgCode() << " " << pvec[fpart]->TrackId() << std::endl;
      

    int trackid ;

	trackid =  mc->GetParticle(0).Daughter(0);
        std::cout << "####### NDaughters: " << trackid << " "<<mc->GetParticle(0).NumberDaughters() <<  " "<<mc->GetParticle(1).NumberDaughters() <<  std::endl;

	for(int xx=0;xx<mc->GetParticle(0).NumberDaughters();xx++)
		{
      trackid =  mc->GetParticle(0).Daughter(xx);
        std::cout << "####### is PDG, trackid: " << trackid << " "<<mc->GetParticle(0).NumberDaughters() << std::endl; 
              } 
	unsigned int jj;
      for(jj = 0; jj < pvec.size(); jj++) // Don't look below i.
	    {
	      if (trackid==pvec[jj]->TrackId())
		{
		   std::cout << "daughter particle "<<jj << " " << pvec[jj]->PdgCode() << std::endl; // get the pointer, 
			break; 
            
                 }            


    	   }
     neut=*pvec[fpart];

     } //end foton clause
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
double drifttick=(xyz_vertex[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./.198);

const double origin[3] = {0.};
for(unsigned int p=0;p<fNPlanes;p++)
{
double pos[3];
unsigned int  wirevertex, t;
geom->Plane(p).LocalToWorld(origin, pos);
	//planex[p] = pos[0];
std::cout << "plane X positionp " << p << " " << pos[0] << std::endl;

pos[1]=xyz_vertex[1];
pos[2]=xyz_vertex[2];
 unsigned int channel2 = geom->NearestChannel(pos);
       geom->ChannelToWire(channel2,t,p,wirevertex); 
       
fWire_vertex[p]=wirevertex;
fTime_vertex[p]=drifttick+(pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./.198);
std::cout<<"wirevertex= "<<wirevertex<< " timevertex " << fTime_vertex[p] <<std::endl;
}



  return (void)0;
}




//int cluster::ShowerAngleCluster::Get2Dvariables(float Wire_vertexI_wt, float Wire_vertexC_wt, float Time_I_wt, float Time_C_wt){
void cluster::ShowerAngleCluster::Get2DVariables(art::PtrVector < recob::Hit> hitlist) {  

  art::ServiceHandle<geo::Geometry> geom;
  // only needed for drawing the axis of the shower in the event display
  
  unsigned int channel;
 // double omega_sh, wire_cm, time_cm;


// double AC, BC, omega; 
  double time;
  unsigned int wire,plane;

double a,c;
double wlst,wlend;
double tlst,tlend;

double wire_bar=0,time_bar=0;
int nhits=0;

//double minlength={10000};

///////// enter hits? or just wire coordinates?! 


art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin();
    channel = (*hitIter)->Wire()->RawDigit()->Channel();
 geom->ChannelToWire(channel, tpc, plane, wire);

a=linefit[plane]->GetParameter(1);
lineslope[plane]=a;
 c=linefit[plane]->GetParameter(0);
lineinterc[plane]=c;
 wlst=fWire_vertex[plane];   // temporary - will need to get first wire coordinate for each plane and hitlist.
 tlst=a*wlst+c;

 double aprim=0;
      if(a)	
	{aprim=-1./a;
	}

// second loop, where the last points are found...
// for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
//     art::Ptr<recob::Hit> theHit = (*hitIter);
//     time = theHit->PeakTime();  
//     //time_C -= (presamplings+10.1);
//     art::Ptr<recob::Wire> theWire = theHit->Wire();
//     channel = theWire->RawDigit()->Channel();
//     geom->ChannelToWire(channel, tpc, plane, wire);
// 
// 
// 
//     }


  std::cout << "========= line params, plane: a,c " << plane << " " << a << " " << slope[plane] << " " << c << std::endl;


////////////// temporary - assuming that line is in direction of shower - need to be smarter in the future. 


double extreme_intercept_end=0;
double extreme_intercept_start=999999;
int multiplier=1;   // +1 for positive angles, -1 for negative angles. to compensate that we are looking for either the highest (omega >0 ) or lowest (omega<0) intercept.

if(fOmega_Mean[plane]>0)
	{//extreme_intercept=-99999;
	multiplier=1;
	}
else if(fOmega_Mean[plane]<0)
	{//extreme_intercept=99999;
	multiplier=-1;
	}

for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    time = theHit->PeakTime() ;  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, tpc, plane, wire);
 
   wire_bar+=wire;
   time_bar+=time;	
   nhits++;
    
      if(a)	
	{
// if omega> 0 looking for larger intercept. if omega<0 looking for smaller intercept as endpoint.
       

	double intercept=time-aprim*(double)wire;
 //	double wire_on_line=(intercept - c)/(a-aprim);
   //     double time_on_line=a*wire_on_line+c;  

      //std::cout << "[[----]] intercepts check, wire,time " << wire<< " " <<time << " int " << intercept << " end:  " << multiplier*extreme_intercept_end << " start: " << multiplier*extreme_intercept_start << std::endl;

	if(fabs(intercept) > multiplier*extreme_intercept_end ) 	
		{
		extreme_intercept_end=intercept;
		//wire_end[plane]=wire;
		
		}

	if(fabs(intercept) < multiplier*extreme_intercept_start ) 	
		{
		extreme_intercept_start=intercept;
		//wire_start[plane]=wire;
		
		}


       



	}




}   // end of HitIter loop


std::cout << ":::::::::: mean wire and time " << wire_bar << " " << time_bar << " " <<nhits << " " << wire_bar/nhits << " " << time_bar/nhits << std::endl; 

wire_bar/=nhits;
time_bar/=nhits;

wlst=(extreme_intercept_start - c)/(a-aprim);
tlst=a*wlst+c;

wlend=(extreme_intercept_end - c)/(a-aprim);
tlend=a*wlend+c;

double Sxx=0,Sxy=0,Syy=0;

std::cout << "^^^^^^^^^ a^prim + max and min intercept " << aprim << " " << extreme_intercept_end << " " << extreme_intercept_start << std::endl;

double min_length_from_start=99999;//,min_length_from_end=99999;
//int wire_online_end=(extreme_intercept_end+multiplier*100 - c)/(a-aprim);
int wire_online_begin=(extreme_intercept_start-multiplier*100 - c)/(a-aprim);
//double time_online_end=a*wire_online_end+c;
double time_online_begin=a*wire_online_begin+c;

//std::cout << "^^^^^^^^^ wire_time_online_begin and end " <<  wire_online_begin << " " << time_online_begin << " " << wire_online_end << " " << time_online_end << std::endl;

//third loop to find first and last points (in theory)
for(art::PtrVector<recob::Hit>::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    time = theHit->PeakTime() ;  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, tpc, plane, wire);

	double intercept=time-aprim*(double)wire;
 	double wire_on_line=(intercept - c)/(a-aprim);
        double time_on_line=a*wire_on_line+c;  

//	double dist=TMath::Sqrt( pow(wire_on_line-(double)wire,2)+pow(time_on_line-time,2) );
//	double dist_from_start=TMath::Sqrt( pow(wire_on_line-wlst,2)+pow(time_on_line-tlst,2) );


    double dist_begin=TMath::Sqrt( pow(wire_online_begin-(int)wire,2)+pow(time_online_begin-time,2) );

    Sxx+=(wire-wire_bar)*(wire-wire_bar);
    Syy+=(time-time_bar)*(time-time_bar);
    Sxy+=(time-time_bar)*(wire-wire_bar);

    if(dist_begin<min_length_from_start)
	{
	wire_start[plane]=wire;
	time_start[plane]=time;
	min_length_from_start=dist_begin;
	}	

 //  double dist_end=TMath::Sqrt( pow(wire_online_end-(int)wire,2)+pow(time_online_end-time,2) );

    //if(dist_end<min_length_from_end)
	//{
	//wire_end[plane]=wire;
	//min_length_from_end=dist_end;
	//}


//std::cout << "((----- )) point at " << wire << " " << time << "dist_b,_e : " << dist_begin << " " << dist_end << " " <<   min_length_from_start << " " << min_length_from_end  << " ot dist "<< dist << " " << dist_from_start <<std::endl;


	fShowerPosition2D[plane].push_back(TMath::Sqrt( pow(wire_on_line-wlst,2)+pow(time_on_line-tlst,2) ));  
	fShowerWidthProfile2D[plane].push_back(TMath::Sqrt( pow(wire_on_line-(double)wire,2)+pow(time_on_line-time,2) ));
	fShowerChargeProfile2D[plane].push_back(theHit->Charge()); 
	



}


 wlst=wire_start[plane];
 tlst=a*wlst+c;
 //time_start[plane]=tlst; 

 wire_end[plane]=wlend;   // temporary - will need to get last wire coordinate for each plane and hitlist.
 tlend=a*wlend+c;

time_end[plane]=tlend;

//std::cout << "======== start and end positions for plane" << plane << " " << wlst << " " << tlst << " " << wlend << " " << tlend << std::endl; 



std::cout << ":::::::: Sxx, Syy, Sxy " << Sxx << " " << Syy << " " << Sxy << " slope " << Sxy/Sxx << " corr vars " << Sxx/nhits<< " " << Syy/nhits <<" " << Sxy/nhits<< std::endl;


//ok. now have, presumable start and endpoint positions for cluster. 




  return (void)0;
}

