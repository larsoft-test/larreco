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
#include "art/Framework/Core/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Core/TFileDirectory.h" 
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

  /**Get TFileService and define output Histograms*/
  art::ServiceHandle<art::TFileService> tfs;

  /** Create Histos names*/
 
  char tit_h_theta[128] = {0};
  char tit_h_theta_wt[128] = {0};
  char sh_long_tit[128] = {0};
  char sh_tit[128] = {0};
  char shT_tit[128] = {0};
  
  int nbins;
  
  
  
  
  for(unsigned int i=0;i<fNPlanes;++i){


    //    sprintf(&tit_dedx[0],"fh_dedx_%.4i_%.4i_%i",i);
   

    /**Histos for the angular distribution theta of the shower*/
    sprintf(&tit_h_theta[0],"fh_theta_%i",i);
    fh_theta[i] = tfs->make<TH1F>(tit_h_theta,"Theta distribution",720,-180., 180.);

    
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


 //  ftree_cluster->Branch("Pitch","std::vector<double>", &fPitch);

// this should be temporary - until the omega is sorted out.
  // ftree_cluster->Branch("fh_omega2_evt","std::vector<TH1F *>", &fh_omega2_evt);


    ftree_cluster->Branch("omega_2d","std::vector<double>", &fOmega_Mean);
    ftree_cluster->Branch("omega_2d_RMS","std::vector<double>", &fOmega_RMS);

    ftree_cluster->Branch("omega_2d_reb","std::vector<double>", &fOmega_Mean_reb);
    ftree_cluster->Branch("omega_2d_reb_RMS","std::vector<double>", &fOmega_RMS_reb);
    ftree_cluster->Branch("omega_2d_mean","std::vector<double>", &fOmega_Mean_Mean);
   


   ftree_cluster->Branch("Eventangleposition","std::vector<std::vector<double>>",&fSingleEvtAngle);
ftree_cluster->Branch("Eventanglepositionval","std::vector<std::vector<double>>",&fSingleEvtAngleVal);

  // ftree_cluster->Branch("fslope_2d"," std::vector<double>", &fSlope_2d);
  // ftree_cluster->Branch("fintercept_2d","std::vector<double>", &fIntercept_2d);
//   

  



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
  
    	fWire_vertex.resize(0);  // wire coordinate of vertex for each plane
    	fTime_vertex.resize(0);  // time coordinate of vertex for each plane
	fWire_last.resize(0);  // wire coordinate of vertex for each plane
    	fTime_last.resize(0);  // time coordinate of vertex for each plane

        fOmega_Mean.resize(0);    // Mean value of the 2D angular distribution (1=Ind - 0=Coll) cm,cm
        fOmega_RMS.resize(0);;     // RMS of the 2D angular distribution  (1=Ind - 0=Coll) cm, cm
	
        fOmega_Mean_reb.resize(0);    // Mean value of the 2D angular Rebinned by 4
        fOmega_RMS_reb.resize(0);     // RMS of the 2D angular distribution  Rebinned by 4
        fOmega_Mean_Mean.resize(0);    // Mean value of the 2D angular use mean instead of maximum
        

        fOmega_wt_Mean.resize(0);; // Mean value of the angular distribution (1=Ind - 0=Coll) wire,time
        fOmega_wt_RMS.resize(0);;  // RMS of the angular distribution  (1=Ind - 0=Coll) wire,time
        fChannel_vertex.resize(0);  // wire coordinate of vertex for each plane
         fChannel_last.resize(0);  // wire coordinate of vertex for each plane



	//fPitch.resize(0);  // Pitch calculated the old way
    	

 
fSingleEvtAngle.resize(fNPlanes); 
fSingleEvtAngleVal.resize(fNPlanes); 

 for(int ii=0;ii<fNPlanes;ii++)
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


        fOmega_Mean.resize(fNPlanes);    // Mean value of the 2D angular distribution (1=Ind - 0=Coll) cm,cm
        fOmega_RMS.resize(fNPlanes);     // RMS of the 2D angular distribution  (1=Ind - 0=Coll) cm, cm
        fOmega_wt_Mean.resize(fNPlanes); // Mean value of the angular distribution (1=Ind - 0=Coll) wire,time
        fOmega_wt_RMS.resize(fNPlanes);  // RMS of the angular distribution  (1=Ind - 0=Coll) wire,time
        fOmega_Mean_reb.resize(fNPlanes);    // Mean value of the 2D angular Rebinned by 4
        fOmega_RMS_reb.resize(fNPlanes);     // RMS of the 2D angular distribution  Rebinned by 4
        fOmega_Mean_Mean.resize(fNPlanes);    // Mean value of the 2D angular use mean instead of maximum
      


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

      for(art::PtrVectorItr<recob::Hit> a = hitlist.begin(); a != hitlist.end();  a++) //loop over cluster hits
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














   for(int i=0;i<fNPlanes;i++)
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

for(int iplane=0;iplane<fNPlanes;iplane++)
{
 clusterHits.clear();
	int maxlength=0,maxpos=0;
 	for( i = 0; i < clusterListHandle->size(); ++i){
   		art::Ptr<recob::Cluster> prod(clusterListHandle, i);

		if((prod->View()-1)!=iplane)
			continue;

	//	std::cout << "----- looping on clusters " << i << " " << prod->View()-1 << " " <<prod->Hits().size() << std::endl;
	//	std::cout << "----- lenght and pos " << i << " " << maxlength << " "  <<maxpos << std::endl;

   		if(  (prod->Hits().size() > maxlength) )
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
if(fUseMCVertex)
{
wirevert=fWire_vertex[iplane];
timevert=fTime_vertex[iplane];
}

recob::Cluster temp(clusterHits,wirevert, prod->SigmaStartPos()[0],timevert, prod->SigmaStartPos()[1], prod->EndPos()[0], prod->SigmaEndPos()[0],prod->EndPos()[1], prod->SigmaEndPos()[1], slope[iplane], slope[iplane]*0.05, -999.,-999., iplane);

std::cout << "######## in plane loop filling clusters " << std::endl; 

ShowerAngleCluster->push_back(temp);
}













//for(unsigned int tt=0;tt<prodvec.size();tt++)
//ShowerAngleCluster->push_back(prodvec[tt]);

//get direction cosines and set them for the shower




//
  /**Fill the output tree with all information */
    ftree_cluster->Fill();

//for(unsigned int iplane = 0; iplane < fNPlanes; ++iplane)
  //fh_theta[iplane]->Write(Form("fh_theta_%d_%d",iplane,evt.id().event()));
  // This needs work, clearly.  
  //for(int p=0;p<2;p++)Shower3DVector->push_back(shower);
  evt.put(ShowerAngleCluster);

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



for(int i=0;i<180;i++)
{fSingleEvtAngleVal[iplane][i]=fh_omega_evt_reb[iplane]->GetBinContent(i);
//fSingleEvtAngle[iplane][i]=fh_omega_evt[iplane]->GetBinContent(i);
fSingleEvtAngle[iplane][i]=(double)i*2-180;
}
  //}
  double  Low_th  = fOmega_Mean[iplane]-(alpha*fOmega_RMS[iplane]);
  double  High_th = fOmega_Mean[iplane]+(alpha*fOmega_RMS[iplane]);


slope[iplane] = tan((fOmega_Mean[iplane]*pi/180))*fMean_wire_pitch/(ftimetick*fdriftvelocity);

}

// ***************** //


// ***************** //
//int cluster::ShowerAngleCluster::Get2Dvariables(float Wire_vertexI_wt, float Wire_vertexC_wt, float Time_I_wt, float Time_C_wt){
void cluster::ShowerAngleCluster::Get2DVariables(art::PtrVector < recob::Hit> hitlist) {  

  art::ServiceHandle<geo::Geometry> geom;
  // only needed for drawing the axis of the shower in the event display
  
  unsigned int channel;
  double omega_sh, wire_cm, time_cm;


 double AC, BC, omega; 
  double time;
  unsigned int wire,plane;

for(art::PtrVectorItr<recob::Hit> hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    time = theHit->PeakTime() ;  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, tpc, plane, wire);
    //    if(time_C<1020)continue;
    wire_cm = wire * fMean_wire_pitch; //in cm
    time_cm = time *ftimetick*fdriftvelocity; //in cm
   
 /*  if(hitIter == hitlist.begin())
	{ fWire_vertex[plane] = fWire_vertex[plane] ; //in cm
          fTime_vertex[plane] = fTime_vertex[plane] ; //in cmm        
        }  */  

   
    // moving to polar coordinates
    BC = (wire_cm - fWire_vertex[plane]* fMean_wire_pitch)+fMean_wire_pitch; //in cm
    AC = (time_cm - fTime_vertex[plane]*ftimetick*fdriftvelocity); // in cm 
    omega = asin(AC/sqrt(pow(AC,2)+pow(BC,2)));
    omega = 180*omega/3.14; // in deg
    //std::cout << " WireI1=" << wireI1 << " BI= " << BI << "    ThetaI = " << thetaI <<std::endl;
       
    if( (omega>(fOmega_Mean[plane]-1.0))&&(omega<(fOmega_Mean[plane]+1.0)) ){
      fWire_last[plane] = wire;
      fChannel_last[plane]=  channel;  // wire coordinate of vertex for each plane
      fTime_last[plane] = time;
    }

}   // end of HitIter loop

 
  return (void)0;
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

  // this should changed on the loop on the cluster of the shower
  for(art::PtrVectorItr<recob::Hit> hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    time = theHit->PeakTime();  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geom->ChannelToWire(channel, tpc, plane, wire);
    
   
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

  return (void)0;
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




//for(int j = 0; j < mc->NParticles(); ++j){
  //    simb::MCParticle part(mc->GetParticle(j));


        //event_has_pi_plus=1;
      

//double vertex[3]={0,0,0};
//for( unsigned int i = 0; i < mclist.size(); ++i )
  //  art::Ptr<simb::MCTruth> mc(mclist[0]);

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





