////////////////////////////////////////////////////////////////////////
/// \file  ShowerAngleCluster.h
/// \brief
///
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author:  biagio, ART port: echurch, detector agnostic + branch into cluster module: andrzejs
////////////////////////////////////////////////////////////////////////
//#ifndef SHOWERANGLECLUSTER_H
//#define SHOWERANGLECLUSTER_H

#include "art/Framework/Core/EDProducer.h" // include the proper bit of the framework
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include <vector>
#include <string>

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "TH1F.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TH2F.h"

#include "Geometry/Geometry.h"
#include "Utilities/LArProperties.h"
#include "Utilities/GeometryUtilities.h"
#include "Utilities/DetectorProperties.h"
#include "RecoBase/Cluster.h"


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
#include "art/Framework/Core/ModuleMacros.h"

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

// LArSoft includes
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "Utilities/AssociationUtil.h"
#include "Geometry/PlaneGeo.h"
#include "SimulationBase/MCTruth.h"
#include "Utilities/GeometryUtilities.h"

// ***************** //


namespace cluster {

  class ShowerAngleCluster : public art::EDProducer {
    
  public:

    /**METHODS global*/
    explicit ShowerAngleCluster(fhicl::ParameterSet const& pset);/**Constructor*/
    virtual ~ShowerAngleCluster();                               /**Destructor*/
    void beginJob();                                     
    void beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& pset);
    void produce(art::Event& evt);                       /**Routine that finds the cluster and sets the dTdW of the 2D shower*/
   
    void   Find2DAxisRough(unsigned int NClust,std::vector < art::Ptr < recob::Hit> > hitlist); /**Calculate 2D angle histograms, provided vertex is know */ 
    void   Get2DVariables(unsigned int nClust,std::vector < art::Ptr < recob::Hit> > hitlist);   /** Calculate 2D variables to be saved into 	  */ 
    void   FitAngularDistributions(std::vector < art::Ptr < recob::Hit> > hitlist);  /** Get actual 2D omega angle from the histograms  */
  

    void   GetVertexN(art::Event& evt); /**Get vertex from MCtruth information. Temporary. */


    
    //    float LastWire[2];  // last wire of the shower
    //    float LastTime[2];  // last t_hit of the shower
    
     
    

  private:

    
    double fWiretoCm,fTimetoCm,fWireTimetoCmCm;
    
    std::vector< unsigned int >fNWires;
    double fNTimes;
    
    
    
    art::ServiceHandle<util::DetectorProperties> detp;
    art::ServiceHandle<util::LArProperties> larp;
    art::ServiceHandle<geo::Geometry> geo;
    util::GeometryUtilities gser;
   
    float fTimeTick; // time sample in us
    float fPresamplings;
    float fDriftVelocity;
    const static int    alpha           = 5;     // parameter (how many RMs (of the anglular distribution) is large the cone of the shower)
   
   
    int fRun,fEvent,fSubRun;
   

    
    std::vector< unsigned int > fMinWire,fMaxWire;
    std::vector < double > fMinTime,fMaxTime;
    
    //input parameter labels:
 
    std::string fClusterModuleLabel;
    std::string fVertexCLusterModuleLabel;
    std::string fMCGeneratorLabel;
    std::string fLarGeantlabel;
    int fUseMCVertex;
  
  
    // 2D slope and intercept of the shower axis
    std::vector <double> slope;       // in wire, time
    std::vector <double> slope_cm;       // in cm, cm
    std::vector <double> xangle;       // in cm, cm
   
    std::vector <double> lineslope;   // in wire, time

    std::vector <double>  lineinterc;   

   
    std::vector<double> xyz_vertex;
    std::vector<double> xyz_vertex_fit;

 
    double fWirePitch ;   // wire pitch in cm

    std::vector<double> fOmega_Mean;    // Mean value of the 2D angular distribution (0=Ind - 1=Coll) cm,cm
    std::vector<double> fOmega_RMS;     // RMS of the 2D angular distribution  (0=Ind - 1=Coll) cm, cm
    
    //std::vector<double> fOmega_Mean_reb;    // Mean value of the 2D angular Rebinned by 4
    //std::vector<double> fOmega_RMS_reb;     // RMS of the 2D angular distribution  Rebinned by 4
    //std::vector<double> fOmega_Mean_Mean;    // Mean value of the 2D angular use mean instead of maximum
 
   // std::vector<double> fOmega_Mean_line;    // Mean value of the 2D angular distribution obtained from linear fit
   // std::vector<double> fOmega_RMS_line;     // RMS of the 2D angular distribution obtained from linear fit

   
   // std::vector<double> fOmega_wt_Mean; // Mean value of the angular distribution (0=Ind - 1=Coll) wire,time
   // std::vector<double> fOmega_wt_RMS;  // RMS of the angular distribution  (0=Ind - 1=Coll) wire,time

  
    std::vector<std::vector<double> > fSingleEvtAngle;  // vector to show the plane omega distributions 
    std::vector<std::vector<double> > fSingleEvtAngleVal;  //vector to show the plane omega distributions


    std::vector<std::vector<double> > fShowerWidthProfile2D;  // vector to show the plane shower Width distribution 
    std::vector<std::vector<double> > fShowerChargeProfile2D;  //vector to show the plane shower Charge distribution
    std::vector<std::vector<double> > fShowerPosition2D;  //vector to store the positions of hit values stored in the previous two vectors.

    std::vector<unsigned int> fWireVertex;  // wire coordinate of vertex for each plane
    std::vector<double> fTimeVertex;  // time coordinate of vertex for each plane

    std::vector<unsigned int> fWireLast;  // wire coordinate of last point for each plane
    std::vector<double> fTimeLast;  // time coordinate of last point for each plane

    
    std::vector<double>  test_wire_start,test_time_start;

    std::vector<double> fRMS_wire;
    std::vector<double> fRMS_time;
    std::vector<double> fChisq;
    std::vector<double> fcorrelation;
    std::vector<double> fcovariance;



    // unsigned int tpc;    //tpc type
    unsigned int fNPlanes; // number of planes  

    TH1F* fh_theta[3]; 
   // std::vector< TH1F * > fh_omega_evt_reb;
    std::vector< TH1F * > fh_omega_evt;
    TH1F *  fh_omega_single;
    TH1F* fh_theta_wt[3]; 
    TTree* ftree_cluster;
    std::vector < TF1 * > linefit_cm;
    std::vector < TF1 * > linefit2_cm;	
    std::vector < TH2F * > tgx;
    std::vector < TH2F * > tgx2;

    // TH1F * hithist[3];
    void CalculateAxisParameters(unsigned nClust, std::vector < art::Ptr < recob::Hit> >  hitlist,double wstart,double tstart,double wend,double tend);
    
    void Find2DStartPoints(int nClust,std::vector< art::Ptr < recob::Hit> > hitlist);

    void Find2DBestPlanes(std::vector<int> &best_planes);
  
    void Find_Extreme_Intercepts(std::vector < art::Ptr<recob::Hit> > hitlist,
			     double perpslope,
			     double &inter_high,
			     double &inter_low);  
    
    art::Ptr<recob::Hit> FindClosestHit(std::vector < art::Ptr < recob::Hit> > hitlist,
					unsigned int wire,
					double time);  
    
    void RefineStartPoints(unsigned int nClust, std::vector< art::Ptr < recob::Hit> >  hitlist, double  wire_start,double time_start);
    
    double Get2DAngleForHit( art::Ptr<recob::Hit> starthit,std::vector < art::Ptr < recob::Hit> > hitlist);
    
    double Get2DAngleForHit( unsigned int wire, double time,std::vector < art::Ptr < recob::Hit> > hitlist);
    
    
    int GetPlaneAndTPC(art::Ptr<recob::Hit> a,unsigned int &p,unsigned int &cs,unsigned int &t,unsigned int &w);
    
    void FindMinMaxWiresTimes(unsigned int nClust,std::vector < art::Ptr < recob::Hit> >  hitlist);
    
    void ClearandResizeVectors(unsigned int nClusters);
    
 //   double Get2DDistance(double wire1,double time1,double wire2,double time2);
  
 //   double Get2DPitchDistance(double angle,double inwire,double wire);
    
 //   int GetPointOnLine(double slope,double intercept,double wire1,double time1,double &wireout,double &timeout);
    
 //   int GetPointOnLine(double slope,double wirestart,double timestart,double wire1,double time1,double &wireout,double &timeout);
    
 //   int GetPointOnLineWSlopes(double slope,double intercept,double ort_intercept,double &wireout,double &timeout);
    
  
    
    //temporary
    int mcpdg;
    double mcenergy;
    double mcphi;
    double mctheta;
    std::vector< unsigned int> mcwirevertex;  // wire coordinate of vertex for each plane
    std::vector< double> mctimevertex;  // time coordinate of vertex for each plane

    std::vector<bool> startflag;
    bool endflag;
    bool matchflag;
    
    
  }; // class ShowerAngleCluster

//   struct SortByWire 
// {
//   bool operator() (recob::Hit const& h1, recob::Hit const& h2) const 
//   { return 
//       h1.Wire()->RawDigit()->Channel() < h2.Wire()->RawDigit()->Channel() ;
//   }
// };

  
}





//------------------------------------------------------------------------------
cluster::ShowerAngleCluster::ShowerAngleCluster(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
  produces< std::vector<recob::Cluster> >();
  produces< art::Assns<recob::Cluster, recob::Hit>  >(); 
  produces< std::vector < art::PtrVector <recob::Cluster> >  >();
}


void cluster::ShowerAngleCluster::reconfigure(fhicl::ParameterSet const& pset) 
{
  fClusterModuleLabel 		=pset.get< std::string >("ClusterModuleLabel");
  fVertexCLusterModuleLabel	=pset.get<std::string > ("VertexClusterModuleLabel");
  fMCGeneratorLabel		=pset.get<std::string > ("MCGeneratorLabel");
  fLarGeantlabel		=pset.get<std::string > ("LarGeantlabel");     
 }

// ***************** //
cluster::ShowerAngleCluster::~ShowerAngleCluster()
{
}

//____________________________________________________________________________
 void cluster::ShowerAngleCluster::beginRun(art::Run& run)
  {
    
  //  std::cout << "------------- In SHowANgle preBeginRun"<<  larp->Efield() << std::endl;
  //  


    return;
  }







// ***************** //
void cluster::ShowerAngleCluster::beginJob()
{

  //temporary:
  unsigned int tpc=0;
  

  
  // this will not change on a run per run basis.
  fNPlanes = geo->Nplanes();
  fWirePitch = geo->WirePitch(0,1,0);    //wire pitch in cm
  fTimeTick=detp->SamplingRate()/1000.; 
 
  /**Get TFileService and define output Histograms*/
  art::ServiceHandle<art::TFileService> tfs;

  
    fNTimes=geo->DetHalfWidth(tpc)*2/(fTimetoCm);
    fNWires.resize(fNPlanes);
    
  for(unsigned int i=0;i<fNPlanes;++i){
   
     fNWires[i]=geo->Plane(i,tpc).Nwires();
     
  

    /**Histos for the angular distribution theta of the shower*/
    //\todo - which of these histograms are needed?
    fh_theta[i] = tfs->make<TH1F>(Form("fh_theta_%i",i),"Theta distribution",720,-180., 180.);

    
    	
    
    /**Histos for the angular distribution theta of the shower*/
    fh_omega_evt.push_back( tfs->make<TH1F>(Form("fh_omega_evt_%i",i),
					"Theta distribution per event",720,-180., 180.) );
       
  //  fh_omega_evt_reb.push_back( tfs->make<TH1F>(Form("fh_omega_evt_reb_%i",i),
  //					  "Theta distribution per event, rebinned",180,-180., 180.) );
   
    
   

    /**Histos for the angular distribution theta wire of the shower*/
   fh_theta_wt[i] = tfs->make<TH1F>(Form("ftheta_wire_%i",i),
				     "Theta wire distribution",720,-180., 180.);
  }  // end loop on planes
  
  ftree_cluster =tfs->make<TTree>("ShowerAngleCluster","Results");/**All-knowing tree with reconstruction information*/
   fh_omega_single= tfs->make<TH1F>("fh_omega_single","Theta distribution Hit",720,-180., 180.) ;
    
    ftree_cluster->Branch("run",&fRun,"run/I");
    ftree_cluster->Branch("subrun",&fSubRun,"subrun/I");
    ftree_cluster->Branch("event",&fEvent,"event/I");
    ftree_cluster->Branch("nplanes",&fNPlanes,"nplanes/I");
  
    ftree_cluster->Branch("mcpdg",&mcpdg,"mcpdg/I");
    ftree_cluster->Branch("mcenergy",&mcenergy,"mcenergy/D");
   
    ftree_cluster->Branch("mcphi",&mcphi,"mcphi/D");
    ftree_cluster->Branch("mctheta",&mctheta,"mctheta/D");
       
    ftree_cluster->Branch("wire_vertex","std::vector<unsigned int>", &fWireVertex);
    ftree_cluster->Branch("time_vertex","std::vector<double>", &fTimeVertex);

    ftree_cluster->Branch("mcwirevertex","std::vector<unsigned int>", &mcwirevertex);
    ftree_cluster->Branch("mctimevertex","std::vector<double>", &mctimevertex);
      
    ftree_cluster->Branch("wire_last","std::vector<unsigned int>", &fWireLast);
    ftree_cluster->Branch("time_last","std::vector<double>", &fTimeLast);

    ftree_cluster->Branch("test_wire_start","std::vector<double>", &test_wire_start);
    ftree_cluster->Branch("test_time_start","std::vector<double>", &test_time_start);

  //  ftree_cluster->Branch("fitw_last","std::vector<double>", &wire_end);
  //  ftree_cluster->Branch("fitt_last","std::vector<double>", &time_end); 

    ftree_cluster->Branch("xyz_vertex","std::vector<double>", &xyz_vertex);
    ftree_cluster->Branch("xyz_vertex_fit","std::vector<double>", &xyz_vertex_fit);

    ftree_cluster->Branch("omega_2d","std::vector<double>", &fOmega_Mean);
    ftree_cluster->Branch("omega_2d_RMS","std::vector<double>", &fOmega_RMS);

   

    ftree_cluster->Branch("slope","std::vector<double>", &slope);		
    ftree_cluster->Branch("lineslope","std::vector<double>", &lineslope);
   
    ftree_cluster->Branch("RMS_wire","std::vector<double>", &fRMS_wire);
    ftree_cluster->Branch("RMS_time","std::vector<double>", &fRMS_time);

    ftree_cluster->Branch("Chisq","std::vector<double>", &fChisq);
    ftree_cluster->Branch("MinWire","std::vector<unsigned int>", &fMinWire);
    ftree_cluster->Branch("MaxWire","std::vector<unsigned int>", &fMaxWire);
    ftree_cluster->Branch("MinTime","std::vector<double>", &fMinTime);
    ftree_cluster->Branch("MaxTime","std::vector<double>", &fMaxTime);
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

// ************************************* //
void cluster::ShowerAngleCluster::ClearandResizeVectors(unsigned int nClusters) {
  
  ///////////////
  fMinWire.clear();
  fMaxWire.clear();
  fMinTime.clear();
  fMaxTime.clear();
  startflag.clear();
  
 
  // delete histograms maybe?
  
 // tgx.clear();
 // tgx2.clear();
  
  ///////////////////
  
  fRMS_wire.clear();
  fRMS_time.clear();
  
  fChisq.clear();
  fcorrelation.clear();
  fcovariance.clear();
  
  lineslope.clear();
  lineinterc.clear();
  
  
  /////
  
  
  
  //fPitch.resize(0);  // Pitch calculated the old way
  fShowerWidthProfile2D.clear(); ;  // vector to show the plane shower Width distribution 
  fShowerChargeProfile2D.clear(); ;  //vector to show the plane shower Charge distribution
  fShowerPosition2D.clear(); ;  //vector to store the positions of hit values stored in the previous two vectors.
  fSingleEvtAngle.clear(); 
  fSingleEvtAngleVal.clear();

  fWireVertex.clear();
  fTimeVertex.clear();
  fWireLast.clear();
  fTimeLast.clear();
  
  slope.clear();
  xangle.clear();
  slope_cm.clear();
  
  
 
  fSingleEvtAngle.resize(nClusters); 
  fSingleEvtAngleVal.resize(nClusters); 
  fShowerWidthProfile2D.resize(nClusters); ;  // vector to show the plane shower Width distribution 
  fShowerChargeProfile2D.resize(nClusters); ;  //vector to show the plane shower Charge distribution
  fShowerPosition2D.resize(nClusters); ;  //vector to store the positions of hit values stored in the previous two

  for(unsigned int ii=0;ii<nClusters;ii++){   
    fSingleEvtAngle[ii].resize(180); 
    fSingleEvtAngleVal[ii].resize(180); 
    fShowerWidthProfile2D[ii].resize(0); ;  // vector to show the plane shower Width distribution 
    fShowerChargeProfile2D[ii].resize(0); ;  //vector to show the plane shower Charge distribution
    fShowerPosition2D[ii].resize(0); ;  //vector to store the positions of hit values stored in the
  }


  // fPitch.resize(fNPlanes); 
 	 
 
  mcwirevertex.resize(fNPlanes);  // wire coordinate of vertex for each plane 
  mctimevertex.resize(fNPlanes);  // time coordinate of vertex for each plane
  
  
  ///// ??? ???
  xyz_vertex.resize(3);
  xyz_vertex_fit.resize(3);	

  test_wire_start.resize(fNPlanes);
  test_time_start.resize(fNPlanes);

  
  
  


  fOmega_Mean.resize(fNPlanes);    // Mean value of the 2D angular distribution (1=Ind - 0=Coll) cm,cm
  fOmega_RMS.resize(fNPlanes);     // RMS of the 2D angular distribution  (1=Ind - 0=Coll) cm, cm
 

  
}
  
  
  
// ***************** //
void cluster::ShowerAngleCluster::produce(art::Event& evt)
{ 
  

  
  
 

  //%%%%% this goes into ana module.
  //Find run, subrun and event number:
  fRun = evt.id().run();
  fSubRun = evt.id().subRun();
  fEvent = evt.id().event();

  /**Get Clusters*/
  fDriftVelocity=larp->DriftVelocity(larp->Efield(),larp->Temperature());
  fWiretoCm=fWirePitch;
  fTimetoCm=fTimeTick*fDriftVelocity;
  fWireTimetoCmCm=(fTimeTick*fDriftVelocity)/fWirePitch;
 
  
 
  
  
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterModuleLabel,clusterListHandle);

  
  art::Handle< std::vector<art::PtrVector < recob::Cluster> > > clusterAssociationHandle;
  matchflag=true;  
  try{
//     art::FindManyP<recob::Cluster>  fmclust= art::FindManyP<recob::Cluster>(clusterListHandle, evt, fClusterModuleLabel);
   
    evt.getByLabel(fClusterModuleLabel,clusterAssociationHandle);
  
  std::cout << " cluster Assoc Handle size: " << clusterAssociationHandle->size() << std::endl;

 //    std::cout << " CLusters associated " << fmclust.at(0).size() << " " <<std::endl; 
      }
  catch(cet::exception e) {
      mf::LogWarning("ShowerAngleCluster") << "caught exception \n"
                                   << e;
      matchflag=false;				   
      }
  
  
  art::FindManyP<recob::Hit> fmh(clusterListHandle, evt, fClusterModuleLabel);

  //art::FindManyP<recob::Cluster> fmclust;  
  

  
 // std::vector< art::PtrVector < recob::Hit> > hitlist_all;
  //hitlist_all.resize(fNPlanes);
 
  art::PtrVector<recob::Cluster> clusters;


  std::cout << " ++++ Clusters received " << clusterListHandle->size() << " +++++ " << std::endl;
  if(clusterListHandle->size() ==0 )
  {
    std::cout << " no clusters received! exiting " << std::endl;
    return;
  }
  ClearandResizeVectors(clusterListHandle->size());
  // resizing once cluster size is known.
  
  
  
  endflag=false;
 // GetVertexN(evt);
  
  
  for(unsigned int iClust = 0; iClust < clusterListHandle->size(); iClust++){

    art::Ptr<recob::Cluster> cl(clusterListHandle, iClust);
    std::vector< art::Ptr<recob::Hit> > hitlist = fmh.at(iClust);
    
   // std::cout << "++++ hitlist size " << hitlist.size() << std::endl;
    
    unsigned int p(0),w(0), t(0),cs(0); //c=channel, p=plane, w=wire
    GetPlaneAndTPC(hitlist[0],p,cs,t,w);
    
   // if(hitlist.size()>=5){
      clusters.push_back(cl);

      std::vector< double > spos=cl->StartPos();
      std::vector< double > sposerr=cl->SigmaStartPos();
      
      std::vector< double > epos=cl->EndPos();
      std::vector< double > eposerr=cl->SigmaEndPos();
      
     
      if(spos[0]!=0 && spos[1]!=0 && sposerr[0]==0 && sposerr[1]==0 ){
	fWireVertex.push_back(spos[0]);
	fTimeVertex.push_back(spos[1]);
	startflag.push_back(true);
	//std::cout << "setting external starting points " << spos[0] << " " << spos[1] <<" " << sposerr[0] <<" "<< sposerr[1] << std::endl;
      }
	  else
	startflag.push_back(false); 
	
      if(epos[0]!=0 && epos[1]!=0 && eposerr[0]==0 && eposerr[1]==0 ){
	fWireLast.push_back(epos[0]);
	fTimeLast.push_back(epos[1]);
	endflag=true;
	//std::cout << "setting external ending points " << epos[0] << " " << epos[1] << std::endl;
      }
	
    //  for(std::vector < art::Ptr < recob::Hit > >::const_iterator a = hitlist.begin(); a != hitlist.end();  a++){ //loop over cluster hits
//	GetPlaneAndTPC(*a,p,cs,t,w);
//          hitlist_all.push_back(hitlist);
     //use MC vertex:
     /////////// ----------- to remove - testing purposes
//      if(fWireVertex.size()<iClust && iClust<fNPlanes){
//        startflag[iClust]=true;
//        fWireVertex.push_back(mcwirevertex[iClust]);
//      }
//      else{
//        fWireVertex[iClust]=mcwirevertex[iClust];
//         startflag[iClust]=true;
//      }
// 	
//      if(fTimeVertex.size()<iClust && iClust<fNPlanes)
//       fTimeVertex.push_back(mctimevertex[iClust]);
//      else
//        fTimeVertex[iClust]=mctimevertex[iClust];

     /////////////////------------- end of to remove
     

	  
  // find min and max wires for this set of clusters: 
      fMinWire.push_back(9999),fMaxWire.push_back(0);
      fMinTime.push_back(99999),fMaxTime.push_back(0);
      FindMinMaxWiresTimes(iClust,hitlist);

      // Get rough estimate of 2D angle.    
      Find2DAxisRough(iClust,hitlist); 	  
  
      Find2DStartPoints(iClust, hitlist);
      
      xangle.push_back(Get2DAngleForHit( fWireVertex[iClust],fTimeVertex[iClust], hitlist));
      slope_cm.push_back(tan((xangle[iClust]*TMath::Pi()/180)));
      slope.push_back(slope_cm[iClust]/fWireTimetoCmCm);
      
      Get2DVariables(iClust, hitlist);
      
      
    } // End loop on clusters.
  

   
//    // this goes into the ana module.
    

  
 
  // make an art::PtrVector of the clusters
  std::unique_ptr<std::vector<recob::Cluster> > ShowerAngleCluster(new std::vector<recob::Cluster>);
  std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);

  for(unsigned int iplane=0;iplane<fNPlanes;iplane++){
    std::vector< art::Ptr<recob::Hit> > hitlist = fmh.at(iplane);
    if(hitlist.size()>0) {
      
    double wverror=fWireVertex[iplane]*0.05,tverror=fTimeVertex[iplane]*0.05;
    
    if(startflag[iplane])
    {
     wverror=0;
     tverror=0;
    }
      
      //std::cout << " saving cluster with errors, plane: "<< iplane <<" " << wverror << " " << tverror << std::endl;
    
    recob::Cluster temp(fWireVertex[iplane], wverror,
                         fTimeVertex[iplane], tverror,
                         fWireLast[iplane], fWireLast[iplane]*0.05,
                         fTimeLast[iplane], fTimeLast[iplane]*0.05,  
		    xangle[iplane], xangle[iplane]*0.05, lineslope[iplane],lineinterc[iplane],5.,
		    geo->Plane(iplane,0,0).View(),
		    iplane);
  
      std::cout << " Saving cluster for plane: " << iplane << " w,t " << fWireVertex[iplane] << " " << fTimeVertex[iplane] << " 2D angle: " <<  xangle[iplane] << " lslope " << lineslope[iplane] << std::endl;

    
    ShowerAngleCluster->push_back(temp);
    // associate the hits to this cluster
    util::CreateAssn(*this, evt, *(ShowerAngleCluster.get()), hitlist, *(assn.get()));
  // std::cout << "######## in plane loop filling clusters " << std::endl;
    }

    
  }

  /////////////////////////////////////////////
  std::unique_ptr< std::vector < art::PtrVector < recob::Cluster > > > classn(new std::vector < art::PtrVector < recob::Cluster > >);
  
  if(!matchflag)
     {
      //matching code here.
       
    art::PtrVector < recob::Cluster > cvec;
	
    for(unsigned int ip=0;ip<fNPlanes;ip++)  {
	art::ProductID aid = this->getProductID< std::vector < recob::Cluster > >(evt);
	art::Ptr< recob::Cluster > aptr(aid, ip, evt.productGetter(aid));
	cvec.push_back(aptr);
      }
     classn->push_back(cvec); 
     }
  else
    {
    for(unsigned int i=0;i<clusterAssociationHandle->size();i++)  
    classn->push_back((*clusterAssociationHandle)[i]);
    }

  /**Fill the output tree with all information */
  ftree_cluster->Fill();

  evt.put(std::move(ShowerAngleCluster));
  evt.put(std::move(assn));
  evt.put(std::move(classn));
}


// ******************************* //
int cluster::ShowerAngleCluster::GetPlaneAndTPC(art::Ptr<recob::Hit> a,
						unsigned int &p,
						unsigned int &cs,
						unsigned int &t,
						unsigned int &w)
{
  unsigned int c = a->Wire()->RawDigit()->Channel(); 
  geo->ChannelToWire(c,cs,t,p,w);
    
  return 0;
}




//fMinWire and fMaxWire should be accessible elsewhere.
// \todo


void cluster::ShowerAngleCluster::FindMinMaxWiresTimes(unsigned int nClust,std::vector< art::Ptr < recob::Hit> > hitlist){

  unsigned int wire,tpc, cstat;
  unsigned int plane;

  
   for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    double time = (*hitIter)->PeakTime();  
    //time_C -= (presamplings+10.1);
    GetPlaneAndTPC((*hitIter),plane,cstat,tpc,wire);
    
    if(time>fMaxTime[nClust])
	fMaxTime[nClust]=time;

    if(time<fMinTime[nClust])
	fMinTime[nClust]=time;
    
    if(wire>fMaxWire[nClust])
	fMaxWire[nClust]=wire;

    if(wire<fMinWire[nClust])
	fMinWire[nClust]=wire;

  }
// std::cout << " +++ maximums " << fMinWire[nClust] << " " << fMaxWire[nClust] << " " <<fMinTime[nClust] << " " << fMaxTime[nClust] << std::endl;
}


//*********************************//
// Angular distribution of the energy of the shower - Collection view
// 
///////////////////////////////

void cluster::ShowerAngleCluster::Find2DAxisRough(unsigned int nClust,std::vector < art::Ptr < recob::Hit> > hitlist){
 
  double time;
  unsigned int wire,tpc, cstat;
  unsigned int plane;

  art::ServiceHandle<art::TFileService> tfs;
  
  if(hitlist.size()==0)
    return;
  
 
 // padding of the selected TGraph in wires and time. 
  int wirepad=20;
  int timepad=wirepad*fWiretoCm/fTimetoCm+0.5;
 
  int nbinsx= (fMaxWire[nClust]-fMinWire[nClust]+2*wirepad)*fWiretoCm;  // nbins to have 
  int nbinsy= (fMaxTime[nClust]-fMinTime[nClust]+2.*(double)timepad)*fTimetoCm;  // nbins to have 
 
 

//   std::cout << " +++ maximums " << fMinWire[nClust] << " " << fMaxWire[nClust] << " " <<fMinTime[nClust] << " " << fMaxTime[nClust] << std::endl;
//  std::cout << "------ NClust" << nClust+1 << " " << tgx.size() << " if: " << ((nClust+1) > tgx.size()) << std::endl;
//  
//  std::cout << nbinsx<< " " << ((double)fMinWire[nClust]-(double)wirepad)*fWiretoCm<< " " << 			((double)fMaxWire[nClust]+(double)wirepad)*fWiretoCm<< " " << nbinsy<< " " << 
// 			  (fMinTime[nClust]-timepad)*fTimetoCm<< " " << (fMaxTime[nClust]+timepad)*fTimetoCm <<std::endl;
//  
  // argh. Stupid unsigned int.
  if((nClust+1)>tgx.size()) {
    std::cout << " inside true if " << std::endl;
   tgx.push_back(tfs->make<TH2F>(Form("charge distrib_%i",nClust),"charge distribution per wires",
			  nbinsx,((double)fMinWire[nClust]-(double)wirepad)*fWiretoCm,			((double)fMaxWire[nClust]+(double)wirepad)*fWiretoCm,nbinsy,
			  (fMinTime[nClust]-timepad)*fTimetoCm,(fMaxTime[nClust]+timepad)*fTimetoCm));
    
  tgx2.push_back(tfs->make<TH2F>(Form("hit distrib_%i",nClust),"Hit distribution per wires",
			    nbinsx,((double)fMinWire[nClust]-(double)wirepad)*fWiretoCm,			((double)fMaxWire[nClust]+(double)wirepad)*fWiretoCm,nbinsy,
			    (fMinTime[nClust]-timepad)*fTimetoCm,(fMaxTime[nClust]+timepad)*fTimetoCm));
  
  linefit_cm.push_back(tfs->make<TF1>(Form("linefit_cm_%d",nClust),"pol1",
				      ((double)fMinWire[nClust]-(double)wirepad)*fWiretoCm,					((double)fMaxWire[nClust]+(double)wirepad)*fWiretoCm));  
  linefit2_cm.push_back(tfs->make<TF1>(Form("linefit_2_cm_%d",nClust),"pol1",
				       ((double)fMinWire[nClust]-(double)wirepad)*fWiretoCm,					((double)fMaxWire[nClust]+(double)wirepad)*fWiretoCm));
  
  }
  else {
    tgx[nClust]->Reset();
    tgx[nClust]->SetBins(nbinsx,((double)fMinWire[nClust]-(double)wirepad)*fWiretoCm,			((double)fMaxWire[nClust]+(double)wirepad)*fWiretoCm,nbinsy,
			  (fMinTime[nClust]-timepad)*fTimetoCm,(fMaxTime[nClust]+timepad)*fTimetoCm);

    tgx2[nClust]->Reset();
    tgx2[nClust]->SetBins(nbinsx,((double)fMinWire[nClust]-(double)wirepad)*fWiretoCm,			((double)fMaxWire[nClust]+(double)wirepad)*fWiretoCm,nbinsy,
			  (fMinTime[nClust]-timepad)*fTimetoCm,(fMaxTime[nClust]+timepad)*fTimetoCm);
    linefit_cm[nClust]->SetRange(((double)fMinWire[nClust]-(double)wirepad)*fWiretoCm,										((double)fMaxWire[nClust]+(double)wirepad)*fWiretoCm);
    linefit2_cm[nClust]->SetRange(((double)fMinWire[nClust]-(double)wirepad)*fWiretoCm,										((double)fMaxWire[nClust]+(double)wirepad)*fWiretoCm);
   }
  
  
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != 				hitlist.end();  hitIter++){
    
    time =  (*hitIter)->PeakTime();  
    GetPlaneAndTPC((*hitIter),plane,cstat,tpc,wire);
  
   // std::cout << " -- wire, time, charge" << wire*fWiretoCm << " " << time*fTimetoCm << " " << (*hitIter)->Charge() << std::endl;
    
    tgx[nClust]->Fill((double)wire*fWiretoCm,
		     time*fTimetoCm,(*hitIter)->Charge());
    tgx2[nClust]->Fill((double)wire*fWiretoCm,time*fTimetoCm);		
  }
  
 
  tgx[nClust]->Fit(Form("linefit_cm_%d",nClust),"QMRNCFrob=0.8");
  tgx2[nClust]->Fit(Form("linefit_2_cm_%d",nClust),"QMRNCFrob=0.95");


//   std::cout << "{{{-----}}}  histo stats: rms w,t " << tgx[nClust]->GetRMS(1) << " " << 
//   tgx[plane]->GetRMS(2) << " chisq " << linefit_cm[nClust]->GetChisquare()/linefit_cm[nClust]->GetNDF()
//   << " max, min wires and times " <<  fMinWire[nClust] << " " <<fMaxWire[nClust] << " " <<  fMinTime[nClust] << " " << fMaxTime[nClust] << std::endl;


  fRMS_wire.push_back(tgx[nClust]->GetRMS(1));
  fRMS_time.push_back(tgx[nClust]->GetRMS(2));
  
  fChisq.push_back(linefit_cm[nClust]->GetChisquare()/linefit_cm[nClust]->GetNDF());
  fcorrelation.push_back(tgx[nClust]->GetCorrelationFactor());
  fcovariance.push_back(tgx[nClust]->GetCovariance());
  
  lineslope.push_back(linefit_cm[nClust]->GetParameter(1)/fWireTimetoCmCm);
  lineinterc.push_back(linefit_cm[nClust]->GetParameter(0)/fTimetoCm);

  
  
}








// ***************** //
// void cluster::ShowerAngleCluster::FitAngularDistributions(art::PtrVector < recob::Hit> hitlist){
//   /** Fit function of the angular distribution (cm,cm)*/
//  
//  
//   double time;
//   unsigned int wire;
//   double BC,AC;
//   double omega;
//   unsigned int channel,iplane,plane,tpc,cstat;
// 
//   if(hitlist.size()==0)
//     return;
// 
//     art::Ptr<recob::Hit> theHit = (*hitlist.begin());
//     time = theHit->PeakTime();  
//     //time_C -= (presamplings+10.1);
//     art::Ptr<recob::Wire> theWire = theHit->Wire();
//     channel = theWire->RawDigit()->Channel();
//     geo->ChannelToWire(channel, cstat, tpc, iplane, wire);
// 
// 
// 
//   
// 
//   // this should changed on the loop on the cluster of the shower
//   for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
//     art::Ptr<recob::Hit> theHit = (*hitIter);
//     time = theHit->PeakTime();  
//     //time_C -= (presamplings+10.1);
//     art::Ptr<recob::Wire> theWire = theHit->Wire();
//     channel = theWire->RawDigit()->Channel();
//     geo->ChannelToWire(channel, cstat, tpc, plane, wire);
//   
//       BC = ((double)wire - fWireVertex[plane])*fWiretoCm; // in cm
//       AC = ((double)time - fTimeVertex[plane])*fTimetoCm; //in cm 
//       omega = asin(  AC/sqrt(pow(AC,2)+pow(BC,2)) );
//  
//      double omx=gser.Get2Dangle(((double)wire - fWireVertex[plane]),((double)time - fTimeVertex[plane]));
//       
//       
//       if(BC<0)  // for the time being. Will check if it works for AC<0
// 	  { 
// 	  if(AC!=0)
// 	  omega= AC/fabs(AC)*TMath::Pi()-omega;  //subtract when negative, add when positive
// 	  else    
// 	  omega=TMath::Pi();
// 	  } 
// 
//       omega = 180*omega/TMath::Pi();
//       
// 
//       std::cout << " omega and omx "<< omega << " "<< omx << std::endl;
//       
//       fh_theta[plane]->Fill(omega, theHit->Charge()); // Filling the histo (angle, energy of the hit)
//       fh_omega_evt[plane]->Fill(omega, theHit->Charge());
// 
//   
//       //fh_omega_evt_reb[plane]->Fill(omega, theHit->Charge());
// 
//   }
//   //for(unsigned int iplane = 0; iplane < fNPlanes; ++iplane){
//   
// 
//     fOmega_Mean[iplane] =
//     fh_omega_evt[iplane]->GetBinCenter(fh_omega_evt[iplane]->GetMaximumBin());// Mean value of the fit
//     fOmega_RMS[iplane] = fh_omega_evt[iplane]->GetRMS(); // RMS of the fit of the angular distribution in deg
//    
// 
// // \todo this should be a parameter not hardcoded.
// for(int i=0;i<720;i+=720/180)
// {
//   double locsum=0;
//   for(int j=0;j<720/180;j++)
//     locsum+=fh_omega_evt[iplane]->GetBinContent(i+j);
//   
//   fSingleEvtAngleVal[iplane][i]=locsum;
// //fSingleEvtAngle[iplane][i]=fh_omega_evt[iplane]->GetBinContent(i);
//   fSingleEvtAngle[iplane][i]=(double)i/4-180;
// }
//   //}
// 
// xangle[iplane]=fOmega_Mean[iplane]*TMath::Pi()/180;
// slope[iplane] = tan((fOmega_Mean[iplane]*TMath::Pi()/180))/fWireTimetoCmCm;
// 
// 
// 
// 
// }



// this goes to the one above?

double cluster::ShowerAngleCluster::Get2DAngleForHit( unsigned int swire,double stime,std::vector < art::Ptr < recob::Hit> > hitlist) {
  
  fh_omega_single->Reset();
  
  unsigned int channel,plane,tpc,cstat,wire;
  // this should changed on the loop on the cluster of the shower
   for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
      art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime();  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = theHit->Wire();
    channel = theWire->RawDigit()->Channel();
    geo->ChannelToWire(channel, cstat, tpc, plane, wire);
  
    double omx=gser.Get2Dangle((double)wire,(double)swire,time,stime);
      
      
    fh_omega_single->Fill(180*omx/TMath::Pi(), theHit->Charge());
  
   }
  //for(unsigned int iplane = 0; iplane < fNPlanes; ++iplane){
  
double omega = fh_omega_single->GetBinCenter(fh_omega_single->GetMaximumBin());// Mean value of the fit
   
return omega; // in degrees.
  
  
  
  
  
}



double cluster::ShowerAngleCluster::Get2DAngleForHit( art::Ptr<recob::Hit> starthit,std::vector < art::Ptr < recob::Hit> > hitlist)
{
//  fh_omega_single->Reset();
  
    unsigned int channel,iplane,tpc,cstat,swire;
    double stime = starthit->PeakTime();  
    //time_C -= (presamplings+10.1);
    art::Ptr<recob::Wire> theWire = starthit->Wire();
    channel = theWire->RawDigit()->Channel();
    geo->ChannelToWire(channel, cstat, tpc, iplane, swire);

return Get2DAngleForHit( swire,stime,hitlist);  // in degrees

}




//to be moved to an _ana module
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


double drifttick=(xyz_vertex[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./fTimeTick);

const double origin[3] = {0.};
for(unsigned int iplane=0;iplane<fNPlanes;iplane++)
{
double pos[3];
 unsigned int  wirevertex, t, cs;
unsigned int p;
geo->Plane(iplane).LocalToWorld(origin, pos);
	//planex[p] = pos[0];
std::cout << "plane X positionp " << iplane << " " << pos[0] << std::endl;

pos[1]=xyz_vertex[1];
pos[2]=xyz_vertex[2];
 unsigned int channel2 = geo->NearestChannel(pos,iplane);
 geo->ChannelToWire(channel2,cs,t,p,wirevertex); 
       
if(iplane!=p)
	{std::cout << " error - planes don't match " << iplane << " " << p << std::endl;
	return ;
	}

mcwirevertex[p]=wirevertex;  // wire coordinate of vertex for each plane
mctimevertex[p]=drifttick-(pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./fTimeTick)+detp->TriggerOffset();  // time coordinate of vertex for each plane

//fWireVertex[p]=wirevertex;
//fTimeVertex[p]=drifttick-(pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./fTimeTick)+60;
std::cout<<"wirevertex= "<<wirevertex<< " timevertex " << mctimevertex[p] << " correction "<< (pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./fTimeTick) << " " << pos[0] <<std::endl;
 

}



  return (void)0;
}




//int cluster::ShowerAngleCluster::Get2Dvariables(float Wire_vertexI_wt, float Wire_vertexC_wt, float Time_I_wt, float Time_C_wt){
void cluster::ShowerAngleCluster::Get2DVariables(unsigned int nClust,std::vector < art::Ptr < recob::Hit> > hitlist)  {  

    
//get parameters of the slope obtained by searching for the maximum and start points
double tst=fTimeVertex[nClust];
double wst= fWireVertex[nClust];
double slp=slope[nClust];

double wireend=fWireLast[nClust];
//double timeend=time_end[iplane];

//double intercept=tst-slp*(double)wst;



 // std::cout << "========= line params, inside 2d variables, plane: a,c " << nClust <<" " << slp << " " << intercept << std::endl;


double projlength=gser.Get2DPitchDistance(xangle[nClust],wst,wireend);


int nbins = (hitlist.size() > 1000) ? 100 : hitlist.size()/10;

//std::cout << "---- projected length for plane: " << nClust << " " << projlength << " nbins " << nbins <<  std::endl;

TH1F * hithist=new TH1F(Form("hithist_ev_%d_pl_%d",fEvent,nClust),Form("hithist_ev_%d_pl_%d",fEvent,nClust),nbins,0.,projlength);  
  
TH1F * hithist2=new TH1F(Form("hithist_ev_%d_pl_%d_w",fEvent,nClust),Form("hithist_ev_%d_pl_%d",fEvent,nClust),nbins,0.,projlength);  

// start loop to calculate the profile of the shower along the shower axis.

//double extreme_intercept_end=-999999;
//double extreme_intercept_start=999999;

//double extr_wire_pos=0,extr_time_pos=0;

// int multiplier=1;   // +1 for positive angles, -1 for negative angles. to compensate that we are looking for either the highest (omega >0 ) or lowest (omega<0) intercept.
// 

  

for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    unsigned int wire,cstat, tpc,plane;
    GetPlaneAndTPC(theHit,plane,cstat,tpc, wire);
    
    double wire_on_line,time_on_line;
    
    gser.GetPointOnLine(slp,wst,tst,wire,time,wire_on_line,time_on_line);
    
  
    double linedist=gser.Get2DDistance(wire_on_line,time_on_line,wst,tst);

    double ortdist=gser.Get2DDistance(wire_on_line,time_on_line,wire,time);
    
   
    
    hithist->Fill(linedist);
    hithist2->Fill(linedist,ortdist);
    
    fShowerPosition2D[nClust].push_back(linedist);  
    fShowerWidthProfile2D[nClust].push_back(ortdist);
    fShowerChargeProfile2D[nClust].push_back(theHit->Charge()); 
  }
 

 TFile * trh=new TFile("histos.root","UPDATE");
 
hithist->Write();
hithist2->Write();
 
 delete trh;
 delete hithist;
 

  return (void)0;
}




void   cluster::ShowerAngleCluster::CalculateAxisParameters(unsigned nClust, std::vector < art::Ptr < recob::Hit> >  hitlist,double wstart,double tstart,double wend,double tend)
{
//need slope and intercept. Calculate Sum of vertical distances and distribution histogram of points along the axis.
unsigned int plane,cs,t,w;
GetPlaneAndTPC(hitlist[0],plane,cs,t,w);
  
//get parameters of the slope obtained by searching for the maximum and start points
double wst=wstart;
double tst= tstart;

double slp=lineslope[nClust];
double slp_cm=lineslope[nClust]*fWireTimetoCmCm;
double wireend=wend;
//double timeend=time_end[iplane];

double intercept=lineinterc[nClust];


//get slope of lines orthogonal to those found crossing the shower.
//  double aprim=0;
//       
// 	if(slp)	
// 	{
// 	aprim=-1./slp;
// 	}

 // std::cout << "========= line params, inside calculate parameters, nClust, plane: a,c " <<nClust<< " " << plane <<" " << slp << " " << intercept << " |||| wst,tst  "<<wst << " " << tst << std::endl;


//calculate from max intercepts also to check
double projlength=gser.Get2DPitchDistanceWSlope(slp_cm,wst,wireend);
//std::cout << "========= projlength, " << projlength << std::endl;

//int nbins = (hitlist.size() > 1000) ? 100 : hitlist.size()/10;
int nbins = 10;

//std::cout << "---- projected length for plane: " << " "<<nClust << " " << plane << " " << projlength << " nbins " << nbins <<  std::endl;

TH1F * hhist=new TH1F(Form("hhist_ev_%d_pl_%d",fEvent,nClust),Form("hhist_ev_%d_pl_%d",fEvent,nClust),nbins,0.,projlength);  
  
TH1F * hhist2=new TH1F(Form("hhist_ev_%d_pl_%d_w",fEvent,nClust),Form("hhist_ev_%d_pl_%d",fEvent,nClust),nbins,0.,projlength);  


  

for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    unsigned int wire,cstat, tpc,plane;
    GetPlaneAndTPC(theHit,plane, cstat, tpc,  wire);
 
        
    double wire_on_line,time_on_line;
    
    gser.GetPointOnLine(slp,intercept,wire,time,wire_on_line,time_on_line);
    
    double linedist=gser.Get2DDistance(wire_on_line,time_on_line,wst,tst);
    
    double ortdist=gser.Get2DDistance(wire_on_line,time_on_line,wire,time);
    
    hhist->Fill(linedist);
    hhist2->Fill(linedist,ortdist);
    
  //  fShowerPosition2D[plane].push_back(linedist);  
  //  fShowerWidthProfile2D[plane].push_back(ortdist);
  //  fShowerChargeProfile2D[plane].push_back(theHit->Charge()); 
  }
 

 TFile * trh=new TFile("histos_start.root","UPDATE");
 
hhist->Write();
hhist2->Write();
 
 delete trh;
 delete hhist;
 delete hhist2;
 
  
}



//////////////////////////removed from 2D StartPoints due to reorganization.
//// find for which planes the correlation factor is largest (and preferably larger than 0.6)
// std::vector<double>  wire_start,wire_end;
//   std::vector<double>  time_start,time_end;
//  
//   wire_start.resize(fNPlanes);wire_end.resize(fNPlanes);
//   time_start.resize(fNPlanes);time_end.resize(fNPlanes);
//   
//   std::vector< int > best_planes;
// 
//   Find2DBestPlanes( best_planes);
//
//

//   //error checking:
//   if(iplane!=plane){
//     std::cout << " error: planes mismatch  " << iplane << " "<< plane << std::endl;
//     return;
//     }
//
 // temporary cross-check  
//   for(unsigned int ii=0;ii<best_planes.size();ii++)
// 	std::cout << " ----++++ determined start wire points per planes " << best_planes[ii] << " "  << wire_start[best_planes[ii]] << " " << time_start[best_planes[ii]] <<std::endl;  
// 
//   //sort the best planes in order of correlation factor
//   //first sort
//   if(fabs(tgx[best_planes[0]]->GetCorrelationFactor())<fabs(tgx[best_planes[1]]->GetCorrelationFactor()))
//     std::swap(best_planes[0],best_planes[1]);
// 
//   //second sort
//   if(best_planes.size()>2 && 
//     (fabs(tgx[best_planes[1]]->GetCorrelationFactor())<fabs(tgx[best_planes[2]]->GetCorrelationFactor())))
// 	std::swap(best_planes[1],best_planes[2]);
// 
//   // sanity check - see if times of the points found are close enough.
// 	/////////////////
//   for(unsigned int ii=0;ii<fNPlanes;ii++)	
// 	std::cout << " ---- plane parameters " << ii << " " << tgx[ii]->GetCorrelationFactor() << " RMS ratio " << tgx[ii]->GetRMS(1) /tgx[ii]->GetRMS(2) << " Chisq " << linefit_cm[ii]->GetChisquare()/linefit_cm[ii]->GetNDF() << " covariance " << tgx[ii]->GetCovariance() << std::endl;
// 	
// 	
//   const double origin[3] = {0.};
//   std::vector <std::vector <  double > > position;
// 
//   // get starting positions for all planes
//   for(unsigned int xx=0;xx<fNPlanes;xx++){
//     double pos1[3];
//     geo->Plane(xx).LocalToWorld(origin, pos1);
//     std::vector <double > pos2;
//     pos2.push_back(pos1[0]);
//     pos2.push_back(pos1[1]);
//     pos2.push_back(pos1[2]);
//     position.push_back(pos2);
//   }
//////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//   //loop to check time discrepancies between points found
//   for(unsigned int xx=0;xx<bestplanes.size()-1;xx++){  // best_planes.size()-1 because we want to always find a next plane
//     for(unsigned int yy=xx+1;yy<best_planes.size();yy++){
// 	std::cout << "**** difference between planes X position, planes: "<< best_planes[xx] << " "<< best_planes[yy] << " "<< position[best_planes[xx]][0] << " " << position[best_planes[yy]][0] << " " << fabs(position[best_planes[xx]][0]-position[best_planes[yy]][0]) << " " << " " << time_start[best_planes[xx]] << " " << time_start[best_planes[yy]] << " " << fabs(time_start[best_planes[xx]]-time_start[best_planes[yy]]) << std::endl;
// 			
// 	if(fabs(time_start[best_planes[xx]]-time_start[best_planes[yy]]) > 1.5*fabs(position[best_planes[xx]][0]-position[best_planes[yy]][0])){
// 		std::cout << " !!!! Something's wrong in the wire determination " << fabs(time_start[best_planes[xx]]-time_start[best_planes[yy]]) << " " << 1.5*fabs(position[xx][0]-position[yy][0]) << std::endl;
// 	}	
//     } // end inner yy loop
//   } // end outer xx loop
// 
// 	
// 	
// 
// 	if((fabs(time_start[best_planes[0]]-time_start[best_planes[1]]) > 1.5*fabs(position[0][0]-position[1][0])) && best_planes.size() > 2) //time discrepancy of best correlation factor pair and we have a choice.
// 		{std::cout << " time discrepancy of best correlation factor pair 0 and 1 " << time_start[best_planes[0]] << " " << time_start[best_planes[1]] << " "<< position[0][0] << " " << position[1][0] << std::endl;
// 		if(!(fabs(time_start[best_planes[0]]-time_start[best_planes[2]]) > 2.5*fabs(position[0][0]-position[2][0]))) //0,1 is bad but 0,2 is ok:
// 			{
// 			std::cout << " but 0,2 is ok " << time_start[best_planes[0]] << " " << time_start[best_planes[2]] << " "<< position[0][0] << " " << position[2][0] << std::endl;
// 			std::swap(best_planes[1],best_planes[2]);
// 			}
// 		else   //0,1 is not ok and 0,2 is not ok
// 			{
// 			std::cout << " 0,1 and 0,2 is not ok " << std::endl;
// 			if(!(fabs(time_start[best_planes[1]]-time_start[best_planes[2]]) > 2.5*fabs(position[1][0]-position[2][0]))) //0,1 and 0,2 is bad but 1,2 is ok.
// 				{
// 				std::cout << " but 1,2 is ok " << time_start[best_planes[1]] << " " << time_start[best_planes[2]] << " "<< position[1][0] << " " << position[2][0] << std::endl;
// 				std::swap(best_planes[0],best_planes[1]);
// 				std::swap(best_planes[1],best_planes[2]);  // shift zero to last position.
// 				}
// 			}
// 
// 		}
//////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

//  // Assuming there is no problem ( and we found the best pair that comes close in time )
//   // we try to get the Y and Z coordinates for the start of the shower. 
// 	int chan1=geo->PlaneWireToChannel(best_planes[0],wire_start[best_planes[0]], 0);
// 	int chan2=geo->PlaneWireToChannel(best_planes[1],wire_start[best_planes[1]], 0);
// 
// 	double y,z;
// 	bool wires_cross = geo->ChannelsIntersect(chan1,chan2,y,z);
// 
// 	
// 	xyz_vertex_fit[1]=y;
// 	xyz_vertex_fit[2]=z;
// 	xyz_vertex_fit[0]=(time_start[best_planes[0]]-detp->TriggerOffset()) *fDriftVelocity*fTimeTick+position[0][0];
// 
// 
// 	std::cout << ":::::: found x,y,z vertex " << wires_cross << " " << xyz_vertex_fit[0] << " " << y << " " << z << std::endl;
// 
// 
// 	// assume some condition - probably that correlation factor is too small. Then project the found vertex in x,y,z into wire, time coordinates in the last plane.
// ///////////////// Only do the following part if there are 3 planes.
//////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
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

	geo->Plane(worst_plane).LocalToWorld(origin, pos);
	//planex[p] = pos[0];
	std::cout << "plane X positionp " << worst_plane << " " << pos[0] << std::endl;*/

	
	
	///////////////////////////////////////
	// geometry test:
// 	double width  = 2.*geo->TPC(0).HalfWidth();  //notice the geometry gives the 1/2 width, so multiply by 2
// 	double height = 2.*geo->TPC(0).HalfHeight(); //notice the geometry gives the 1/2 height, so multiply by 2
// 	double length =    geo->TPC(0).Length();     //notice the geometry gives the total length
// 	
// 	
// 	std::cout << "-------- height " << height << " " << length << " " <<width << std::endl;
// 	
// 	//display first and last wires:
// 	
// 	for(unsigned int iplane=0;iplane<fNPlanes;iplane++)
// 	    { std::cout << " +++ pl" << iplane << geo->Plane(iplane).Nwires()  << std::endl;
// 	      double wire1_Start[3]={0},wire1_End[3]={0};
// 	    geo->WireEndPoints(0,iplane,0,wire1_Start,wire1_End);
// 	      std::cout << "++++++++ wire positions: w nr " << 0 << " " << wire1_Start[0] << " " << wire1_Start[1] << " " << wire1_Start[2] << " "<< wire1_End[0] << " " << wire1_End[1] << " "<< wire1_End[2] << " " << std::endl;
// 	
// 	      geo->WireEndPoints(0,iplane,geo->Plane(iplane).Nwires()-1,wire1_Start,wire1_End);
// 	      std::cout << "++++++++ wire positions: w nr " << geo->Plane(iplane).Nwires()-1 << " " << wire1_Start[0] << " " << wire1_Start[1] << " " << wire1_Start[2] << " "<< wire1_End[0] << " " << wire1_End[1] << " "<< wire1_End[2] << " " << std::endl;
// 	
// 	      
// 	
// 	    }
//////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

	
//     pos[1]=xyz_vertex_fit[1];
//     pos[2]=xyz_vertex_fit[2];
//     unsigned int channel2 = geo->NearestChannel(pos,worst_plane);
//     geo->ChannelToWire(channel2,cstat,t,worst_plane,wirevertex); 
// 
// 
//     double drifttick=(xyz_vertex_fit[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./fTimeTick);
// 
// 
//     double timestart=drifttick-(pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./fTimeTick)+detp->TriggerOffset();
// 	std::cout << " worst plane " << worst_plane <<" wirevertex= "<<wirevertex<< " timevertex " << timestart << " correction " << detp->TriggerOffset() << " " << (pos[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./fTimeTick) << " "<<pos[0] <<std::endl;
// 
// 
// 	double min_dist=999999.;
// 
// 	for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist_all[worst_plane].begin(); hitIter != hitlist_all[worst_plane].end();  hitIter++){
//     		art::Ptr<recob::Hit> theHit = (*hitIter);
//     		time = theHit->PeakTime() ;  
//     		unsigned int plane;
//     		GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
// 		
// 	
//     		double dist_begin=gser.Get2DDistance(wirevertex,timestart,wire,time);
//     		//TMath::Sqrt( pow((double)((int)wirevertex-(int)wire)*fWirePitch,2)+pow((timestart-time)*fDriftVelocity*fTimeTick,2) );	
// 
// 		//std::cout << "=== min_dist " << wire << " " << time <<" " << dist_begin << " " << pow((double)((int)wirevertex-(int)wire)*fWirePitch,2) << " " << ((int)wirevertex-(int)wire)*fWirePitch << " " << min_dist << std::endl; 
// 		
// 		if(dist_begin<min_dist)
// 			{
// 			min_dist=dist_begin;
// 			wire_start[worst_plane]=wire;
// 			time_start[worst_plane]=time;
// 
// 			}
// 
// 
// 	} // end loop on hits.
// 
// 
// } // end big if(fNPlanes >= 3)





void   cluster::ShowerAngleCluster::Find2DStartPoints(int nClust,std::vector< art::Ptr < recob::Hit> > hitlist)
{


  
  double  wire_start=0.,wire_end=0.;
  double  time_start=0.,time_end=0.;
 
//  double time;
  unsigned int wire,plane,tpc,cstat;

  double a,c;
 
  
  
  // loop on selected planes
 
  if(hitlist.size()==0) // this should never happen.
   {
   if(!startflag[nClust]){
    fWireVertex.push_back(wire_start);
    fTimeVertex.push_back(time_start);
    }
   if(!endflag){
    fWireLast.push_back(wire_end);
    fTimeLast.push_back(time_end);
    } 
     
     
   }
  
  GetPlaneAndTPC((*hitlist.begin()),plane,cstat,tpc,wire);
   

  //get paramters of the straight line fit. (and rescale them to cm/cm)
   a=lineslope[nClust]*fWireTimetoCmCm;
   c=lineinterc[nClust]*fTimetoCm;
 
  //get slope of lines orthogonal to those found crossing the shower.
  double aprim=0;
  if(a){
    aprim=-1./a;
  }

  //std::cout << "========= line params, nClust, plane: a,c " << nClust << " " << nClust << " " << a << " " << lineslope[nClust] << " " << c << std::endl;

  // find extreme intercepts. For the time being we don't care which one is the start point and which one is the endpoint.
  
  double extreme_intercept_high,extreme_intercept_low;
  
  Find_Extreme_Intercepts(hitlist,aprim,extreme_intercept_high,extreme_intercept_low);
    
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
  
  double wire_online_end,wire_online_begin,time_online_end,time_online_begin;
  
  gser.GetPointOnLineWSlopes(a,c,extreme_intercept_end,wire_online_end,time_online_end);
  gser.GetPointOnLineWSlopes(a,c,extreme_intercept_start,wire_online_begin,time_online_begin);
  
 // std::cout << " :::::::: wire_online begin point " << wire_online_begin << " " << time_online_begin << std::endl;

  //calculate the first and last cluster points on the through line:
  
  art::Ptr<recob::Hit> startHit=FindClosestHit(hitlist, (unsigned int)wire_online_begin,time_online_begin);
  art::Ptr<recob::Hit> endHit=FindClosestHit(hitlist, (unsigned int)wire_online_end,time_online_end);
  
  GetPlaneAndTPC(startHit,plane,cstat,tpc,wire);
  wire_start=wire;
  time_start=startHit->PeakTime();
  
  GetPlaneAndTPC(endHit,plane,cstat,tpc,wire);
  wire_end=wire;
  time_end=endHit->PeakTime();
  
 
  CalculateAxisParameters(nClust,hitlist,wire_start,time_start,wire_end,time_end);
  
   // \todo save rough start hits here. 
   // andthen refined as a separate entry in TTree.
  
///////////////// Clean up done up to here.
// Refine start points
//
  RefineStartPoints(nClust,hitlist, wire_start,time_start);

  std::cout << " ----++++ determined start wire points Cluster " << nClust << " "  << wire_start << " " << time_start <<std::endl;  



    
    if(!startflag[nClust]){
    fWireVertex.push_back(wire_start);
    fTimeVertex.push_back(time_start);
    }
    if(!endflag){
    fWireLast.push_back(wire_end);
    fTimeLast.push_back(time_end);
    }
   
   
}

  
  
  
// ------------------ select the clusters for the best two planes. This should make sure it's looking at one TPC in the future.


void cluster::ShowerAngleCluster::Find2DBestPlanes(std::vector<int> &best_planes){
   
   
  // first loop: select planes where the shower has a large enough correlation factor
  for (unsigned int iplane=0;iplane<fNPlanes;iplane++){
    if(fabs(tgx[iplane]->GetCorrelationFactor()) > 0.6 )	
      best_planes.push_back(iplane);
  //  std::cout << " correlation factors in 0.6 search" << tgx[iplane]->GetCorrelationFactor() << std::endl;
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

   // std::cout << "pushing back " << maxplane << std::endl;
    used_plane=maxplane;  // to cut out redundancy.
    best_planes.push_back(maxplane);

    } // end while loop

  /////test values:
  //for(unsigned int ii=0;ii<best_planes.size();ii++)
   // std::cout << "======+++==== " << best_planes[ii] << " " << 
    //tgx[best_planes[ii]]->GetCorrelationFactor() << std::endl; 


   
  return;
}



void cluster::ShowerAngleCluster::RefineStartPoints(unsigned int nClust,std::vector< art::Ptr < recob::Hit> > hitlist, double  wire_start,double time_start)
{
  //parameters to be set later?
  double linearlimit=12; 
  double ortlimit=15;
  double difflimit=1;
  
  //first select a subset of points that is close to the selected start point:
  std::vector < art::Ptr<recob::Hit> > hitlistlocal;
  std::vector < art::Ptr<recob::Hit> > hitlistlocal_refined;  //these in principle should only be hits on the line.
  
 for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    	art::Ptr<recob::Hit> theHit = (*hitIter);
    	double time = theHit->PeakTime() ;  
    	unsigned int plane,cstat,tpc,wire;
	GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
	
	double wonline,tonline;
	gser.GetPointOnLine(lineslope[nClust],lineinterc[nClust],wire,time,wonline,tonline);
	
	//calculate linear distance from start point and orthogonal distance from axis
	double lindist=gser.Get2DDistance(wire,time,wire_start,time_start);
	double ortdist=gser.Get2DDistance(wire,time,wonline,tonline);
	
	//std::cout << bestplanes[xx] << " " << lineslope[bestplanes[xx]] << " " << lineinterc[bestplanes[xx]] << " w,t: " << wire << " " << time << " ws,ts " << wire_start[bestplanes[xx]]<< " " << time_start[bestplanes[xx]] <<" "<< lindist << " " << ortdist << std::endl;
	
	if(lindist<linearlimit && ortdist<ortlimit)
	  hitlistlocal.push_back(theHit);
    
    
    }
   // std::cout << "hitlist_local size " << nClust << " " <<  hitlistlocal.size() << std::endl;
   // end old loop on best planes 
    
    
    
  //ok. now I have a sub selection of hits that are close to the perceived start of the shower.
 
    //std::cout << "hitlist_local size " << 
 //nwire_start.resize(fNPlanes);
 // ntime_start.resize(fNPlanes);
   int  nwire_start(fNPlanes);
   double ntime_start(fNPlanes);

  double distdiff=0;
  art::Ptr<recob::Hit> newFirstHit;
 
 double mindistdiff=999999;
 
 for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlistlocal.begin(); hitIter != hitlistlocal.end();  hitIter++){
    	art::Ptr<recob::Hit> theHit = (*hitIter);
    	double time = theHit->PeakTime() ;  
    	unsigned int plane,cstat,tpc,wire;
	GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
        distdiff=0;
	// First determine the local angle from the hit.
	double angle=Get2DAngleForHit(theHit,hitlist);
	
	//std::cout << "angle " << angle << " " << cos(angle) << " " << cos(TMath::Pi()*angle/180) << std::endl;
	//loop over all the hits from the same plane and local rectangle and calculate the distance
	
	// maybe recalculate the local rectangle for each hit?!
	for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIterin = hitlistlocal.begin(); hitIterin != hitlistlocal.end();  hitIterin++){
	      if(hitIterin==hitIter)
		  continue;
	      
	art::Ptr<recob::Hit> theHitin = (*hitIterin);
	double intime = theHitin->PeakTime() ;  
	unsigned int inplane,inwire;
	GetPlaneAndTPC(theHitin,inplane,cstat,tpc,inwire);
	    
	// then calculate the distances resulting from omega: dist = slope*2Dwirepitch and dist=dt^2+dw^2
	double stdist=gser.Get2DDistance(inwire,intime,wire,time);
	    
	double pitchdist=gser.Get2DPitchDistance(angle,inwire,wire);
	    
	distdiff+=fabs(pitchdist-stdist);
	     
	    // if(wire==wire_start[bestplanes[xx]])
	
	     
	  }
	
	if(distdiff<mindistdiff){
	  mindistdiff=distdiff;
	  newFirstHit=theHit;
	  nwire_start=wire;
	  ntime_start=time;
	  }
	 
	 // std::cout <<"plane:" << plane <<  " HIT (w,t) "<<wire <<","<<time<< " minimum hit difference" << mindistdiff << " ddiff " << distdiff << " scaled limit: " << distdiff/hitlistlocal[bestplanes[xx]].size()  << std::endl;
	 
	 
	 if(distdiff/hitlistlocal.size()<difflimit)
	    hitlistlocal_refined.push_back(theHit);
	// set the points with minimum difference into new collection. 
	
  
    }
   
   mf::LogVerbatim("ShowerAngleCluster") << "new wire and time starts, before last step (TEST): " << nClust << " " << nwire_start << " " << ntime_start << " old vals: " << wire_start<< " " << time_start << std::endl;
      
   // here was the end of the second bsetplanes loop.
  //

  //loop on refined hit collection to find points without holes - on the 
 
  
 
   int counter=0;
   hitlistlocal_refined.size();
  // hitlistlocal_refined.sort(cluster::SortByWire());   // later do a sort by pitch, need slope  and direction?
   
   //   std::cout << " counter " <<  counter <<" " << hitlistlocal_refined.size() << " nClust: " <<nClust << " "  << std::endl;
      int lsublclustsize=0;
      int sizecounter=0;
      
//       // going down the list??
//    if(hitlistlocal_refined.size()>1) {
//       for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlistlocal_refined.end()-1; hitIter != hitlistlocal_refined.begin();  hitIter--){
// 	
//       art::Ptr<recob::Hit> theHit = (*(hitIter));
//      //	double time = theHit->PeakTime() ;  
//       unsigned int plane,cstat,tpc,wire;
// 	GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
// 	//std::cout << "final hits " << plane << " " << wire << " " << time << std::endl;
// 	}
//        }
       
    if(hitlistlocal_refined.size()>1) {
	for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlistlocal_refined.end()-1; hitIter != hitlistlocal_refined.begin();  hitIter--){
	  //std::cout << " counter " <<  counter <<" " << std::endl;
	  art::Ptr<recob::Hit> theHit = (*(hitIter));
	  double time = theHit->PeakTime() ;  
	  unsigned int plane,cstat,tpc,wire;
	  GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
  
	  unsigned int plane2,cstat2,tpc2,wire2;
	  art::Ptr<recob::Hit> HitBefore = (*(hitIter-1));
	  GetPlaneAndTPC(HitBefore,plane2,cstat2,tpc2,wire2);
	
	 // std::cout << "wires " << wire2 << " " << wire << " " << fabs((double)wire2-(double)wire) << std::endl;
	
	
	  sizecounter++;
	  if(sizecounter>=lsublclustsize)
	    {
	    
	    nwire_start=wire;  
	    ntime_start=time;
	 //   std::cout << "sizecounter, subclustize" << lsublclustsize << " " << sizecounter << " " << wire << " " << time << std::endl;
	    lsublclustsize=sizecounter;  
	    }
	  
	  if((fabs((double)wire2-(double)wire)>1)) //presumably end of shower start, and following are comptons somewhere.
	  {
	  sizecounter=0;
	  }
	  
	  counter++;
	
	//go from back and stop at hit that has no hole after wires.
	
	} // end for loop
      }
   // end of third bestplanes loop
   
//   for(unsigned int xx=0;xx<bestplanes.size();xx++){
//   std::cout << "new wire and time starts: " << nClust << " " << nwire_start << " " << ntime_start << std::endl;
   
   //}
   

}



void cluster::ShowerAngleCluster::Find_Extreme_Intercepts(std::vector < art::Ptr<recob::Hit> > hitlist,
			     double perpslope,
			     double &inter_high,
			     double &inter_low)  
{
  
  inter_high=-999999;
  inter_low=999999;

  unsigned int plane,tpc,wire,cstat;
  
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
    
    //wire_bar+=wire;
    //time_bar+=time;	
    //nhits++;
    
    double intercept=time*fTimetoCm-perpslope*(double)wire*fWiretoCm;
    
    if(intercept > inter_high ){
      inter_high=intercept;
    }
    if(intercept < inter_low ){
      inter_low=intercept;
    }  

    

  }   // end of first HitIter loop, at this point we should have the extreme intercepts 

}




art::Ptr<recob::Hit> cluster::ShowerAngleCluster::FindClosestHit(std::vector < art::Ptr < recob::Hit > > hitlist,
			     unsigned int wire_online,
			     double time_online)
{
  
  double min_length_from_start=99999;
  art::Ptr<recob::Hit> nearHit;
   
  unsigned int plane,tpc,wire,cstat;
   
   
  for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
    art::Ptr<recob::Hit> theHit = (*hitIter);
    double time = theHit->PeakTime() ;  
    GetPlaneAndTPC(theHit,plane,cstat,tpc,wire);
    
    double dist_mod=gser.Get2DDistance(wire_online,time_online,wire,time);
    //TMath::Sqrt( pow(((double)wire_online-(double)wire*fWirePitch),2)+pow((time_online-time*fDriftVelocity*fTimeTick),2) );	

    if(dist_mod<min_length_from_start){
	//wire_start[plane]=wire;
	//time_start[plane]=time;
	nearHit=(*hitIter);
	min_length_from_start=dist_mod;
	}	

  } 
  
return nearHit;    
}

namespace cluster {

DEFINE_ART_MODULE(ShowerAngleCluster);

}

//#endif // SHOWERANGLECLUSTER_H


