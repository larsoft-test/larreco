////////////////////////////////////////////////////////////////////////
/// \file  ShowerAngleCluster_module.cc
/// \brief Create a Cluster with an angle 
///
/// \version $Id: AngleCluster.cxx,v 0.1 19/07/2011 12:45:16 PM  andrzejs $
/// \author andrzej.szelc@yale.edu
///  
////////////////////////////////////////////////////////////////////////
// This class solves the following problem:
//
// Create a cluster saving its slope angle

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

#ifndef SHOWERANGLECLUSTER_H
#define SHOWERANGLECLUSTER_H

#include "art/Framework/Core/EDProducer.h" // include the proper bit of the framework
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include <vector>
#include <string>


#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "Utilities/LArProperties.h"
#include "Utilities/GeometryUtilities.h"
#include "Utilities/DetectorProperties.h"
//#include "RecoAlg/HoughBaseAlg.h"
#include "RecoAlg/ClusterParamsAlg.h"

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>


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

//#include "TMatrixD.h"
//#include "TVectorD.h"
//#include "TDecompSVD.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1F.h"



// LArSoft includes
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"
#include "SimulationBase/MCTruth.h"
#include "Geometry/geo.h"
#include "RecoBase/Cluster.h"
#include "Utilities/AssociationUtil.h"
#include "Geometry/PlaneGeo.h"

//#include "SimulationBase/simbase.h"
//#include "RawData/RawDigit.h"
//#include "SummaryData/summary.h"
#include "CLHEP/Random/JamesRandom.h"
#include "Utilities/SeedCreator.h"




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
   
      

  private:

    double Get2DAngleForHit( art::Ptr<recob::Hit> starthit,std::vector < art::Ptr < recob::Hit> > hitlist);
    double Get2DAngleForHit( unsigned int wire, double time,std::vector < art::Ptr < recob::Hit> > hitlist);
    void ClearandResizeVectors(unsigned int nClusters);

 
    //HoughBaseAlg fHBAlg; 
    ClusterParamsAlg fCParAlg;
 //   double fWiretoCm,fTimetoCm,fWireTimetoCmCm;
    
    std::vector< unsigned int >fNWires;
    double fNTimes;
    
    
    
    
    
    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detp;
    art::ServiceHandle<util::LArProperties> larp;
    art::ServiceHandle<geo::Geometry> geo;
    util::GeometryUtilities gser;
   
    float fTimeTick; // time sample in us
    float fPresamplings;
    float fDriftVelocity;
     
    int fRun,fEvent,fSubRun;
    bool fForceRightGoing;

    //input parameter labels:
 
    std::string fClusterModuleLabel;
  
    std::vector <double> xangle;       // in cm, cm
   
    std::vector <double> lineslopetest;   // in wire, time
    std::vector <double>  lineinterctest;   
    std::vector<double>  fWireVertex,fTimeVertex;
    std::vector<double>  fWireEnd,fTimeEnd;
    std::vector<double> fVerticalness;

    unsigned int fNPlanes; // number of planes  

    TH1F *  fh_omega_single;
    TTree* ftree_cluster;
    std::vector<bool> startflag;
    bool endflag;
    bool matchflag;
  
  }; // class ShowerAngleCluster

}

#endif // SHOWERANGLECLUSTER_H


////////////////////////////////////////////////// Module definition

// ***************** //

//------------------------------------------------------------------------------
cluster::ShowerAngleCluster::ShowerAngleCluster(fhicl::ParameterSet const& pset)
 :// fHBAlg(pset.get< fhicl::ParameterSet >("HoughBaseAlg")),
   fCParAlg(pset.get< fhicl::ParameterSet >("ClusterParamsAlg"),pset.get< std::string >("module_type"))
{
  this->reconfigure(pset);
  produces< std::vector<recob::Cluster> >();
  produces< art::Assns<recob::Cluster, recob::Hit>  >(); 
 // produces< art::Assns<recob::Cluster, recob::Hit>  >(); 
  //produces< art::Assns<recob::Cluster, recob::Cluster>  >(); 
  produces< std::vector < art::PtrVector <recob::Cluster> >  >();
  
    // Create random number engine needed for PPHT
  createEngine(SeedCreator::CreateRandomNumberSeed(),"HepJamesRandom");

  
}


void cluster::ShowerAngleCluster::reconfigure(fhicl::ParameterSet const& pset) 
{
  fClusterModuleLabel 		=pset.get< std::string >("ClusterModuleLabel");
  fCParAlg.reconfigure(pset.get< fhicl::ParameterSet >("ClusterParamsAlg"));
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




//-----------------------------------------------
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

    
  // this will not change on a run per run basis.
  fNPlanes = geo->Nplanes();
 // fWirePitch = geo->WirePitch(0,1,0);    //wire pitch in cm
  fTimeTick=detp->SamplingRate()/1000.; 
 
  /**Get TFileService and define output Histograms*/
  art::ServiceHandle<art::TFileService> tfs;

  
   // fNTimes=geo->DetHalfWidth(tpc)*2/(fTimetoCm);
    fNWires.resize(fNPlanes);
    
  for(unsigned int i=0;i<fNPlanes;++i){
   
     fNWires[i]=geo->Nwires(i); //Plane(i,tpc).Nwires();
  }
  

  
  
  
  
  ftree_cluster =tfs->make<TTree>("ShowerAngleCluster","Results");/**All-knowing tree with reconstruction information*/
   fh_omega_single= tfs->make<TH1F>("fh_omega_single","Theta distribution Hit",720,-180., 180.) ;
    
    ftree_cluster->Branch("run",&fRun,"run/I");
    ftree_cluster->Branch("subrun",&fSubRun,"subrun/I");
    ftree_cluster->Branch("event",&fEvent,"event/I");
    ftree_cluster->Branch("nplanes",&fNPlanes,"nplanes/I");
  
    ftree_cluster->Branch("wire_vertex","std::vector<double>", &fWireVertex);
    ftree_cluster->Branch("time_vertex","std::vector<double>", &fTimeVertex);
      
    ftree_cluster->Branch("wire_last","std::vector<double>", &fWireEnd);
    ftree_cluster->Branch("time_last","std::vector<double>", &fTimeEnd);
    
    
}

// ************************************* //
void cluster::ShowerAngleCluster::ClearandResizeVectors(unsigned int nClusters) {
  
  ///////////////

  
   fVerticalness.clear();
    
  
  
   startflag.clear();
   lineslopetest.clear();
   lineinterctest.clear();
  
   xangle.clear();
   fWireVertex.clear();
   fTimeVertex.clear();
   fWireEnd.clear();
   fTimeEnd.clear();
   xangle.resize(nClusters); 


    fWireVertex.resize(nClusters);
    fTimeVertex.resize(nClusters);
    fWireEnd.resize(nClusters);
    fTimeEnd.resize(nClusters);
    
    fVerticalness.resize(nClusters);
    lineslopetest.resize(nClusters);
    lineinterctest.resize(nClusters);

    
   
    
  
  
}
  


  
  
  
// ***************** //
void cluster::ShowerAngleCluster::produce(art::Event& evt)
{ 
  

  
  
 mf::LogWarning("ShowerAngleCluster") << "In produce module " ; 

  //%%%%% this goes into ana module.
  //Find run, subrun and event number:
  fRun = evt.id().run();
  fSubRun = evt.id().subRun();
  fEvent = evt.id().event();

  
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
 
  
  
  for(unsigned int iClust = 0; iClust < clusterListHandle->size(); iClust++){

    art::Ptr<recob::Cluster> cl(clusterListHandle, iClust);
   std::vector< art::Ptr<recob::Hit> > hitlist = fmh.at(iClust);
    
  
    
    ///////////////////////////////////

      clusters.push_back(cl);

      std::vector< double > spos=cl->StartPos();
      std::vector< double > sposerr=cl->SigmaStartPos();
      
      std::vector< double > epos=cl->EndPos();
      std::vector< double > eposerr=cl->SigmaEndPos();
      
     // Start positions are determined elsewhere and accepted here
//       if(spos[0]!=0 && spos[1]!=0 && sposerr[0]==0 && sposerr[1]==0 ){
// 	fWireVertex.push_back(spos[0]);
// 	fTimeVertex.push_back(spos[1]);
// 	startflag.push_back(true);
// 	//std::cout << "setting external starting points " << spos[0] << " " << spos[1] <<" " << sposerr[0] <<" "<< sposerr[1] << std::endl;
//       }
// 	  else
// 	startflag.push_back(false); 
	
      ///! Change for accepting DBCluster and cheatcluster, so that it doesn't get fooled.
      startflag.push_back(false); 
      
      
      if(epos[0]!=0 && epos[1]!=0 && eposerr[0]==0 && eposerr[1]==0 ){
	fWireEnd.push_back(epos[0]);
	fTimeEnd.push_back(epos[1]);
	endflag=true;
	//std::cout << "setting external ending points " << epos[0] << " " << epos[1] << std::endl;
      }
	
    std::cout << " hitlist size: " << hitlist.size() << std::endl;
    if(hitlist.size()<=15 )
	continue;

      
    double lineslope, lineintercept,goodness,wire_start,time_start,wire_end,time_end;
    int nofshowerclusters=0;
   // std::cout << "++++ hitlist size " << hitlist.size() << std::endl;
    
    //////////////////////////////////
   // fCParAlg.Find2DAxisRoughHighCharge(lineslope,lineintercept,goodness,hitlist);
   //  std::cout << "%%%%%%%% lineslope, intercept " << lineslope << " "<< lineintercept << std::endl;

    
    
    fCParAlg.Find2DAxisRough(lineslope,lineintercept,goodness,hitlist);
    fVerticalness[iClust]=goodness;
    //if(hitlist_high.size()<=3 )
	//continue;
    
    fCParAlg.Find2DStartPointsHighCharge( hitlist,wire_start,time_start,wire_end,time_end);
    
    fWireVertex[iClust]=wire_start;
    fTimeVertex[iClust]=time_start;
   

    
    
    double wstn=0,tstn=0,wendn=0,tendn=0;
    fCParAlg.FindTrunk(hitlist,wstn,tstn,wendn,tendn,lineslope,lineintercept);
    int fDirection = (wstn<wendn)  ? 1 : -1 ;     // if current starting point is less then end point then direction is 1.
    
    
    double HiBin,LowBin,invHiBin,invLowBin;
    fCParAlg.FindDirectionWeights(lineslope,wstn,tstn,wendn,tendn,hitlist,HiBin,LowBin,invHiBin,invLowBin); 
    
    if(invHiBin+invLowBin> 1000)
      nofshowerclusters++;
    
    fDirection=fCParAlg.DecideClusterDirection(hitlist,lineslope,wstn,tstn,wendn,tendn);
     std::cout << "%%%%%%%% direction start points: (" << wstn<<","<<tstn<<"), end: ( "<<wendn << ","<<tendn <<")" << "Direction: " << fDirection << std::endl;
    wire_start=wstn;
    time_start=tstn;
    wire_end=wendn;
    time_end=tendn;
    fCParAlg.RefineStartPointsHough(hitlist, wire_start,time_start,wire_end,time_end,fDirection); 
     std::cout << "%%%%%%%% Hough line refine start points: (" << wire_start<<","<<time_start<<"), end: ( "<<wire_end << ","<<time_end <<")" << "Direction: " << fDirection << std::endl; 
    fWireVertex[iClust]=wire_start;
    fTimeVertex[iClust]=time_start;
    fWireEnd[iClust]=wire_end;
    fTimeEnd[iClust]=time_end; 
    lineslopetest[iClust]=lineslope; 
    lineinterctest[iClust]=lineintercept;
      
	  
     xangle[iClust]=Get2DAngleForHit( fWireVertex[iClust],fTimeVertex[iClust], hitlist);
    } // End loop on clusters.
  

   

  ////ugly temp fix
  std::vector< double > errors;
  errors.resize(clusterListHandle->size());
  
  for(unsigned int i=0;i<clusterListHandle->size();i++)
  {
   if(fVerticalness[i]<1) 
    errors[i]=0.1;
   else
    errors[i]=10; 
  }
 
  // make an art::PtrVector of the clusters
  std::unique_ptr<std::vector<recob::Cluster> > ShowerAngleCluster(new std::vector<recob::Cluster>);
  std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);
  
  for(unsigned int iplane=0;iplane<clusterListHandle->size();iplane++){
  //for(unsigned int iplane=0;iplane<fNPlanes;iplane++){
    std::cout << " in save loop, \"iplane\" " << iplane <<  std::endl;
    std::vector< art::Ptr<recob::Hit> > hitlist = fmh.at(iplane);
    std::cout << " hitlist size: " << hitlist.size() <<  std::endl;
    if(hitlist.size()>15) {
      
//     double wverror=fWireVertex[iplane]*0.05,tverror=fTimeVertex[iplane]*0.05;
//     
//     if(startflag[iplane])
//     {
//      wverror=0;
//      tverror=0;
//     }
      
    unsigned int xplane;
    xplane=hitlist[0]->WireID().Plane;
    geo::View_t viewfix = hitlist[0]->View();
    
    std::cout << " houghwire " << fWireVertex[iplane] << std::endl;
     std::cout << " errors " << errors[iplane] << std::endl;
     std::cout << " houghtime " << fTimeVertex[iplane] << std::endl;
     std::cout << " wirelast " << fWireEnd[iplane] << std::endl;
     std::cout << " timelast " <<  fTimeEnd[iplane] << std::endl; 
     std::cout << " xangle " <<  xangle[iplane] << std::endl; 
     std::cout <<"lineslope " <<  lineslopetest[iplane] << std::endl; 
       std::cout <<"lineinterc " << lineinterctest[iplane] <<std::endl;		  
    recob::Cluster temp(fWireVertex[iplane], errors[iplane],
                         fTimeVertex[iplane], errors[iplane],
                         fWireEnd[iplane], fWireEnd[iplane]*0.05,
                         fTimeEnd[iplane], fTimeEnd[iplane]*0.05,  
		    xangle[iplane], xangle[iplane]*0.05, lineslopetest[iplane],lineinterctest[iplane],5.,
		    viewfix,
		    xplane);
    
      std::cout << " Saving cluster for plane: " << iplane << " w,t " << fWireVertex[iplane] << " " << fTimeVertex[iplane] << " 2D angle: " <<  xangle[iplane] << " lslope " << lineslopetest[iplane] << std::endl;

    
    ShowerAngleCluster->push_back(temp);
  
    util::CreateAssn(*this, evt, *(ShowerAngleCluster.get()), hitlist, *(assn.get()));
    }

    
  }

  /////////////////////////////////////////////
  std::unique_ptr< std::vector < art::PtrVector < recob::Cluster > > > classn(new std::vector < art::PtrVector < recob::Cluster > >);
  
  if(!matchflag)
     {
      //matching code here.
      //// declare,  matching table
      
    //  for(unsigned int ii=0;ii<clusterListHandle->size();ii++){
	//for(unsigned int jj=ii+1;jj<clusterListHandle->size();jj++){
	//  // don't try to match clusters in the same view.
	//  if(ShowerAngleCluster[ii]->View()==ShowerAngleCluster[jj]->View())
	 //   continue;
	  
	 // std::vector<double > xpos1=ShowerAngleCluster[ii]->StartPos();
	 // std::vector<double > xpos2=ShowerAngleCluster[jj]->StartPos();
	
	  
	//} //end jj loop
     // } // end ii loop
       
       
       
       
       //
       
    art::PtrVector < recob::Cluster > cvec;
	
    for(unsigned int ip=0;ip<fNPlanes;ip++)  {
	art::ProductID aid = this->getProductID< std::vector < recob::Cluster > >(evt);
	art::Ptr< recob::Cluster > aptr(aid, ip, evt.productGetter(aid));
	cvec.push_back(aptr);
      }
     if((*ShowerAngleCluster).size()>0)   //make sure to make associations to something that exists.
      classn->push_back(cvec); 
     }
  else
    {
      std::cout << " the else: " << clusterAssociationHandle->size() << " " << (*ShowerAngleCluster).size() << " condition: " << 
      (0<clusterAssociationHandle->size() && 0<(*ShowerAngleCluster).size()) << std::endl;
    for(unsigned int i=0;i<clusterAssociationHandle->size() && i<(*ShowerAngleCluster).size();i++)   // do not save associations if clusters are not saved
      classn->push_back((*clusterAssociationHandle)[i]);
    }

  /**Fill the output tree with all information */
  ftree_cluster->Fill();

  evt.put(std::move(ShowerAngleCluster));
  evt.put(std::move(assn));
  evt.put(std::move(classn));

  
}



////////////////////////////////////////////////////////////////////////////////
// Method to get the 2D angle ogf a Cluster based on its starting wire and time.
////////////////////////////////////////////////////////////////////////////////

double cluster::ShowerAngleCluster::Get2DAngleForHit( unsigned int swire,double stime,std::vector < art::Ptr < recob::Hit> > hitlist) {
  
  fh_omega_single->Reset();
  
  unsigned int wire;
  // this should changed on the loop on the cluster of the shower
   for(std::vector < art::Ptr < recob::Hit > >::const_iterator hitIter = hitlist.begin(); hitIter != hitlist.end();  hitIter++){
     // art::Ptr<recob::Hit> theHit = (*hitIter);
      double time = (*hitIter)->PeakTime();  
      wire=(*hitIter)->WireID().Wire; 
      double omx=gser.Get2Dangle((double)wire,(double)swire,time,stime);
      fh_omega_single->Fill(180*omx/TMath::Pi(),(*hitIter)->Charge());
     }
    
  double omega = fh_omega_single->GetBinCenter(fh_omega_single->GetMaximumBin());// Mean value of the fit
   
  return omega; // in degrees.
}






namespace cluster {

  DEFINE_ART_MODULE(ShowerAngleCluster)

}
