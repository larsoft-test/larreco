////////////////////////////////////////////////////////////////////////
//
// LineMerger class
//
// maddalena.antonello@lngs.infn.it
// ornella.palamara@lngs.infn.it
// biagio.rossi@lhep.unibe.ch
// msoderbe@syr.edu
//
// This algorithm is designed to merge 2D lines with similar slope and intercept 
//  
////////////////////////////////////////////////////////////////////////

//Framework includes:
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//LArSoft includes:
#include "ClusterFinder/LineMerger.h"
#include "Geometry/geo.h"
#include "Geometry/WireGeo.h"
#include "Geometry/VolumeUtility.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>


namespace cluster{

  //-------------------------------------------------
  LineMerger::LineMerger(edm::ParameterSet const& pset) : 
    fClusterModuleLabel(pset.getParameter<std::string>("ClusterModuleLabel")),
    fSlope(pset.getParameter<double>("Slope")),
    fIntercept(pset.getParameter<double>("Intercept"))
  {
    produces< std::vector<recob::Cluster> >();
  }

  //-------------------------------------------------
  LineMerger::~LineMerger()
  {
  }

  //-------------------------------------------------
  void LineMerger::beginJob(edm::EventSetup const&)
  {
    //this doesn't do anything now, but it might someday
  }
    
  //------------------------------------------------------------------------------------//
  void LineMerger::produce(edm::Event& evt, edm::EventSetup const&)
  { 
    // Get a Handle for the input Cluster object(s).
    edm::Handle< std::vector<recob::Cluster> > clusterVecHandle;
    evt.getByLabel(fClusterModuleLabel,clusterVecHandle);

    edm::PtrVector<recob::Cluster> kU_Cls;
    std::vector<int> kU_Cls_matches;
    edm::PtrVector<recob::Cluster> kV_Cls;
    std::vector<int> kV_Cls_matches;

    // loop over the input Clusters
    for(unsigned int i = 0; i < clusterVecHandle->size(); ++i){
      
      //get a edm::Ptr to each Cluster
      edm::Ptr<recob::Cluster> cl(clusterVecHandle, i);
      //std::cout << *cl << std::endl;
      
      switch(cl->View()){
      case geo::kU :
	kU_Cls.push_back(cl);
	kU_Cls_matches.push_back(0);
	break;
      case geo::kV :
	kV_Cls.push_back(cl);
	kV_Cls_matches.push_back(0);
	break;
      default :
	break;
      }
    }

    //std::cout << kU_Cls.size() << " U Hough Clusters ; " << kV_Cls.size() << " V Hough Clusters " << std::endl;

     //////////////////////////////////////////////////////
    // Make a std::auto_ptr<> for the thing you want to put into the event
    // because that handles the memory management for you
    //////////////////////////////////////////////////////
    std::auto_ptr<std::vector<recob::Cluster> > SuperClusters(new std::vector<recob::Cluster>);

    int ikU_Cls=0;
    for(edm::PtrVectorItr<recob::Cluster> clusIter1 = kU_Cls.begin(); clusIter1!=kU_Cls.end();++clusIter1){
      edm::Ptr<recob::Cluster> cl1 = *clusIter1;
      if(kU_Cls_matches[ikU_Cls]==1){
	ikU_Cls++;
	continue;
      }
      SuperClusters->push_back(*cl1);
      kU_Cls_matches[ikU_Cls]=1; 
      SuperClusters->back().SetID(SuperClusters->size()-1);
      recob::Cluster SCl= SuperClusters->back();
  
      int jkU_Cls=0;
      for(edm::PtrVectorItr<recob::Cluster> clusIter2 = kU_Cls.begin(); clusIter2!=kU_Cls.end();++clusIter2){
	edm::Ptr<recob::Cluster> cl2 = *clusIter2;
	if(kU_Cls_matches[jkU_Cls]==1){
	  jkU_Cls++;
	  continue;
	}
	bool sameSlope = SlopeCompatibility(SCl.Slope(),cl2->Slope());
	bool sameIntercept = InterceptCompatibility(SCl.Intercept(),cl2->Intercept());
	bool AreCompatible = sameSlope && sameIntercept;
	if(AreCompatible){
	  SuperClusters->back() = SuperClusters->back() + *cl2;
	  kU_Cls_matches[jkU_Cls]=1;       
	}
	jkU_Cls++;
      }
      ikU_Cls++;
    }

    int ikV_Cls=0;
    for( edm::PtrVectorItr<recob::Cluster> clusIter1 = kV_Cls.begin(); clusIter1!=kV_Cls.end();++clusIter1){
      edm::Ptr<recob::Cluster> cl1 = *clusIter1;
      if(kV_Cls_matches[ikV_Cls]==1){
	ikV_Cls++;
	continue;
      }
      SuperClusters->push_back(*cl1);
      kV_Cls_matches[ikV_Cls]=1; 
      SuperClusters->back().SetID(SuperClusters->size()-1);
      recob::Cluster SCl=SuperClusters->back();
  
      int jkV_Cls=0;
      for( edm::PtrVectorItr<recob::Cluster> clusIter2 = kV_Cls.begin(); clusIter2!=kV_Cls.end();++clusIter2){
	edm::Ptr<recob::Cluster> cl2 = *clusIter2;
	if(kV_Cls_matches[jkV_Cls]==1){
	  jkV_Cls++;
	  continue;
	}
	bool sameSlope = SlopeCompatibility(SCl.Slope(),cl2->Slope());
	bool sameIntercept = InterceptCompatibility(SCl.Intercept(),cl2->Intercept());
	bool AreCompatible = sameSlope && sameIntercept;
	if(AreCompatible){
	  SuperClusters->back() = SuperClusters->back() + *cl2;
	  kV_Cls_matches[jkV_Cls]=1;       
	}
	jkV_Cls++;
      }
      ikV_Cls++;
    }

    std::sort(SuperClusters->begin(),SuperClusters->end());//sort before Putting
    std::cout<<"LineMerger found "<< SuperClusters->size() <<" Clusters : "<<std::endl;
    for(std::vector<recob::Cluster>::iterator clusIter1 = SuperClusters->begin(); clusIter1 !=  SuperClusters->end();  clusIter1++) {
      recob::Cluster cl1 = *clusIter1;
      std::cout << cl1 << std::endl;
    }

    evt.put(SuperClusters);
     
    return;

  }

  //------------------------------------------------------------------------------------//
  bool LineMerger::SlopeCompatibility(double slope1, double slope2)
  { 
    double sl1=atan(slope1);
    double sl2=atan(slope2);

    bool comp = fabs(sl1-sl2)<fSlope ? true : false;
    return comp;
  }
  //------------------------------------------------------------------------------------//
  bool LineMerger::InterceptCompatibility(double intercept1, double intercept2)
  { 
    double in1=intercept1;
    double in2=intercept2;
 
    bool comp = fabs(in1-in2)<fIntercept ? true : false;
    return comp;
  }

}//namespace cluster{
