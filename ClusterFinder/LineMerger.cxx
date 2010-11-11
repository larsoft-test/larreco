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

    edm::Service<geo::Geometry> geo;
    int nplanes = geo->Nplanes();

    edm::PtrVector<recob::Cluster> Cls[nplanes];//one PtrVector for each plane in the geometry
    std::vector<int> Cls_matches[nplanes];//vector with indicators for whether a cluster has been merged already

    // loop over the input Clusters
    for(unsigned int i = 0; i < clusterVecHandle->size(); ++i){
      
      //get a edm::Ptr to each Cluster
      edm::Ptr<recob::Cluster> cl(clusterVecHandle, i);
      
      switch(cl->View()){
      case geo::kU :
	Cls[0].push_back(cl);
	Cls_matches[0].push_back(0);
	break;
      case geo::kV :
	Cls[1].push_back(cl);
	Cls_matches[1].push_back(0);
	break;
      case geo::kW :
	Cls[2].push_back(cl);
	Cls_matches[2].push_back(0);
	break;
      default :
	break;
      }
    }

     //////////////////////////////////////////////////////
    // Make a std::auto_ptr<> for the thing you want to put into the event
    // because that handles the memory management for you
    //////////////////////////////////////////////////////
    std::auto_ptr<std::vector<recob::Cluster> > SuperClusters(new std::vector<recob::Cluster>);

    for(int i = 0; i<nplanes; ++i){
      int clustersfound = 0;//how many merged clusters found in each plane
      int clsnum1 = 0;
      for(edm::PtrVectorItr<recob::Cluster> clusIter1 = Cls[i].begin(); clusIter1!=Cls[i].end();++clusIter1){
	edm::Ptr<recob::Cluster> cl1 = *clusIter1;
	if(Cls_matches[i][clsnum1]==1){
	  clsnum1++;
	  continue;
	}
	SuperClusters->push_back(*cl1);
	Cls_matches[i][clsnum1]=1; 
	SuperClusters->back().SetID(clustersfound);//IDs are sequential by plane, starting from 0
	++clustersfound;
	recob::Cluster SCl= SuperClusters->back();
	
	int clsnum2 = 0;
	for(edm::PtrVectorItr<recob::Cluster> clusIter2 = Cls[i].begin(); clusIter2!=Cls[i].end();++clusIter2){
	  edm::Ptr<recob::Cluster> cl2 = *clusIter2;
	  if(Cls_matches[i][clsnum2]==1){
	    clsnum2++;
	    continue;
	  }
	  bool sameSlope = SlopeCompatibility(SCl.Slope(),cl2->Slope());
	  bool sameIntercept = InterceptCompatibility(SCl.Intercept(),cl2->Intercept());
	  if(sameSlope && sameIntercept){
	    SuperClusters->back() = SuperClusters->back() + *cl2;
	    Cls_matches[i][clsnum2]=1;       
	  }
	  clsnum2++;
	}
	clsnum1++;
      }
    }

    std::sort(SuperClusters->begin(),SuperClusters->end());//sort before Putting
    std::cout<<"LineMerger found "<< SuperClusters->size() <<" Clusters : "<<std::endl;
    for(std::vector<recob::Cluster>::iterator clusIter1 = SuperClusters->begin(); clusIter1 !=  SuperClusters->end();  clusIter1++) {
      recob::Cluster cl1 = *clusIter1;
      std::cout << cl1 << std::endl;
      //cl1.PrintHits();
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
