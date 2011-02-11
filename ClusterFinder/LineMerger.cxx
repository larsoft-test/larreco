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
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//LArSoft includes:
#include "ClusterFinder/LineMerger.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>


namespace cluster{

  //-------------------------------------------------
  LineMerger::LineMerger(fhicl::ParameterSet const& pset) : 
    fClusterModuleLabel(pset.get<std::string>("ClusterModuleLabel")),
    fSlope             (pset.get<double     >("Slope")),
    fIntercept         (pset.get<double     >("Intercept"))
  {
    produces< std::vector<recob::Cluster> >();
  }

  //-------------------------------------------------
  LineMerger::~LineMerger()
  {
  }

  //-------------------------------------------------
  void LineMerger::beginJob()
  {
    //this doesn't do anything now, but it might someday
  }
    
  //------------------------------------------------------------------------------------//
  void LineMerger::produce(art::Event& evt)
  { 
    // Get a Handle for the input Cluster object(s).
    art::Handle< std::vector<recob::Cluster> > clusterVecHandle;
    evt.getByLabel(fClusterModuleLabel,clusterVecHandle);

    art::ServiceHandle<geo::Geometry> geo;
    int nplanes = geo->Nplanes();

    art::PtrVector<recob::Cluster> Cls[nplanes];//one PtrVector for each plane in the geometry
    std::vector<int> Cls_matches[nplanes];//vector with indicators for whether a cluster has been merged already

    // loop over the input Clusters
    for(unsigned int i = 0; i < clusterVecHandle->size(); ++i){
      
      //get a art::Ptr to each Cluster
      art::Ptr<recob::Cluster> cl(clusterVecHandle, i);
      
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
      for(art::PtrVectorItr<recob::Cluster> clusIter1 = Cls[i].begin(); clusIter1!=Cls[i].end();++clusIter1){
	art::Ptr<recob::Cluster> cl1 = *clusIter1;
	if(Cls_matches[i][clsnum1]==1){
	  clsnum1++;
	  continue;
	}
	SuperClusters->push_back(*cl1);
	Cls_matches[i][clsnum1]=1; 
	//SuperClusters->back().SetID(clustersfound);//IDs are sequential by plane, starting from 0
	++clustersfound;
	recob::Cluster SCl= SuperClusters->back();
	
	int clsnum2 = 0;
	for(art::PtrVectorItr<recob::Cluster> clusIter2 = Cls[i].begin(); clusIter2!=Cls[i].end();++clusIter2){
	  art::Ptr<recob::Cluster> cl2 = *clusIter2;
	  if(Cls_matches[i][clsnum2]==1){
	    clsnum2++;
	    continue;
	  }

	  //check that the slopes are the same
	  bool sameSlope = SlopeCompatibility(SCl.dTdW(),cl2->dTdW());

	  double sclint = SCl.StartPos()[1] - SCl.dTdW()*SCl.StartPos()[0];
	  double cl2int = cl2->StartPos()[1] - cl2->dTdW()*cl2->StartPos()[0];
	  bool sameIntercept = InterceptCompatibility(sclint, cl2int);

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

    mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    mf::LogVerbatim("Summary") << "LineMerger Summary:";
    for(int i = 0; i<SuperClusters->size(); ++i) mf::LogVerbatim("Summary") << SuperClusters->at(i) ;

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
