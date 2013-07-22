////////////////////////////////////////////////////////////////////////
//
// \file SmallClusterFinder.cxx
//
// andrzej.szelc@yale.edu 
// corey.adams@yale.edu
//
// This algorithm is designed to fid small clusters that could correspond to gama or low energy e
//
/*	There are two parameters that matter from the fcl file:
		fNHitsInClust is the number of hits that should be in these small clusters
			^-- Gamma events seem to rarely have more than 4 hits in the cluster
			^-- SN events are unclear.  Should this even be used for SN?
		fRadiusSizePar is the distance (in cm) between the small clusters and any other hits.
	
	This algorithm sorts the hits by plane, and then looks at each hit individually.  If
	there is a hit within RadiusSizePar, it gets added to a local list.  All other hits
	are ignored.  Then, if the number of hits that got added to the local list is greater 
	then NHitsInClust, the original hit is ignored.  If it's less, the original hit is 
	presumed to be part of a very small (or single hit) cluster.  So its added to the list
	of hits in the small cluster.
	
	This algorithm dumps all the hits that could be part of a small cluster into one cluster
	There is no grouping done between individual hits or anything like that.
	
	The scheme for cluster identification was:
	Clusters with ID 0 to nPlanes-1 are the gammas on planes 0, 1, ... nPlanes-1
	Clusters with ID nPlanes to 2*nPlanes -1 are the leftover hits on the planes 0, 1, ... nPlanes-1
	
	But has been updated to:
	Clusters are ID'd with numbers like 107, meaning the 7th cluster on plane 1.  The formula
	is ID = 100*iPlane + Cluster on that plane.  
	
	By Convention, cluster 000, 100, 200, etc. are the leftover hits that aren't gammas.
	
	-Corey
*/
//
// 
///////////////////////////////////////////////////////////////////////

#ifndef SMALLCLUSTERFILTER_H
#define SMALLCLUSTERFILTER_H

#include "art/Framework/Core/EDFilter.h" // include the proper bit of the framework
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include <vector>
#include <string>

#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/GeometryUtilities.h"
#include "Utilities/DetectorProperties.h"



//class recob::Hit;
// class TH1F;
// class TF1;
// class TTree;



namespace cluster {

  class SmallClusterFilter : public art::EDFilter {
    
  public:

    /**METHODS global*/
    explicit SmallClusterFilter(fhicl::ParameterSet const& pset);/**Constructor*/
    virtual ~SmallClusterFilter();                               /**Destructor*/
    void beginJob();                                     
    bool beginRun(art::Run& run);
    void reconfigure(fhicl::ParameterSet const& pset);
    bool filter(art::Event& evt);                       /**Routine that finds the cluster and sets the dTdW of the 2D shower*/


  private:
    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detp;
    art::ServiceHandle<util::LArProperties> larp;
    art::ServiceHandle<geo::Geometry> geo;
    util::GeometryUtilities gser;
   
    std::vector < std::vector< art::Ptr<recob::Hit> > > hitlistbyplane;  //list of all hits on each plane

    std::vector< unsigned int >fNWires;			// Number of wires on each plane
    
    //input parameters
    std::string fHitFinderModuleLabel; 
    std::vector<int> fMaxHitsByPlane;
    int fMaxTotalHits;
    							// Remember, this is the *small* cluster finder

    unsigned int fNPlanes; // number of planes  


    void ClearandResizeVectors(unsigned int nHits);
    
    int GetPlaneAndTPC(	art::Ptr<recob::Hit> a,
    					unsigned int &p,
    					unsigned int &cs,
    					unsigned int &t,
    					unsigned int &w);
    
    
  }; // class SmallAngleFinder

}

#endif // SHOWERANGLECLUSTER_H

