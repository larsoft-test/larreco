////////////////////////////////////////////////////////////////////////
//
// \file SmallClusterFinder.cxx
//
// corey.adams@yale.edu
// andrzej.szelc@yale.edu 
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
	
	-Corey
*/
//
// 
///////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

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

// LArSoft includes
//#include "Simulation/sim.h"
//#include "Simulation/SimListUtils.h"
#include "Filters/SmallClusterFilter.h"
#include "Geometry/geo.h"
#include "Utilities/AssociationUtil.h"
//#include "SimulationBase/simbase.h"
//#include "RawData/RawDigit.h"
//#include "SummaryData/summary.h"
//#include "CLHEP/Random/JamesRandom.h"
//#include "Utilities/SeedCreator.h"

// ***************** //

//------------------------------------------------------------------------------
cluster::SmallClusterFilter::SmallClusterFilter(fhicl::ParameterSet const& pset)
{
	this->reconfigure(pset);
	produces< std::vector<recob::Cluster> >();				//This code makes clusters
	produces< art::Assns<recob::Cluster, recob::Hit>  >();  //Matches clusters with hits
  
	// Create random number engine needed for PPHT
 	// createEngine(SeedCreator::CreateRandomNumberSeed(),"HepJamesRandom");

  
}


void cluster::SmallClusterFilter::reconfigure(fhicl::ParameterSet const& pset) 
{
  	fHitFinderModuleLabel 	=pset.get< std::string 	> ("HitFinderModuleLabel"); 
  	fMaxHitsByPlane			=pset.get< std::vector<int> > ("MaxHitsByPlane");
  	fMaxTotalHits			=pset.get<	 int		> ("MaxTotalHits");
  	
  	//fRadiusSizePar is used to exclude hits from a cluster outside of a certain size
  	//fNHitsInClust ensures the clusters don't get too big
  	//max hits by plane filters this event if the hits on that plane is too big
  	//max hits total checks against the sum of all hits found.
 }

// ***************** //
cluster::SmallClusterFilter::~SmallClusterFilter()
{
	//Nothing to do in the destructor
}

//____________________________________________________________________________
 bool cluster::SmallClusterFilter::beginRun(art::Run& run)
  {
  	//nothing to do at beginRun()
    return true;
  }

//-----------------------------------------------

// ***************** //
void cluster::SmallClusterFilter::beginJob()
{
   	// this will not change on a run per run basis.
   	fNPlanes = geo->Nplanes(); 				//get the number of planes in the TPC
 	
   	/**Get TFileService and define output Histograms*/
    art::ServiceHandle<art::TFileService> tfs;

	return;
}



// ************************************* //
// Clear and resize - exactly what it sounds like?
// Don't know why it takes the number of clusters...
void cluster::SmallClusterFilter::ClearandResizeVectors(unsigned int nClusters) {

	  ///////////////
	//   fMinWire.clear();
	//   fMaxWire.clear();
	//   fMinTime.clear();
	//   fMaxTime.clear();
	//   
	//   fRMS_wire.clear();
	//   fRMS_time.clear();
	//   
	//   mcwirevertex.resize(fNPlanes);  // wire coordinate of vertex for each plane 
	//   mctimevertex.resize(fNPlanes);  // time coordinate of vertex for each plane
	//   
	//   fRMS_wire.resize(fNPlanes);
	//   fRMS_time.resize(fNPlanes); 
	//   
	  ////////////////
	hitlistbyplane.clear();
	hitlistbyplane.resize(fNPlanes);    
	return;
}
  
  
  
// ***************** //
// This method actually makes the clusters.
bool cluster::SmallClusterFilter::filter(art::Event& evt)
{ 
	//Check the size of the maxHitsByPlane vector against fNPlanes
	if (fMaxHitsByPlane.size() != fNPlanes) return false;

  /**Get Clusters*/
 

 	
 	//Get the hits for this event:
 	art::Handle< std::vector<recob::Hit> > HitListHandle;
	evt.getByLabel(fHitFinderModuleLabel,HitListHandle);  

	//A vector to hold hits, not yet filled:
	std::vector< art::Ptr<recob::Hit> > hitlist;
	
	//How many hits in this event?  Tell user:
	std::cout << " ++++ Hitsreceived received " << HitListHandle->size() << " +++++ " << std::endl;

  	//Catch the case were there are no hits in the event:
  	if(HitListHandle->size() ==0 )
  	{
    	std::cout << " no hits received! exiting " << std::endl;
    	return false;
  	}
  	if (HitListHandle->size() > (unsigned int) fMaxTotalHits)
  	{
  		std::cout << "Not an empty event, exiting." << std::endl;
  		return false;
  	}
  	
  	
  	
  	ClearandResizeVectors(HitListHandle->size());
  	// resizing once cluster size is known.
  
  	art::Ptr< recob::Hit> theHit;
  
	bool collFound = false;
	//add all of the hits to the hitlist, and sort them into hits by plane
  	for(unsigned int iHit = 0; iHit < HitListHandle->size(); iHit++){

    	theHit = art::Ptr< recob::Hit>(HitListHandle, iHit);
    	//std::cout << "Just read in this hit " << iHit << std::endl;


    	unsigned int p(0),w(0), t(0),cs(0); //c=channel, p=plane, w=wire, but what is t?
    	GetPlaneAndTPC(theHit,p,cs,t,w);  //Find out what plane this hit is on.
    
    	//Do a check to catch crazy hits:
    	if (theHit -> Charge() > 500 ) continue;
    
   		//add this hit to the total list
      	hitlist.push_back(theHit);

     	//add this hit to the list specific to this plane
     	hitlistbyplane[p].push_back(theHit);
     	//hitlistleftover[p].push_back(theHit);
 		//Just for searching for Ar39:

		if (theHit -> SignalType() == geo::kCollection) collFound = true;
      
    } // End loop on hits.
  	
  
  	
  	//Check against each plane:
  	for (unsigned int i = 0; i < fNPlanes; i ++){
  		if (hitlistbyplane[i].size() > (unsigned int) fMaxHitsByPlane[i]) return false;
  	}
  	
  	//check that there is at least 1 hit on collection:
  	if (!collFound) return false;
   	
   	std::cout << "\nPassing with " << hitlistbyplane[fNPlanes-1].size() << " hit(s) on collection.\n" << std::endl;
   	
   	return true;
}

// ******************************* //
// ******************************* //
int cluster::SmallClusterFilter::GetPlaneAndTPC(art::Ptr<recob::Hit> a, //the hit
						unsigned int &p, //plane
 						unsigned int &cs,  //cryostat
						unsigned int &t, //time
						unsigned int &w) //wire
{
	art::ServiceHandle<geo::Geometry> geom;
	unsigned int channel = a->Wire()->RawDigit()->Channel(); 
	geom->ChannelToWire(channel);
	p = a -> WireID().Plane;
 	t = a -> PeakTime();
	w = a -> WireID().Wire;
	return 0;
}

