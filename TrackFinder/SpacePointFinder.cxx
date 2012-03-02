//
// Name: SpacePointFinder.cxx
//
// Purpose: Implementation file for module SpacePointFinder.
//
// Created: 15-Dec-2011  H. Greenlee
//

#include <vector>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TrackFinder/SpacePointFinder.h"
#include "TrackFinder/SpacePointService.h"
#include "Geometry/geo.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Prong.h"
#include "Utilities/AssociationUtil.h"

namespace trkf {

  SpacePointFinder::SpacePointFinder(const fhicl::ParameterSet& pset) :
    //
    // Purpose: Constructor.
    //
    // Arguments: pset - Module parameters.
    //
    fFilter(true),
    fMerge(false),
    fNumEvent(0),
    fNumProng2(0),
    fNumProng3(0)
  {
    reconfigure(pset);
    produces<std::vector<recob::Prong> >();
    produces<art::Assns<recob::Prong, recob::Cluster> >();

    // Report.

    mf::LogInfo("SpacePointFinder") 
      << "SpacePointFinder configured with the following parameters:\n"
      << "  ClusterModuleLabel = " << fClusterModuleLabel << "\n"
      << "  Filter = " << fFilter << "\n"
      << "  Merge = " << fMerge;
  }

  SpacePointFinder::~SpacePointFinder()
  //
  // Purpose: Destructor.
  //
  {}

  void SpacePointFinder::reconfigure(fhicl::ParameterSet const& pset)
  //
  // Purpose: Reconfigure method.
  //
  // Arguments: pset - Configuration parameters.
  //
  {
    fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
    fFilter = pset.get<bool>("Filter");
    fMerge = pset.get<bool>("Merge");
  }

  void SpacePointFinder::beginJob()
  {}

  void SpacePointFinder::produce(art::Event& evt)
  //
  // Purpose: Produce method.
  //
  // Arguments: event - Art event.
  //
  {
    ++fNumEvent;

    // Get Services.

    art::ServiceHandle<trkf::SpacePointService> sptsvc;
    art::ServiceHandle<geo::Geometry> geom;

    // Get clusters.

    art::Handle< std::vector<recob::Cluster> > clusterh;
    evt.getByLabel(fClusterModuleLabel, clusterh);

    // Make a double or triple loop over clusters in distinct views
    // (depending on minimum number of views configured in SpacePointService).

    if(clusterh.isValid()) {

      // Make a collection of prongs that will eventually be inserted into the event.

      std::auto_ptr<std::vector<recob::Prong> > prongs(new std::vector<recob::Prong>);
      std::auto_ptr< art::Assns<recob::Prong, recob::Cluster> > assn(new art::Assns<recob::Prong, recob::Cluster>);
    
      // Make a hit vector which will be used to store hits to be passed
      // to SpacePointService.

      art::PtrVector<recob::Hit> hits;      

      // Loop over first cluster.

      int nclus = clusterh->size();
      for(int iclus = 0; iclus < nclus; ++iclus) {
	art::Ptr<recob::Cluster> piclus(clusterh, iclus);
	geo::View_t iview = piclus->View();

	// Test first view.

	if((iview == geo::kU && sptsvc->enableU()) ||
	   (iview == geo::kV && sptsvc->enableV()) ||
	   (iview == geo::kW && sptsvc->enableW())) {

	  // Store hits from first view into hit vector.

	  art::PtrVector<recob::Hit> ihits = piclus->Hits();
	  unsigned int nihits = ihits.size();
	  hits.clear();
	  hits.reserve(nihits);
	  for(art::PtrVector<recob::Hit>::const_iterator i = ihits.begin();
	      i != ihits.end(); ++i)
	    hits.push_back(*i);
	  
	  // Loop over second cluster.

	  for(int jclus = 0; jclus < iclus; ++jclus) {
	    art::Ptr<recob::Cluster> pjclus(clusterh, jclus);
	    geo::View_t jview = pjclus->View();

	    // Test second view.

	    if(((jview == geo::kU && sptsvc->enableU()) ||
		(jview == geo::kV && sptsvc->enableV()) ||
		(jview == geo::kW && sptsvc->enableW()))
	       && jview != iview) {

	      // Store hits from second view into hit vector.

	      art::PtrVector<recob::Hit> jhits = pjclus->Hits();
	      unsigned int njhits = jhits.size();
	      assert(hits.size() >= nihits);
	      //hits.resize(nihits);
	      while(hits.size() > nihits)
		hits.pop_back();
	      assert(hits.size() == nihits);
	      hits.reserve(nihits + njhits);
	      for(art::PtrVector<recob::Hit>::const_iterator j = jhits.begin();
		  j != jhits.end(); ++j)
		hits.push_back(*j);
	  
	      // If two-view space points are allowed, make them here.

	      if(sptsvc->minViews() <= 2) {
		std::vector<recob::SpacePoint> spts;
		sptsvc->makeSpacePoints(hits, spts,
					fFilter, fMerge, 0., 0.);

		// If we found some space points, make a prong.

		if(spts.size() > 0) {
		  art::PtrVector<recob::Cluster> clusters;
		  clusters.reserve(2);
		  clusters.push_back(piclus);
		  clusters.push_back(pjclus);

		  prongs->push_back(recob::Prong(clusters, spts));
		  util::CreateAssn(*this, evt, *(prongs.get()), clusters,*(assn.get()));
		  ++fNumProng2;
		}
	      }

	      // Loop over third cluster.

	      for(int kclus = 0; kclus < jclus; ++kclus) {
		art::Ptr<recob::Cluster> pkclus(clusterh, kclus);
		geo::View_t kview = pkclus->View();

		// Test third view.

		if(((kview == geo::kU && sptsvc->enableU()) ||
		    (kview == geo::kV && sptsvc->enableV()) ||
		    (kview == geo::kW && sptsvc->enableW()))
		   && kview != iview && kview != jview) {

		  // Store hits from third view into hit vector.

		  art::PtrVector<recob::Hit> khits = pkclus->Hits();
		  unsigned int nkhits = khits.size();
		  assert(hits.size() >= nihits + njhits);
		  //hits.resize(nihits + njhits);
		  while(hits.size() > nihits + njhits)
		    hits.pop_back();
		  assert(hits.size() == nihits + njhits);
		  hits.reserve(nihits + njhits + nkhits);
		  for(art::PtrVector<recob::Hit>::const_iterator k = khits.begin();
		      k != khits.end(); ++k)
		    hits.push_back(*k);

		  // Make three-view space points.

		  std::vector<recob::SpacePoint> spts;
		  sptsvc->makeSpacePoints(hits, spts,
					  fFilter, fMerge, 0., 0.);

		  // If we found some space points, make a prong.

		  if(spts.size() > 0) {
		    art::PtrVector<recob::Cluster> clusters;
		    clusters.reserve(3);
		    clusters.push_back(piclus);
		    clusters.push_back(pjclus);
		    clusters.push_back(pkclus);
		    prongs->push_back(recob::Prong(clusters, spts));
		    ++fNumProng3;
		  }
		}
	      }
	    }
	  }
	}
      }

      // Add prongs to event.

      evt.put(prongs);
      evt.put(assn);
    }
  }

  void SpacePointFinder::endJob()
  //
  // Purpose: Print summary.
  //
  {
    mf::LogInfo("SpacePointFinder") 
      << "SpacePointFinder statistics:\n"
      << "  Number of events = " << fNumEvent << "\n"
      << "  Number of 2-view prongs created = " << fNumProng2 << "\n"
      << "  Number of 3-view prongs created = " << fNumProng3;
  }
}
