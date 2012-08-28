//
// Name: BezierTrackerModule.cxx
//
// Purpose: Implementation file for module BezierTrackerModule.
//
// Ben Jones, MIT
//

#include <vector>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TrackFinder/BezierTrackerModule.h"
#include "TrackFinder/SeedFinderAlgorithm.h"
#include "Geometry/geo.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "RecoBase/Hit.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Track.h"
#include "RecoBase/Prong.h"
#include "TrackFinder/BezierTrack.h"
#include "Utilities/AssociationUtil.h"
#include "TrackFinder/BezierTrackerAlgorithm.h"

namespace trkf {

  BezierTrackerModule::BezierTrackerModule(const fhicl::ParameterSet& pset)
  {
    reconfigure(pset);
    produces< std::vector<recob::Track> >();
    produces< art::Assns<recob::Track, recob::Hit> >();
    
  }

  BezierTrackerModule::~BezierTrackerModule()
  {
  }

  void BezierTrackerModule::reconfigure(fhicl::ParameterSet const& pset)
  {
    fSeedModuleLabel   = pset.get<std::string>("SeedModuleLabel");
    fClusterModuleLabel= pset.get<std::string>("ClusterModuleLabel");
    fHitModuleLabel    = pset.get<std::string>("HitModuleLabel");
    fTrackMode         = pset.get<double>("TrackMode");
    fMakeHitAssns      = pset.get<bool>("MakeHitAssns");
    
    fBTrackAlg = new BezierTrackerAlgorithm(pset.get<fhicl::ParameterSet>("BezierTrackerAlgorithm"));
    
}

  void BezierTrackerModule::beginJob()
  {}


  void BezierTrackerModule::produce(art::Event& evt)
  {
 
    // Extract hits PtrVector from event

    art::Handle< std::vector<recob::Hit> > hith;
    evt.getByLabel(fHitModuleLabel, hith);

    std::vector<art::Ptr<recob::Hit> > HitVec;
    for(unsigned int i=0; i < hith->size(); ++i)
      {
	art::Ptr<recob::Hit> hit(hith,i);
	HitVec.push_back(hit);
      }
    

    // Declare products to store

    std::auto_ptr< std::vector<recob::Track > > btracks ( new std::vector<recob::Track>);
    std::auto_ptr< art::Assns<recob::Track, recob::Hit > > assn( new art::Assns<recob::Track, recob::Hit>);
   
    std::vector<trkf::BezierTrack * >         BTracks;
    std::vector<art::PtrVector<recob::Hit> >  HitsForAssns;
    
    
    if(fTrackMode==1)
      {
	// Look for track-like features in seed collections
	art::Handle< std::vector<recob::Seed> > seedh;
        evt.getByLabel(fSeedModuleLabel, seedh);

	std::vector<recob::Seed> TrackSeeds(*seedh);

	BTracks = fBTrackAlg->MakeBezierTracksFromSeeds(TrackSeeds);
	// Insert hit collecting code here
      }

   
    else if(fTrackMode==2)
      {
	// Find tracks in amorphous hit collection
	BTracks = fBTrackAlg->MakeBezierTracksFromHits(HitVec, HitsForAssns);
      }

    else if(fTrackMode==3)
      {
	// Find tracks from cluster combinations
	std::vector<std::vector<recob::SpacePoint> > PtsFromClusters = 
	  GetSpacePointsFromClusters(fClusterModuleLabel, evt);
	for(size_t i=0; i!=PtsFromClusters.size(); ++i)
	  {
	    std::vector<std::vector<recob::SpacePoint> > CataloguedSPs;
	    std::vector<recob::Seed> Seeds = fBTrackAlg->GetSeedFinderAlgorithm()->FindSeeds(PtsFromClusters.at(i),CataloguedSPs);
	    BTracks = fBTrackAlg->MakeBezierTracksFromSeeds(Seeds);
	  }
      }
    else if(fTrackMode==4)
      {
	// Made tracks from 

      }

    
    for(size_t i=0; i!=BTracks.size(); ++i)
      {
	recob::Track ToStore = BTracks.at(i)->GetBaseTrack();
	btracks->push_back(ToStore);
	util::CreateAssn(*this, evt, *(btracks.get()), HitsForAssns.at(i), *(assn.get()));
      }
  
    mf::LogInfo("BezierTrackerAlgorithm")<<"Storing in evt"<<std::endl;
    evt.put(btracks);
    evt.put(assn);
  }



  std::vector<std::vector<recob::SpacePoint> > BezierTrackerModule::GetSpacePointsFromClusters(std::string ClusterModuleLabel, art::Event& evt)
  {
    // Get Services.
    trkf::SpacePointAlg Sptalg = fBTrackAlg->GetSeedFinderAlgorithm()->GetSpacePointAlg();

    art::ServiceHandle<geo::Geometry> geom;

    std::vector<art::Ptr<recob::Cluster> > Clusters;

    art::Handle< std::vector<recob::Cluster> > clusterh;
    evt.getByLabel(ClusterModuleLabel, clusterh);

    if(clusterh.isValid()) {
      art::fill_ptr_vector(Clusters, clusterh);
    }
    art::FindManyP<recob::Hit> fm(clusterh, evt, ClusterModuleLabel);


    std::vector<std::vector<recob::SpacePoint> > SpacePointVectors;

    // Make a double or triple loop over clusters in distinct views
    // (depending on minimum number of views configured in SpacePointAlg).
    art::PtrVector<recob::Hit> hits;

    // Loop over first cluster.

    int nclus = Clusters.size();
    for(int iclus = 0; iclus < nclus; ++iclus) {
      art::Ptr<recob::Cluster> piclus = Clusters.at(iclus);
      geo::View_t iview = piclus->View();

      // Test first view.

      if((iview == geo::kU && Sptalg.enableU()) ||
         (iview == geo::kV && Sptalg.enableV()) ||
         (iview == geo::kW && Sptalg.enableW())) {

        // Store hits from first view into hit vector.

	std::vector< art::Ptr<recob::Hit> > ihits = fm.at(iclus);
        unsigned int nihits = ihits.size();
	hits.clear();
        hits.reserve(nihits);
	for(std::vector< art::Ptr<recob::Hit> >::const_iterator i = ihits.begin();
            i != ihits.end(); ++i)
          hits.push_back(*i);

        // Loop over second cluster.

        for(int jclus = 0; jclus < iclus; ++jclus) {
	  art::Ptr<recob::Cluster> pjclus = Clusters.at(jclus);
	  geo::View_t jview = pjclus->View();

	  // Test second view.

	  if(((jview == geo::kU && Sptalg.enableU()) ||
	      (jview == geo::kV && Sptalg.enableV()) ||
	      (jview == geo::kW && Sptalg.enableW()))
	     && jview != iview) {

	    // Store hits from second view into hit vector.
	    std::vector< art::Ptr<recob::Hit> > jhits = fm.at(jclus);
	    unsigned int njhits = jhits.size();
	    assert(hits.size() >= nihits);
	    //hits.resize(nihits);
	    while(hits.size() > nihits)
	      hits.pop_back();
	    assert(hits.size() == nihits);
	    hits.reserve(nihits + njhits);
	    for(std::vector< art::Ptr<recob::Hit> >::const_iterator j = jhits.begin();
		j != jhits.end(); ++j)
	      hits.push_back(*j);


	    // Loop over third cluster.

	    for(int kclus = 0; kclus < jclus; ++kclus) {
	      art::Ptr<recob::Cluster> pkclus = Clusters.at(kclus);
	      geo::View_t kview = pkclus->View();
	      // Test third view.

	      if(((kview == geo::kU && Sptalg.enableU()) ||
		  (kview == geo::kV && Sptalg.enableV()) ||
		  (kview == geo::kW && Sptalg.enableW()))
		 && kview != iview && kview != jview) {

		// Store hits from third view into hit vector.

		std::vector< art::Ptr<recob::Hit> > khits = fm.at(kclus);
		unsigned int nkhits = khits.size();
		assert(hits.size() >= nihits + njhits);
		//hits.resize(nihits + njhits);
		while(hits.size() > nihits + njhits)
		  hits.pop_back();
		assert(hits.size() == nihits + njhits);
		hits.reserve(nihits + njhits + nkhits);
		for(std::vector< art::Ptr<recob::Hit> >::const_iterator k = khits.begin();
		    k != khits.end(); ++k)
		  hits.push_back(*k);
		// Make three-view space points.

		std::vector<recob::SpacePoint> spts;
		Sptalg.makeSpacePoints(hits, spts);

		if(spts.size() > 0) {
		  SpacePointVectors.push_back(spts);
		}
	      }
	    }
	  }
	}
      }
    }


    return SpacePointVectors;
  }





  //----------------------------------------------------------------------
  void BezierTrackerModule::endJob()
  {

  }
}
