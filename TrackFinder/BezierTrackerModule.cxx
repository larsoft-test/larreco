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
	art::Handle< std::vector<recob::Seed> > seedh;
        evt.getByLabel(fSeedModuleLabel, seedh);
	
	std::vector<art::Ptr<recob::Seed> > TrackSeeds;
        if(seedh->size()>0)
          for(unsigned int iseed=0; iseed!=seedh->size(); iseed++)
            {
	      art::Ptr<recob::Seed> theseed(seedh, iseed);
              TrackSeeds.push_back(theseed);
            }	
	BTracks = fBTrackAlg->MakeBezierTracksFromSeeds(TrackSeeds);
	// Insert hit collecting code here
      }

   
    else if(fTrackMode==2)
      {
	BTracks = fBTrackAlg->MakeBezierTracksFromHits(HitVec, HitsForAssns);
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



  //----------------------------------------------------------------------
  void BezierTrackerModule::endJob()
  {

  }
}
