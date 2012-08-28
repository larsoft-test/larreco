#include "art/Persistency/Common/PtrVector.h"

#ifndef BEZIERTRACKERALG_H
#define BEZIERTRACKERALG_H

//
// Name: BezierTrackerAlgorithm.h
//
// Purpose: Header file for module BezierTrackerAlgorithm.  This modules makes
//          bezier tracks out of seed collections
//
// Configuration parameters.
//
// SeedModuleLabel;     // Cluster module label (e.g. "dbcluster").
// HitModuleLabel;      // Hit module label (e.g. "FFTHitFinder")
//
// Ben Jones, MIT
//

#include "art/Framework/Core/EDProducer.h"
#include "TrackFinder/SeedFinderAlgorithm.h"

namespace recob
{
  class Seed;
  class Track;
  class Hit;
}


namespace trkf {

  class BezierTrack;
 
  class BezierTrackerAlgorithm
  {
  public:
 
    // Constructors, destructor

    explicit BezierTrackerAlgorithm(fhicl::ParameterSet const& pset);
    virtual ~BezierTrackerAlgorithm();

    std::vector<trkf::BezierTrack* > MakeBezierTracksFromSeeds(std::vector<recob::Seed> const& TrackSeeds  );

    std::vector<trkf::BezierTrack* > MakeBezierTracksFromHits(std::vector<art::Ptr<recob::Hit> > HitVec, std::vector<art::PtrVector<recob::Hit> >& HitsForAssns );



    std::vector<std::vector<recob::Seed> > OrganizeSeedsIntoTracks(std::vector<recob::Seed >  SeedVector);


    BezierTrack* ProduceTrackFromSeeds(std::vector<recob::Seed> const& Seeds);

    std::vector<int> DetermineNearbyHits(art::PtrVector<recob::Hit> Hits, BezierTrack * BTrack, double HitCollectionDistance, std::vector<double>& SValues);
    
    trkf::SeedFinderAlgorithm * GetSeedFinderAlgorithm() { return fTheSeedFinder;}
    
    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    
    

  private:

    // Fcl Attributes.

    double fMaxKinkAngle;
    double fMaxTrackMissAngle;
    double fMaxJumpDistance;
    double fHitDistance;

    int fTopTrackID;

    SeedFinderAlgorithm * fTheSeedFinder;

  };
}

#endif
