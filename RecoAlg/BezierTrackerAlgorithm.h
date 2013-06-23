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
#include "RecoAlg/SeedFinderAlgorithm.h"

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

    void MakeBezierTracksFromSeeds(std::vector<trkf::BezierTrack>& ReturnVector,
				   std::vector<recob::Seed> const& TrackSeeds  );

    void MakeBezierTracksFromHits(std::vector<trkf::BezierTrack>& ReturnVector, 
				  std::vector<art::Ptr<recob::Hit> > HitVec, 
				  std::vector<art::PtrVector<recob::Hit> >& HitsForAssns );

    std::vector<std::vector<recob::Seed> > OrganizeSeedsIntoTracks(std::vector<recob::Seed > const&  SeedVector);


    std::vector<int> DetermineNearbyHits(art::PtrVector<recob::Hit> const& Hits, 
					 BezierTrack const& BTrack, 
					 std::vector<double>& SValues);
    
    void MakeDirectJoins(std::vector<trkf::BezierTrack>& BTracks, std::vector<art::PtrVector<recob::Hit> > & HitVecs);
    
    void AddPtrVectors(art::PtrVector<recob::Hit>& Receiever, art::PtrVector<recob::Hit> const & ToAdd);
      
      
    
    trkf::SeedFinderAlgorithm * GetSeedFinderAlgorithm() { return fTheSeedFinder;}

    
    
    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    
    

  private:

    // Fcl Attributes.

    double fMaxKinkDThetaDx;
    double fMaxJumpLengths;
    double fHitDistance;
    double fTrackJoinAngle;
 
    SeedFinderAlgorithm * fTheSeedFinder;

  };
}

#endif
