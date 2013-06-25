#include "art/Persistency/Common/PtrVector.h"
#include "TVector3.h"

#ifndef BEZIERTRACKERALG_H
#define BEZIERTRACKERALG_H

//
// Name: BezierTrackerAlgorithm.h
//
// Purpose: Header file for module BezierTrackerAlgorithm.  This 
//  algorithm contains tools for producing and combining bezier
//  tracks made up of seed segments.
////
// Ben Jones, MIT
//

#include "art/Framework/Core/EDProducer.h"
#include "RecoAlg/SeedFinderAlgorithm.h"

namespace recob
{
  class Seed;
  class Track;
  class Hit;
  class Vertex;
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

    void  MakeVertexJoins(std::vector<trkf::BezierTrack>& BTracks, std::vector<recob::Vertex>& Vertices, std::vector<std::vector<int> > Mapping);
    
    
    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    
    

  private:

    // Fcl Attributes.

    

    
    double fMaxJumpLengths;
    double fHitDistance;
    double fDirectJoinDistance;
    double fTrackJoinAngle;

    double fVertexImpactThreshold;
    double fVertexExtrapDistance;

    void GetImpact(TVector3 t1pt, TVector3 t1dir, TVector3 t2pt, TVector3 t2dir, double& ImpactParam, double& Dist1, double& Dist2);

 


    SeedFinderAlgorithm * fTheSeedFinder;

  };
}

#endif
