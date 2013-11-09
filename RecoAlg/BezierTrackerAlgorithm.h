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


    std::vector<int> DetermineNearbyHits(art::PtrVector<recob::Hit> const& Hits, 
					 BezierTrack const& BTrack, 
					 std::vector<double>& SValues);
    
    void FilterOverlapTracks(std::vector<trkf::BezierTrack>& BTracks, std::vector<art::PtrVector<recob::Hit> > & HitVecs);

    void MakeDirectJoins(std::vector<trkf::BezierTrack>& BTracks, std::vector<art::PtrVector<recob::Hit> > & HitVecs);
    
    void AddPtrVectors(art::PtrVector<recob::Hit>& Receiever, art::PtrVector<recob::Hit> const & ToAdd);
      
      
    
    trkf::SeedFinderAlgorithm * GetSeedFinderAlgorithm() { return fTheSeedFinder;}

    void  MakeVertexJoins(std::vector<trkf::BezierTrack>& BTracks, std::vector<recob::Vertex>& Vertices, std::vector<std::vector<int> > Mapping);
    
    void  FilterAndJoin(std::vector<std::vector<std::vector<recob::Seed> > > Seeds, std::vector<std::vector<std::vector<std::vector<int> > > > HitsPerSeed, size_t UEntries, size_t VEntries, size_t WEntries);


    std::vector<trkf::BezierTrack> MakeTracksNew(std::map<geo::View_t, std::vector<art::PtrVector<recob::Hit> > >& SortedHits, std::vector<art::PtrVector<recob::Hit> >& HitAssocs);
     
    void GetTracksForCombo(std::vector<recob::Seed>& Seeds, art::PtrVector<recob::Hit>& UHits, art::PtrVector<recob::Hit>& VHits, art::PtrVector<recob::Hit>& WHits);

    std::vector<std::vector< recob::Seed > > OrganizeSeedsIntoTracksNew(std::vector<recob::Seed >& AllSeeds, std::vector<art::PtrVector<recob::Hit> * >& AllHits, std::vector<art::PtrVector<recob::Hit> >& WhichHitsPerSeed, std::vector<std::map<uint32_t, std::vector<int> >* >& OrgHits, std::vector<std::vector<std::vector<int> > >& WhichHitsPerTrack);

    void GetSeedDirProjected(recob::Seed const& TheSeed, std::vector<double>& WireCoord, std::vector<double>& TimeCoord);

    std::vector<double> GetOccupancy(recob::Seed& Seed1, recob::Seed& Seed2, double dThresh,  std::vector<art::PtrVector<recob::Hit>*>& AllHits,  std::vector<std::map<uint32_t, std::vector<int> >* >& OrgHits,  std::vector<uint32_t>& LowChan, std::vector<uint32_t>& HighChan, std::vector<std::vector<int> >& HitStatus, std::vector<std::vector<int> >& TheseHits);

    void SortTracksByLength(std::vector<trkf::BezierTrack>& BTracks, std::vector<art::PtrVector<recob::Hit> > & HitVecs);
    
    void CalculateGeometricalElements();



    
    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    
    

  private:

    // Fcl Attributes.

    

    
    double fOverlapCut;
    double fDirectJoinDistance;
    double fTrackJoinAngle;
    std::vector<double> fOccupancyThresh;
    double fTrackResolution;

    double fVertexImpactThreshold;
    double fVertexExtrapDistance;

    std::vector<double>   fPitches;
    std::vector<TVector3> fPitchDir;
    std::vector<TVector3> fWireDir;
    std::vector<double>   fWireZeroOffset;
    TVector3              fXDir, fYDir, fZDir;


    void GetImpact(TVector3 t1pt, TVector3 t1dir, TVector3 t2pt, TVector3 t2dir, double& ImpactParam, double& Dist1, double& Dist2);

 


    SeedFinderAlgorithm * fTheSeedFinder;

  };
}

#endif
