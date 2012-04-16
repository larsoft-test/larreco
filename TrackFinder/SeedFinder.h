#ifndef SEEDFINDER_H
#define SEEDFINDER_H

//
// Name: SeedFinder.h
//
// Purpose: Header file for module SeedFinder.  This modules makes
//           seeds from clusters, via spacepoints.
//          Seeds must pass a point multiplicity requirement,
//           and a directionality cut.
//          There are 2 modes: 
//           1: Find one seed (if possible) in each cluster combo 
//           2: Find as many seeds as possible
//
// Configuration parameters.
//
// ClusterModuleLabel;  // Cluster module label (e.g. "dbcluster").
// Filter;              // Filter space points?
// Merge;               // Merge space points?
// fSeedMode;           // 1: Find 1 seed per track at high Z
//                      // 2: Find as many seeds as possible
// fSeedLength;         // Distance from highest Z point to search
//                      //  for hits
// fMinPointsInCluster; // How many SP's required in a cluster before
//                      //  trying to seed?
// fMinPointsInSeed;    // How many nearby points required to define
//                      //  a seed?
// fAngularDev;         // Angular standard deviation in radians
//                      //  for seed to meet directionality requirement

//
// Ben Jones, MIT, April 2012
//   bjpjones@mit.edu
//

#include "art/Framework/Core/EDProducer.h"

namespace recob
{
  class SpacePoint;
  class Cluster;
  class Seed;
  class Hit;
}

namespace trkf {

  class SeedFinder : public art::EDProducer
  {
  public:
 
    // Constructors, destructor

    explicit SeedFinder(fhicl::ParameterSet const& pset);
    virtual ~SeedFinder();

    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void produce(art::Event& evt);
    
    void endJob();
  
    


  private:

   


                   
    std::vector<std::vector<recob::SpacePoint> > GetSpacePointsFromClusters(std::string ClusterModuleLabel, art::Event & evt);
              // Get space points from clusters
              //    (I expect eventually we will get from event, which
              //     will be a separate method here) 

    std::vector<recob::SpacePoint> GetSpacePointsFromHitVector(art::PtrVector< recob::Hit> );
              // Generate space points from a pointer vector of hits 
              //  using the spacepoint service

    art::PtrVector<recob::Hit> GetHitsFromEvent(std::string HitModuleLabel, art::Event & evt);
              // Get hit vector from event 

    std::vector<std::vector<recob::SpacePoint> > GetSpacePointsFromHits(std::string HitModuleLabel, art::Event & evt);
              // Get hit vector from event 
 
    recob::Seed * FindSeedAtEnd(std::vector<recob::SpacePoint>);
              // Find one seed at a high Z end of a collection of spacepoints
              //    with no quality check


    
    std::vector<recob::Seed *> FindSeedExhaustively(std::vector<recob::SpacePoint>);
              // Generate seed at high Z end, check quality, discard if not adequate
              //    repeat until a strong seed found


 
    std::vector<recob::Seed *> FindAsManySeedsAsPossible(std::vector<recob::SpacePoint>);
              // Starting at high Z, find as many adequate seeds in the collection
              //     as possible

    
    // Fcl Attributes.

    std::string fClusterModuleLabel;           // Where to find clusters for this algorithm

    std::string fHitModuleLabel;           // Where to find hits, if we need them
 
    bool            fFilter, fMerge;           // Space point service settings
                                               //  See that service for info.
    int             fSeedMode;                 // 1: Find 1 seed per track at high Z
                                               // 2: Find as many seeds as possible 
    int             fSource;                   // 1: Use Clusters
                                               // 2: Use Hits
    double          fSeedLength;               // Distance from highest Z point to search
                                               //  for hits
    unsigned int    fMinPointsInCluster;       // How many SP's required in a cluster before
                                               //  trying to seed?
    unsigned int    fMinPointsInSeed;          // How many nearby points required to define
                                               //  a seed?
    float           fAngularDev;               // Angular standard deviation in radians
                                               //  for seed to meet directionality requirement


  };
  
}

#endif // SEEDFINDER_H
