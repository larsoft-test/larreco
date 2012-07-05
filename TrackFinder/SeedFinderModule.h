#ifndef SEEDFINDER_H
#define SEEDFINDER_H

//
// Name: SeedFinderModule.h
//
//
// Ben Jones, MIT, April 2012
//   bjpjones@mit.edu
//

#include "art/Framework/Core/EDProducer.h"
#include "TrackFinder/SpacePointAlg.h"
#include "TrackFinder/SeedFinderAlgorithm.h"

namespace recob
{
  class SpacePoint;
  class Cluster;
  class Seed;
  class Hit;
}

namespace trkf {

  class SeedFinderModule : public art::EDProducer
  {
  public:
 
    // Constructors, destructor

    explicit SeedFinderModule(fhicl::ParameterSet const& pset);
    virtual ~SeedFinderModule();

    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void produce(art::Event& evt);
    
    void endJob();
  
    
    std::vector<std::vector<recob::SpacePoint> > GetSpacePointsFromClusters(std::string fInputModuleLabel, art::Event& evt);

    art::PtrVector<recob::Hit>       GetHitsFromEvent(std::string HitModuleLabel, art::Event & evt);



  private:

   
    // Fcl Attributes.

    SeedFinderAlgorithm    fSeedAlg;                  // Algorithm object
    SpacePointAlg          fSptalg;
    std::string            fInputModuleLabel;         // Where to find hits, if we need them
    int                    fInputSource;              // 1: Use Clusters
                                                      // 2: Use Hits
    unsigned int           fMinPointsInCluster;       // How many SP's required in a cluster before
                                                      //  trying to seed
  };
  
}

#endif // SEEDFINDER_H
