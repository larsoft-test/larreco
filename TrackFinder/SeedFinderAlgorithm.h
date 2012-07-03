#ifndef SEEDFINDERALG_H
#define SEEDFINDERALG_H

//
// Name: SeedFinderAlgorithm.h
//
//
// Ben Jones, MIT, April 2012
//   bjpjones@mit.edu
//

#include "art/Framework/Core/EDProducer.h"
#include "TrackFinder/SpacePointAlg.h"

namespace recob
{
  class SpacePoint;
  class Seed;
  class Hit;
}

namespace trkf {

  class SeedFinderAlgorithm
  {
  public:
 
    // Constructors, destructor

    SeedFinderAlgorithm(const fhicl::ParameterSet& pset);
   ~SeedFinderAlgorithm();


    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);

 
    // Seedfinding methods




    recob::Seed *               FindSeedAtEnd(std::vector<recob::SpacePoint>,std::map<int, int>, std::vector<int>&);

    std::vector<recob::Seed *>  FindSeeds(std::vector<recob::SpacePoint>, std::vector<std::vector<recob::SpacePoint> >&);

    void                        RefitSeed(recob::Seed * TheSeed, std::vector<recob::SpacePoint> SpacePoints);
    
    bool                        ExtendSeed(recob::Seed* TheSeed, std::vector<recob::SpacePoint> AllSpacePoints, 
					   std::vector<int> PointsUsed, std::map<int, int> PointStatus);

    std::vector<recob::SpacePoint> GetSpacePointsFromHitVector(art::PtrVector< recob::Hit> );

    SpacePointAlg const& GetSpacePointAlg() const { return fSptalg; }

    
  private:

    double CountHits(std::vector<recob::SpacePoint> Points);

    std::vector<double> GetHitRMS(recob::Seed* TheSeed, std::vector<recob::SpacePoint>);

                       
    // Fcl Attributes.

    SpacePointAlg         fSptalg;            
 
    double                fInitSeedLength;    
                                        
    unsigned int          fMinPointsInSeed;   
                                        
    float                 fAngularDev;        
                                        
    int                   fRefits;            
    
    std::vector<double>   fMaxViewRMS;

    float                 fExtendThresh;

  };
  
}

#endif // SEEDFINDER_H
