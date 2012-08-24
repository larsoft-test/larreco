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

    
    std::vector<recob::Seed* >  FindSeeds(std::vector<std::vector<recob::SpacePoint> > const&, std::vector<std::vector<recob::SpacePoint> > & );
    

    std::vector<recob::Seed *>  FindSeeds(std::vector<recob::SpacePoint> const&, std::vector<std::vector<recob::SpacePoint> >&);

    void                        RefitSeed(recob::Seed * TheSeed, std::vector<recob::SpacePoint> SpacePoints);

    std::vector<recob::SpacePoint> GetSpacePointsFromHitVector(art::PtrVector<recob::Hit> );

    std::vector<double>         GetHitRMS(recob::Seed* TheSeed, std::vector<recob::SpacePoint>);

    double                      CountHits(std::vector<recob::SpacePoint> Points);
    
    SpacePointAlg const&        GetSpacePointAlg() const { return fSptalg; }



  private:

    recob::Seed *               FindSeedAtEnd(std::vector<recob::SpacePoint> const&, std::map<int, int>&, std::vector<int>&);

    bool                        ExtendSeed(recob::Seed* TheSeed, std::vector<recob::SpacePoint> const& AllSpacePoints, 
					   std::map<int, int>& PointStatus, std::vector<int>& PointsUsed);

    std::vector<recob::SpacePoint> ExtractSpacePoints(std::vector<recob::SpacePoint> AllPoints, std::vector<int> IDsToExtract);

    std::vector<int>            DetermineNearbySPs(recob::Seed* TheSeed, std::vector<recob::SpacePoint> AllSpacePoints, 
						   std::map<int, int> PointStatus, double ExtendResolution);


                       
    // Fcl Attributes.

    SpacePointAlg         fSptalg;            
 
    double                fInitSeedLength;    
                                        
    unsigned int          fMinPointsInSeed;   
                                        
    float                 fAngularDev;        
                                        
    int                   fRefits;            
    
    std::vector<double>   fMaxViewRMS;

    float                 fExtendThresh;
    float                 fExtendStep;
    float                 fExtendResolution;
  };
  
}

#endif // SEEDFINDER_H
