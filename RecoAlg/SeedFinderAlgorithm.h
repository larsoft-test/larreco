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
#include "RecoAlg/SpacePointAlg.h"
#include "TVector3.h"
#include "Geometry/Geometry.h"
#include "TTree.h"

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

 
    std::vector<recob::Seed >  FindSeeds(std::vector<recob::SpacePoint> const&, std::vector<std::vector<recob::SpacePoint> >&);

    void                        RefitSeed(recob::Seed& TheSeed, std::vector<recob::SpacePoint> const& SpacePoints);

    std::vector<recob::SpacePoint> GetSpacePointsFromHitVector(art::PtrVector<recob::Hit> const& );

    std::vector<std::vector<recob::Seed> > GetSeedsFromClusterHits(std::map<geo::View_t, std::vector<art::PtrVector<recob::Hit> > > const &);


    std::vector<double>         GetHitRMS(recob::Seed const& TheSeed, std::vector<recob::SpacePoint> const&);

    size_t                      CountHits(std::vector<recob::SpacePoint> const& Points);
    
    SpacePointAlg *             GetSpacePointAlg() const { return fSptalg; }
    
    size_t                      CountHits(std::vector<recob::SpacePoint>  const& SpacePoints, TVector3 CenterPoint, double d);

    void                        GetCenterAndDirection(std::vector<recob::SpacePoint> const& Points, std::vector<int> const& PointsInRange, TVector3& Center, TVector3& Direction, double& Strength, bool Mode);
    
  private:

    recob::Seed                 FindSeedAtEnd(std::vector<recob::SpacePoint> const&, std::map<int, int>&, std::vector<int>&);

    bool                        ExtendSeed(recob::Seed& TheSeed, std::vector<recob::SpacePoint> const& AllSpacePoints, 
					   std::map<int, int>& PointStatus, std::vector<int>& PointsUsed);

    std::vector<recob::SpacePoint> ExtractSpacePoints(std::vector<recob::SpacePoint> const& AllPoints, std::vector<int> IDsToExtract);

    std::vector<int>            DetermineNearbySPs(recob::Seed const& TheSeed, std::vector<recob::SpacePoint> const& AllSpacePoints, 
						   std::map<int, int> PointStatus, double ExtendResolution);


                       
    // Fcl Attributes.

    SpacePointAlg       *  fSptalg;            
 
    double                fInitSeedLength;    
                                        
    int                   fMinPointsInSeed;   
                                        
    float                 fPCAThreshold;
                                        
    int                   fRefits;            
    
    std::vector<double>   fMaxViewRMS;

    float                 fExtendThresh;
    float                 fExtendStep;
    float                 fExtendResolution;

    TTree *   ftMonitoringTree;
    Float_t   ftThetaXZ;
    Float_t   ftThetaYZ;
    Float_t   ftTheta;
    Float_t   ftEigenvalue;
    Int_t     ftNSpts;
    Int_t     ftNUHits;
    Int_t     ftNVHits;
    Int_t     ftNWHits;
    Float_t   ftURMS;
    Float_t   ftVRMS;
    Float_t   ftWRMS;
    Float_t   ftURMSb;
    Float_t   ftVRMSb;
    Float_t   ftWRMSb;
    bool      ftKeep;
  };
  
}

#endif // SEEDFINDER_H
