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

    //--------------------------------------
    // Constructors, destructor, reconfigure
    //--------------------------------------

    SeedFinderAlgorithm(const fhicl::ParameterSet& pset);
   ~SeedFinderAlgorithm();

    void reconfigure(fhicl::ParameterSet const& pset);

 

    //----------------------
    // Seedfinding methods
    //----------------------


    std::vector<std::vector<recob::Seed> > GetSeedsFromSortedHits( std::map<geo::View_t, std::vector<art::PtrVector<recob::Hit> > >  const& SortedHits, 
								   std::vector<std::vector<art::PtrVector<recob::Hit> > >& HitsPerSeed);
                                    // Return a vector of vectors of seeds, one vector for each supplied cluster 
                                    //   combination which has sufficient overlap. The second argument returns
                                    //   the hits sorted by combo and by seed
   
    

    std::vector<recob::Seed>    GetSeedsFromUnSortedHits(art::PtrVector<recob::Hit> const &, std::vector<art::PtrVector<recob::Hit> >&);
                                    // Return a vector of seeds formed from an unstructured collection of hits    



    //----------------------
    // Alg passing 
    //----------------------


    SpacePointAlg *             GetSpacePointAlg() const { return fSptalg; }
                                    // Return the SpacePointAlg, as configured for the Seed Finding 
   
    


    
  private:

    //----------------------
    // Internal methods
    //----------------------


    std::vector<recob::Seed>    FindSeeds( art::PtrVector<recob::Hit> const& HitsFlat, std::vector<art::PtrVector<recob::Hit> >& CataloguedHits);
                                    // Find a collection of seeds, based on the supplied set of hits.
                                    //  The second argument returns the hits catalogued by which
                                    //  seed they fell into (if any) 
   


    recob::Seed                 FindSeedAtEnd(std::vector<recob::SpacePoint> const&, std::vector<char>&, std::vector<int>&);
                                    // Find one seed at high Z from the spacepoint collection given. Latter arguments are 
                                    //  for internal book keeping.


    void                        RefitSeed(recob::Seed& TheSeed, std::vector<recob::SpacePoint> const& SpacePoints);
                                   // Having found a 3D seed, refit it onto its constituent hits to iteratively minimize 
                                   //  the RMS in each view


    std::vector<double>         GetHitRMS(recob::Seed const& TheSeed, std::vector<recob::SpacePoint> const&);
                                   // Get a vector of RMS values indiciating how closely the seed runs to the hits in each view


    size_t                      CountHits(std::vector<recob::SpacePoint> const& Points);
                                   // Counting the number of hits in each view which are associated with a set of SPs


    void                        GetCenterAndDirection(std::vector<recob::SpacePoint> const& Points, std::vector<int> const& PointsInRange, TVector3& Center, TVector3& Direction, double& Strength, int Mode);
                                   // Given a set of spacepoints, find the center and direction to form a seed. 
                                   //  Mode specifies whether to operate on spacepoints (old) or directly onto hits (new)

    
    void                        ConsolidateSeed(recob::Seed& TheSeed, art::PtrVector<recob::Hit> const&, std::vector<char>& HitStatus,
						std::map<geo::View_t, std::map<uint32_t, std::vector<int> > >& OrgHits, bool Extend);

    void                        GetHitDistAndProj( recob::Seed const& ASeed,  art::Ptr<recob::Hit> const& AHit, double& disp, double& s);
    
    std::vector<recob::SpacePoint> ExtractSpacePoints(std::vector<recob::SpacePoint> const& AllPoints, std::vector<int> IDsToExtract);
                                   // Given a spacepoint vector and ID list, return the vector of spacepoints for those IDs
 

    std::vector<int>            DetermineNearbySPs(recob::Seed const& TheSeed, std::vector<recob::SpacePoint> const& AllSpacePoints, 
						   std::map<int, int> PointStatus, double ExtendResolution);

                                   // Find the spacepoints which are within some radius of seed from a provided sed.



                       
    // Fcl Attributes.

    SpacePointAlg       *  fSptalg;            
 
    double                fInitSeedLength;    
                                        
    int                   fMinPointsInSeed;   
                                        
    float                 fPCAThreshold;
                                        
    int                   fRefits;            
    
    std::vector<double>   fMaxViewRMS;

    float                 fHitResolution;

    float                 fOccupancyCut;

    bool                  fExtendSeeds;

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
