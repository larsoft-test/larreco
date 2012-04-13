#include "art/Persistency/Common/PtrVector.h"

#ifndef BEZIERTRACKER_H
#define BEZIERTRACKER_H

//
// Name: BezierTracker.h
//
// Purpose: Header file for module BezierTracker.  This modules makes
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

namespace recob
{
  class Seed;
  class BezierTrackBase;
}

namespace trkf {

  class BezierTracker : public art::EDProducer
  {
  public:
 
    // Constructors, destructor

    explicit BezierTracker(fhicl::ParameterSet const& pset);
    virtual ~BezierTracker();

    std::vector<std::vector<art::Ptr<recob::Seed> > > OrganizeSeedsIntoTracks(std::vector<art::Ptr<recob::Seed> > SeedVector);
    // std::vector<art::PtrVector<recob::Seed> OrganizeSeedsIntoTracks(art::PtrVector<recob::Seed>  SeedVector);
    recob::BezierTrackBase ProduceBaseTrack(std::vector<art::Ptr<recob::Seed> > Seeds);

    
    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void produce(art::Event& evt);
    void endJob();
    
    
    

  private:

    // Fcl Attributes.

    std::string fSeedModuleLabel;
    std::string fHitModuleLabel;
    double fMaxKinkAngle;
    double fMaxTrackMissAngle;
    double fMaxJumpDistance;
    int fTopTrackID;
  };
}

#endif // SPACEPOINTFINDER_H
