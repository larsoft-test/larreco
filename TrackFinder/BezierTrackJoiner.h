#include "art/Persistency/Common/PtrVector.h"

#ifndef BEZIERTRACKEVBUILDER_H
#define BEZIERTRACKEVBUILDER_H

//
// Name: BezierTrackJoiner.h
//
// Purpose: Header file for module BezierTrackJoiner.  This modules makes
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
  class Track;
  class Hit;
}

class TVector3;

namespace trkf {

  class BezierTrack;

  class BezierTrackJoiner : public art::EDProducer
  {
  public:
 
    // Constructors, destructor

    explicit BezierTrackJoiner(fhicl::ParameterSet const& pset);
    virtual ~BezierTrackJoiner();

    std::vector<BezierTrack*> MakeTouchJoins(std::vector<BezierTrack*>, double);
    std::vector<TVector3> MakeVertexCandidates(std::vector<BezierTrack*>, std::vector<std::vector<int> >&);
    std::vector<BezierTrack*> MakeGlancingJoins(std::vector<BezierTrack*>, double);
    double GetTrackAngle(std::vector<trkf::BezierTrack*> BTracks, int i, int j, bool endi, bool endj);
    
    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void produce(art::Event& evt);
    void endJob();
    
    
    

  private:

    // Fcl Attributes.

    std::string fTrackModuleLabel;
    std::string fHitModuleLabel;

    
    double fJoinThreshold;
    double fVertexAngle;
    double fDistanceForAngle;
    double fHitDistance;
  };
}

#endif // SPACEPOINTFINDER_H
