#ifndef SPACEPOINTFINDER_H
#define SPACEPOINTFINDER_H

//
// Name: SpacePointFinder.h
//
// Purpose: Header file for module SpacePointFinder.
//
// Configuration parameters.
//
// ClusterModuleLabel;  // Cluster module label (e.g. "dbcluster").
//
// Created: 15-Dec-2011  H. Greenlee
//

#include "art/Framework/Core/EDProducer.h"
#include "TrackFinder/SpacePointAlg.h"

namespace trkf {

  class SpacePointFinder : public art::EDProducer
  {
  public:
 
    // Constructors, destructor

    explicit SpacePointFinder(fhicl::ParameterSet const& pset);
    virtual ~SpacePointFinder();

    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void produce(art::Event& evt);
    void endJob();

  private:

    // Fcl Attributes.

    SpacePointAlg fSptalg;         // Algorithm object.
    std::string fClusterModuleLabel;
    unsigned int fMinHits;         // Minimum number of hits per cluster.
    bool fClusterAssns;            // Make Cluster-SpacePoint associations.

    // Statistics.

    int fNumEvent;      // Number of events.
    int fNumSpt2;       // Number of 2-view space points.
    int fNumSpt3;       // Number of 3-view space points.
  };
}

#endif // SPACEPOINTFINDER_H
