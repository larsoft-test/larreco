#ifndef SPACEPOINTFINDER_H
#define SPACEPOINTFINDER_H

//
// Name: SpacePointFinder.h
//
// Purpose: Header file for module SpacePointFinder.  This modules makes
//          generic prongs containing space points.
//
// Configuration parameters.
//
// ClusterModuleLabel;  // Cluster module label (e.g. "dbcluster").
// Filter;              // Filter space points?
// Merge;               // Merge space points?
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
    bool fFilter;                  // Filter space points?
    bool fMerge;                   // Merge space points?

    // Statistics.

    int fNumEvent;      // Number of events.
    int fNumProng2;     // Number of 2-view prongs.
    int fNumProng3;     // Number of 3-view prongs.
    int fNumSpt2;       // Number of 2-view space points.
    int fNumSpt3;       // Number of 3-view space points.
  };
}

#endif // SPACEPOINTFINDER_H
