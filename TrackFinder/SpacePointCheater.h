#ifndef SPACEPOINTCHEATER_H
#define SPACEPOINTCHEATER_H

//
// Name: SpacePointCheater.h
//
// Purpose: Header file for module SpacePointCheater.  This modules makes
//          generic prongs containing space points using mc truth information.
//
// Configuration parameters.
//
// ClusterModuleLabel;  // Cluster module label (e.g. "dbcluster").
// G4ModuleLabel;       // For SimChannel (e.g. "largeant").
// Filter;              // Filter space points?

//
// Created: 15-Dec-2011  H. Greenlee
//

#include "art/Framework/Core/EDProducer.h"

namespace trkf {

  class SpacePointCheater : public art::EDProducer
  {
  public:
 
    // Constructors, destructor

    explicit SpacePointCheater(fhicl::ParameterSet const& pset);
    virtual ~SpacePointCheater();

    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void produce(art::Event& evt);
    void endJob();

  private:

    // Fcl Attributes.

    std::string fClusterModuleLabel;
    std::string fG4ModuleLabel;    // For SimChannel.
    bool fFilter;                  // Filter space points?

    // Statistics.

    int fNumEvent;      // Number of events.
    int fNumProng2;     // Number of 2-view prongs.
    int fNumProng3;     // Number of 3-view prongs.
  };
}

#endif // SPACEPOINTCHEATER_H
