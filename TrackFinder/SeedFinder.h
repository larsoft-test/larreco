#ifndef SEEDFINDER_H
#define SEEDFINDER_H

//
// Name: SeedFinder.h
//
// Purpose: Header file for module SeedFinder.  This modules makes
//          seeds using the SeedFinderService
//
// Configuration parameters.
//
// ClusterModuleLabel;  // Cluster module label (e.g. "dbcluster").
// Filter;              // Filter space points?
// Merge;               // Merge space points?
//
// Ben Jones, MIT
//

#include "art/Framework/Core/EDProducer.h"

namespace recob
{
  class SpacePoint;
  class Cluster;
}

namespace trkf {

  class SeedFinder : public art::EDProducer
  {
  public:
 
    // Constructors, destructor

    explicit SeedFinder(fhicl::ParameterSet const& pset);
    virtual ~SeedFinder();

    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void produce(art::Event& evt);
    void endJob();
    
    std::vector<std::vector<recob::SpacePoint> > GetSpacePoints(art::Handle<std::vector<recob::Cluster> >);
    

  private:

    // Fcl Attributes.

    std::string fClusterModuleLabel;
    bool fFilter, fMerge;
    
    
  };
}

#endif // SPACEPOINTFINDER_H
