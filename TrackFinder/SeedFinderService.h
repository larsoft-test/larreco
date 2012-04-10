////////////////////////////////////////////////////////////////////////
///
/// \file   SeedFinderService
///
/// \brief  Service for finding 3D track seeds 
///
/// \author B J P Jones
///
/// This service finds short straight line segments in 3D spacepoint
/// distributions which are subsequently used to seed tracks 
/// fit by a stepping procedure.
///
///


#ifndef SEEDFINDERSERVICE_H
#define SEEDFINDERSERVICE_H

#include <vector>
#include <string>
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "RecoBase/SpacePoint.h"
#include "Simulation/SimChannel.h"

#include "RecoBase/Seed.h"

namespace trkf {



  class SeedFinderService {
  public:

    // Constructor.
    SeedFinderService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);

    // Destructor.
    ~SeedFinderService();

    // Update configuration parameters.
    void reconfigure(const fhicl::ParameterSet& pset);

    recob::Seed * FindSeedAtEnd(std::vector<recob::SpacePoint>);
    std::vector<recob::Seed *> FindSeedExhaustively(std::vector<recob::SpacePoint>);
    std::vector<recob::Seed *> FindAsManySeedsAsPossible(std::vector<recob::SpacePoint>);
    void ProduceSpacePointPlots(std::vector<std::vector<recob::SpacePoint> > Points);
    void SetEventID(int EventID);

  private:

    double          fSeedLength;
    unsigned int    fMinPointsInCluster;
    unsigned int    fMinPointsInSeed;
    float           fMinDirectionStrength;
    int             fEventID;

  };



}

#endif
