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
/// FCL parameters:
///


#ifndef SPACEPOINTSERVICE_H
#define SPACEPOINTSERVICE_H

#include <vector>
#include <string>
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "RecoBase/SpacePoint.h"
#include "Simulation/SimChannel.h"


namespace trkf {

  class Seed;


  class SeedFinderService {
  public:

    // Constructor.
    SeedFinderService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);

    // Destructor.
    ~SeedFinderService();

    // Update configuration parameters.
    void reconfigure(const fhicl::ParameterSet& pset);

    Seed FindSeedAtEnd(std::vector<recob::SpacePoint>);
    Seed FindSeedExhaustively(std::vector<recob::SpacePoint>);
    void ProduceSpacePointPlots(std::vector<std::vector<recob::SpacePoint> > Points);
    void SetEventID(int EventID);

  private:

    double fSeedLength;
    unsigned int    fMinPointsInCluster;
    unsigned int    fMinPointsInSeed;
    float  fMinDirectionStrength;
    int fEventID;

  };




  class Seed
  {
  public:
    Seed();
    Seed(TVector3 Point, TVector3 Direction);
    ~Seed();

    TVector3 GetDirection();
    TVector3 GetPoint();

    void SetDirection(TVector3 Direction);
    void SetPoint(TVector3 Point);



    bool IsValid();
    void SetValidity(bool Validity);


  private:
    TVector3 fSeedPoint;
    TVector3 fSeedDirection;
    bool     fIsValid;
  };





}

#endif
