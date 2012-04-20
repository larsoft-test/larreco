////////////////////////////////////////////////////////////////////////
///
/// \file   KalmanFilterService.h
///
/// \brief  Kalman Filter.
///
/// \author H. Greenlee 
///
////////////////////////////////////////////////////////////////////////

#ifndef KALMANFILTERSERVICE_H
#define KALMANFILTERSERVICE_H

#include "TrackFinder/KHitsTrack.h"
#include "TrackFinder/Propagator.h"
#include "TrackFinder/KHitContainer.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"

namespace trkf {

  class KalmanFilterService {
  public:

    /// Constructor.
    KalmanFilterService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);

    /// Destructor.
    ~KalmanFilterService();

    // Methods.

    /// Add hits to track.
    void buildTrack(KHitsTrack& trh,                       // Starting track.
		    const Propagator* prop,                // Propagator.
		    const Propagator::PropDirection dir,   // Direction.
		    KHitContainer& hits,                   // Candidate measurements.
		    double maxChisq);                      // Incremental chisquare cut.
  };
}

#endif
