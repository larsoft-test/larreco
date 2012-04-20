////////////////////////////////////////////////////////////////////////
///
/// \file   KalmanFilterAlg.h
///
/// \brief  Kalman Filter.
///
/// \author H. Greenlee 
///
////////////////////////////////////////////////////////////////////////

#ifndef KALMANFILTERALG_H
#define KALMANFILTERALG_H

#include "TrackFinder/KHitsTrack.h"
#include "TrackFinder/Propagator.h"
#include "TrackFinder/KHitContainer.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/PtrVector.h"

namespace trkf {

  class KalmanFilterAlg {
  public:

    /// Constructor.
    KalmanFilterAlg(const fhicl::ParameterSet& pset);

    /// Destructor.
    ~KalmanFilterAlg();

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
