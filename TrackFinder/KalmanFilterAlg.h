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

#include "TrackFinder/KGTrack.h"
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

    /// Reconfigure method.
    void reconfigure(const fhicl::ParameterSet& pset);

    // Accessors.

    bool getTrace() const {return fTrace;}      ///< Trace config parameters.
    int getPlane() const {return fPlane;}       ///< Preferred view plane.

    // Modifiers.

    void setTrace(bool trace) {fTrace = trace;} ///< Set trace config parameter.
    void setPlane(int plane) {fPlane = plane;}  ///< Set preferred view plane.

    // Methods.

    /// Add hits to track.
    bool buildTrack(const KTrack& trk,                     // Starting track.
		    KGTrack& trg,                          // Result global track.
		    const Propagator* prop,                // Propagator.
		    const Propagator::PropDirection dir,   // Direction.
		    KHitContainer& hits,                   // Candidate measurements.
		    double maxChisq);                      // Incremental chisquare cut.

    /// Smooth track.
    bool smoothTrack(KGTrack& trg,                         // Global track to be smoothed.
		     const Propagator* prop);              // Propagator.

  private:

    // Fcl parameters.

    bool fTrace;       ///< Trace flag.

    // Other attributes.

    int fPlane;        ///< Preferred view plane.
    TrackError fErr;   ///< Starting error matrix ("infinite" error).
  };
}

#endif
