////////////////////////////////////////////////////////////////////////
///
/// \file   KalmanFilterAlg.h
///
/// \brief  Kalman Filter.
///
/// \author H. Greenlee 
///
/// Configuration parameters:
///
/// Trace        - Trace flag.
/// MaxPErr      - Maximum pointing error for free propagation.
/// GoodPErr     - Pointing error threshold for switching to free propagation.
/// MaxIncChisq  - Maximum incremental chisquare to accept a hit.
/// MinLHits     - Minimum number of hits to turn off linearized propagation.
/// MaxLDist     - Maximum distance for linearized propagation.
/// MaxPredDist  - Maximum prediciton distance to accept a hit.
/// MaxPropDist  - Maximum propagation distance to candidate surface.
/// MinSortDist  - Sort low distance threshold.
/// MaxSortDist  - Sort high distance threshold.
/// MaxSamePlane - Maximum consecutive hits in same plane.
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

    /// Make a new track.
    bool buildTrack(const KTrack& trk,                     // Starting track.
		    KGTrack& trg,                          // Result global track.
		    const Propagator* prop,                // Propagator.
		    const Propagator::PropDirection dir,   // Direction.
		    KHitContainer& hits) const;            // Candidate measurements.

    /// Smooth track.
    bool smoothTrack(KGTrack& trg,                         // Global track to be smoothed.
		     KGTrack* trg1,                        // Result of unidirectional fit.
		     const Propagator* prop) const;        // Propagator.

    /// Add hits to existing track.
    bool extendTrack(KGTrack& trg,                         // Global track.
		     const Propagator* prop,               // Propagator.
		     KHitContainer& hits) const;           // Candidate measurements.

    /// Iteratively smooth a track.
    bool smoothTrackIter(int niter,                        // Number of iterations.
			 KGTrack& trg,                     // Global track.
			 const Propagator* prop) const;    // Propagator.

  private:

    // Fcl parameters.

    bool fTrace;             ///< Trace flag.
    double fMaxPErr;         ///< Maximum pointing error for free propagation.
    double fGoodPErr;        ///< Pointing error threshold for switching to free propagation.
    double fMaxIncChisq;     ///< Maximum incremental chisquare to accept a hit.
    unsigned int fMinLHits;  ///< Minimum number of hits to turn off linearized propagation.
    double fMaxLDist;        ///< Maximum distance for linearized propagation.
    double fMaxPredDist;     ///< Maximum prediciton distance to accept a hit.
    double fMaxPropDist;     ///< Maximum propagation distance to candidate surface.
    double fMinSortDist;     ///< Sort low distance threshold.
    double fMaxSortDist;     ///< Sort high distance threshold.
    int fMaxSamePlane;       ///< Maximum consecutive hits in same plane.

    // Other attributes.

    int fPlane;          ///< Preferred view plane.
  };
}

#endif
