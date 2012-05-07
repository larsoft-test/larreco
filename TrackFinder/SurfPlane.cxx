///////////////////////////////////////////////////////////////////////
///
/// \file   SurfPlane.cxx
///
/// \brief  Base class for Kalman filter planar surfaces.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "TrackFinder/SurfPlane.h"

namespace trkf {

  /// Default constructor.
  SurfPlane::SurfPlane()
  {}

  /// Destructor.
  SurfPlane::~SurfPlane()
  {}

  /// Get pointing error of track.
  ///
  /// Arguments:
  ///
  /// vec - Track parameters.
  /// err - Track error matrix.
  ///
  /// Returns: Pointing error.
  ///
  /// This method calculates the track pointing error based on the
  /// slope track paramers and errors (parameters 2 and 3).
  ///
  double SurfPlane::PointingError(const TrackVector& vec, const TrackError& err) const
  {
    // Get slope parameters and error matrix.

    double mu = vec(2);
    double mv = vec(3);
    double euu = err(2, 2);
    double evv = err(3, 3);
    double euv = err(3, 2);

    // Calculate maximal slope error.

    double dd = euu - evv;
    double delta = std::sqrt(dd*dd + 4.*euv*euv);
    double emax = std::sqrt(0.5 * (euu + evv + delta));

    // Convert slope error to angle error.

    double etheta = emax / (1. + mu*mu + mv*mv);

    // Done.

    return etheta;
  }

  /// Get starting error matrix for Kalman filter.
  ///
  /// Arguments:
  ///
  /// err - Error matrix.
  ///
  void SurfPlane::getStartingError(TrackError& err) const {
    err.resize(5, false);
    err.clear();
    err(0, 0) = 1000.;
    err(1, 1) = 1000.;
    err(2, 2) = 0.25;
    err(3, 3) = 0.25;
    err(4, 4) = 10.;
  }

} // end namespace trkf
