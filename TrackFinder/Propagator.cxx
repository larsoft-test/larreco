///////////////////////////////////////////////////////////////////////
///
/// \file   Propagator.cxx
///
/// \brief  Base class for Kalman filter propagator.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "TrackFinder/Propagator.h"

namespace trkf {

  /// Default constructor.
  Propagator::Propagator()
  {}

  /// Destructor.
  Propagator::~Propagator()
  {}

  /// Propagate with error, but without noise (i.e. reversibly).
  ///
  /// Arguments:
  ///
  /// tre   - Track to propagate.
  /// psurf - Destination surface.
  /// dir   - Propagation direction (FORWARD, BACKWARD, or UNKNOWN).
  ///
  /// Returned value: propagation distance + success flag.
  ///
  boost::optional<double> Propagator::err_prop(KETrack& tre,
					       const boost::shared_ptr<const Surface>& psurf, 
					       PropDirection dir) const
  {
    // Propagate without error, get propagation matrix.

    TrackMatrix prop_matrix;
    boost::optional<double> result = vec_prop(tre, psurf, dir, &prop_matrix);

    // If propagation succeeded, update track error matrix.

    if(!!result) {
      TrackMatrix temp = prod(tre.getError(), trans(prop_matrix));
      TrackError newerr = prod(prop_matrix, temp);
      tre.setError(newerr);
    }

    // Done.

    return result;
  }

  /// Propagate with error and noise.
  ///
  /// Arguments:
  ///
  /// tre   - Track to propagate.
  /// psurf - Destination surface.
  /// dir   - Propagation direction (FORWARD, BACKWARD, or UNKNOWN).
  ///
  /// Returned value: propagation distance + success flag.
  ///
  boost::optional<double> Propagator::noise_prop(KETrack& tre,
						 const boost::shared_ptr<const Surface>& psurf, 
						 PropDirection dir) const
  {
    // Propagate without error, get propagation matrix and noise matrix.

    TrackMatrix prop_matrix;
    TrackError noise_matrix;
    boost::optional<double> result = vec_prop(tre, psurf, dir, &prop_matrix, &noise_matrix);

    // If propagation succeeded, update track error matrix.

    if(!!result) {
      TrackMatrix temp = prod(tre.getError(), trans(prop_matrix));
      TrackError newerr = prod(prop_matrix, temp);
      newerr += noise_matrix;
      tre.setError(newerr);
    }

    // Done.

    return result;
  }
} // end namespace trkf
