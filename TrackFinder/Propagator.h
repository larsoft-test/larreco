////////////////////////////////////////////////////////////////////////
///
/// \file   Propagator.h
///
/// \brief  Base class for Kalman filter track propagator.
///
/// \author H. Greenlee 
///
/// This class provides the general interface for propagating a track
/// (KTrack or KETrack) from its current surface to some destionation
/// surface.
///
/// This class supports various use cases.
///
/// 1.  Propagate without error, short distance (method short_vec_prop).
/// 2.  Propagate without error, long distance (method vec_prop).
/// 3.  Linearized propagation without error, short distance (method lin_prop).
/// 4.  Propagate with error, but without noise (method err_prop).
/// 5.  Propagate with error and noise (method noise_prop).
///
///
/// All methods except short_vec_prop are implemented locally within
/// this class by calling short_vec_prop, which is pure virtual.
///
/// The long distance propagation method (vec_prop) divides
/// propagation into steps if the distance exceeds some threshold
///
/// Linearized propagation uses a linear approximation of the
/// propagation function with respect to some reference trajectory.
/// This type of propagation only makes sense for short distance
/// propagation, so lin_prop is implemented by calling short_vec_prop.
///
/// All of the *vec_prop methods include optional hooks for returning
/// the propagation matrix and the noise matrix.  These hooks provide
/// enough information to propagate the error matrix (methods err_prop
/// and noise_prop) locally within this class in terms of method
/// vec_prop (so long distance propagation with error is supported).
///
/// All propagation methods update the surface and track state vector,
/// provided the propagation is successful.  The error and noise
/// propagation methods additionally update the track error matrix.
///
/// Use case three (propagate with error, but without noise) updates the
/// track error matrix reversibly.
///
/// Use case four (propagate with error and noise) adds irreversible
/// propagation noise to the error matrix.
///
/// All propagation methods allow the direction of propagation to be 
/// specified as forward, backward, or unknown.  If the direction is 
/// specified as unknown, the propagator decides which direction to use.
///
/// All propagation methods return the propagation distance as a value
/// of type boost::optional<double>.  This type of value is equivalent
/// to the contained value (that is, the propagation distance), plus a
/// flag indicating whether the contained value is initialized.  A
/// non-initialized return value means that the propagation failed.
///
////////////////////////////////////////////////////////////////////////

#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include <memory>
#include "TrackFinder/KalmanLinearAlgebra.h"
#include "TrackFinder/KETrack.h"
#include "TrackFinder/Interactor.h"
#include "boost/optional.hpp"

namespace trkf {

  class Propagator
  {
  public:

    /// Propagation direction enum.
    enum PropDirection {FORWARD, BACKWARD, UNKNOWN};

    /// Constructor.
    Propagator(double tcut, const std::shared_ptr<const Interactor>& interactor);

    /// Destructor.
    virtual ~Propagator();

    // Accessors.

    double getTcut() const {return fTcut;}
    const std::shared_ptr<const Interactor>& getInteractor() const {return fInteractor;}

    // Virtual methods.

    /// Clone method.
    virtual Propagator* clone() const = 0;

    /// Propagate without error (short distance).
    virtual boost::optional<double> short_vec_prop(KTrack& trk,
						   const std::shared_ptr<const Surface>& psurf, 
						   PropDirection dir,
						   bool doDedx,
						   TrackMatrix* prop_matrix = 0,
						   TrackError* noise_matrix = 0) const = 0;

    /// Propagate without error (long distance).
    boost::optional<double> vec_prop(KTrack& trk,
				     const std::shared_ptr<const Surface>& psurf, 
				     PropDirection dir,
				     bool doDedx,
				     TrackMatrix* prop_matrix = 0,
				     TrackError* noise_matrix = 0) const;

    /// Linearized propagate without error.
    boost::optional<double> lin_prop(KTrack& trk,
				     const std::shared_ptr<const Surface>& psurf, 
				     PropDirection dir,
				     bool doDedx,
				     KTrack* ref = 0,
				     TrackMatrix* prop_matrix = 0,
				     TrackError* noise_matrix = 0) const;
    /// Propagate with error, but without noise.
    boost::optional<double> err_prop(KETrack& tre,
				     const std::shared_ptr<const Surface>& psurf, 
				     PropDirection dir,
				     bool doDedx,
				     KTrack* ref = 0,
				     TrackMatrix* prop_matrix = 0) const;

    /// Propagate with error and noise.
    boost::optional<double> noise_prop(KETrack& tre,
				       const std::shared_ptr<const Surface>& psurf, 
				       PropDirection dir,
				       bool doDedx,
				       KTrack* ref = 0) const;

    /// Method to calculate updated momentum due to dE/dx.
    boost::optional<double> dedx_prop(double pinv, double mass,
				      double s, double* deriv=0) const;

  private:

    // Attributes.

    double fTcut;                                   ///< Maximum delta ray energy for dE/dx.
    std::shared_ptr<const Interactor> fInteractor;  ///< Interactor (for calculating noise).
  };
}

#endif
