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
/// 1.  Propagate without error (method vec_prop).
/// 2.  Linearized propagate without error (method lin_prop).
/// 3.  Propagate with error, but without noise (method err_prop).
/// 4.  Propagate with error and noise (method noise_prop).
///
/// The second, third, and fourth use cases are implemented locally
/// within this class in terms of the first use case (propagate
/// without error), which is defined by a single pure virtual method,
/// vec_prop.  Linearized propagation uses a linear approximation of
/// the propagation function about some reference trajectory.  Method
/// vec_prop includes optional hooks for returning the propagation
/// matrix and the noise matrix.  The propgation matrix and noise
/// matrix provide enough information to update the track error matrix
/// locally within this class.  Therefore, methods err_prop and
/// noise_prop are not virtual.
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

#include "TrackFinder/KalmanLinearAlgebra.h"
#include "TrackFinder/KETrack.h"
#include "TrackFinder/Interactor.h"
#include "boost/shared_ptr.hpp"
#include "boost/optional.hpp"

namespace trkf {

  class Propagator
  {
  public:

    /// Propagation direction enum.
    enum PropDirection {FORWARD, BACKWARD, UNKNOWN};

    /// Constructor.
    Propagator(double tcut, const boost::shared_ptr<const Interactor>& interactor);

    /// Destructor.
    virtual ~Propagator();

    // Accessors.

    double getTcut() const {return fTcut;}
    const boost::shared_ptr<const Interactor>& getInteractor() const {return fInteractor;}

    // Virtual methods.

    /// Clone method.
    virtual Propagator* clone() const = 0;

    /// Propagate without error.
    virtual boost::optional<double> vec_prop(KTrack& trk,
					     const boost::shared_ptr<const Surface>& psurf, 
					     PropDirection dir,
					     bool doDedx,
					     TrackMatrix* prop_matrix = 0,
					     TrackError* noise_matrix = 0) const = 0;

    // Error and noise propagation methods are not virtual.

    /// Linearized propagate without error.
    virtual boost::optional<double> lin_prop(KTrack& trk,
					     const boost::shared_ptr<const Surface>& psurf, 
					     PropDirection dir,
					     bool doDedx,
					     KTrack* ref = 0,
					     TrackMatrix* prop_matrix = 0,
					     TrackError* noise_matrix = 0) const;
    /// Propagate with error, but without noise.
    boost::optional<double> err_prop(KETrack& tre,
				     const boost::shared_ptr<const Surface>& psurf, 
				     PropDirection dir,
				     bool doDedx,
				     KTrack* ref = 0,
				     TrackMatrix* prop_matrix = 0) const;

    /// Propagate with error and noise.
    boost::optional<double> noise_prop(KETrack& tre,
				       const boost::shared_ptr<const Surface>& psurf, 
				       PropDirection dir,
				       bool doDedx,
				       KTrack* ref = 0) const;

    /// Method to calculate updated momentum due to dE/dx.
    boost::optional<double> dedx_prop(double pinv, double mass,
				      double s, double* deriv=0) const;

  private:

    // Attributes.

    double fTcut;                                     ///< Maximum delta ray energy for dE/dx.
    boost::shared_ptr<const Interactor> fInteractor;  ///< Interactor (for calculating noise).
  };
}

#endif
