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
#include "Utilities/LArProperties.h"

namespace trkf {

  /// Constructor.
  ///
  /// Arguments:
  ///
  /// tcut - Maximum delta ray energy.
  ///
  Propagator::Propagator(double tcut, const boost::shared_ptr<const Interactor>& interactor) :
    fTcut(tcut),
    fInteractor(interactor)
  {}

  /// Destructor.
  Propagator::~Propagator()
  {}

  /// Propagate with error, but without noise (i.e. reversibly).
  ///
  /// Arguments:
  ///
  /// tre    - Track to propagate.
  /// psurf  - Destination surface.
  /// dir    - Propagation direction (FORWARD, BACKWARD, or UNKNOWN).
  /// doDedx - dE/dx enable/disable flag.
  ///
  /// Returned value: propagation distance + success flag.
  ///
  boost::optional<double> Propagator::err_prop(KETrack& tre,
					       const boost::shared_ptr<const Surface>& psurf, 
					       PropDirection dir,
					       bool doDedx,
					       TrackMatrix* prop_matrix) const
  {
    // Propagate without error, get propagation matrix.

    TrackMatrix prop_temp;
    if(prop_matrix == 0)
      prop_matrix = &prop_temp;
    boost::optional<double> result = vec_prop(tre, psurf, dir, doDedx, prop_matrix, 0);

    // If propagation succeeded, update track error matrix.

    if(!!result) {
      TrackMatrix temp = prod(tre.getError(), trans(*prop_matrix));
      TrackError newerr = prod(*prop_matrix, temp);
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
  /// doDedx - dE/dx enable/disable flag.
  ///
  /// Returned value: propagation distance + success flag.
  ///
  boost::optional<double> Propagator::noise_prop(KETrack& tre,
						 const boost::shared_ptr<const Surface>& psurf, 
						 PropDirection dir,
						 bool doDedx) const
  {
    // Propagate without error, get propagation matrix and noise matrix.

    TrackMatrix prop_matrix;
    TrackError noise_matrix;
    boost::optional<double> result = vec_prop(tre, psurf, dir, doDedx, 
					      &prop_matrix, &noise_matrix);

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

  /// Method to calculate updated momentum due to dE/dx.
  ///
  /// Arguments:
  ///
  /// pinv  - Initial inverse momentum (units c/GeV).
  /// mass  - Particle mass (GeV/c^2).
  /// s     - Path distance.
  /// deriv - Pointer to store derivative d(pinv2)/d(pinv1) if nonzero.
  ///
  /// Returns: Final inverse momentum (pinv2) + success flag.
  ///
  /// Failure is returned in case of range out.
  ///
  /// Inverse momentum can be signed (q/p).  Returned inverse momentum
  /// has the same sign as the input.
  ///
  /// In this method, we are solving the differential equation in
  /// terms of energy.
  ///
  /// dE/dx = -f(E)
  ///
  /// where f(E) is the stopping power returned by method
  /// LArProperties::Eloss.
  ///
  /// We expect that this method will be called exclusively for short
  /// distance propagation.  The differential equation is solved using
  /// the midpoint method using a single step, which requires two
  /// evaluations of f(E).
  ///
  /// dE = -s*f(E1)
  /// E2 = E1 - s*f(E1 + 0.5*dE)
  ///
  /// The derivative is calculated assuming E2 = E1 + constant, giving
  ///
  /// d(pinv2)/d(pinv1) = pinv2^3 E2 / (pinv1^3 E1).
  /// 
  ///
  boost::optional<double> Propagator::dedx_prop(double pinv, double mass,
						double s, double* deriv) const
  {
    // Set the default return value to be unitialized with value 0.

    boost::optional<double> result(false, 0.);

    // Get LAr service.

    art::ServiceHandle<util::LArProperties> larprop;

    // Calculate final energy.

    double p1 = 1./std::abs(pinv);
    double e1 = std::sqrt(p1*p1 + mass*mass);
    double de = -0.001 * s * larprop->Eloss(p1, mass, fTcut);
    double emid = e1 + 0.5 * de;
    if(emid > mass) {
      double pmid = std::sqrt(emid*emid - mass*mass);
      double e2 = e1 - 0.001 * s * larprop->Eloss(pmid, mass, fTcut);
      if(e2 > mass) {
	double p2 = std::sqrt(e2*e2 - mass*mass);
	double pinv2 = 1./p2;
	if(pinv < 0.)
	  pinv2 = -pinv2;

	// Calculation was successful, update result.

	result = boost::optional<double>(true, pinv2);

	// Also calculate derivative, if requested.

	if(deriv != 0)
	  *deriv = pinv2*pinv2*pinv2 * e2 / (pinv*pinv*pinv * e1);
      }
    }

    // Done.

    return result;
  }
} // end namespace trkf
