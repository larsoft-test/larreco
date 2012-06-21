///////////////////////////////////////////////////////////////////////
///
/// \file   PropXYZPlane.cxx
///
/// \brief  Propagate between two SurfXYZPlane or SurfYZPlane surfaces.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <cmath>
#include "TrackFinder/PropXYZPlane.h"
#include "TrackFinder/SurfYZPlane.h"
#include "TrackFinder/SurfXYZPlane.h"
#include "TrackFinder/InteractPlane.h"
#include "cetlib/exception.h"

namespace trkf {

  /// Constructor.
  ///
  /// Arguments.
  ///
  /// tcut   - Delta ray energy cutoff for calculating dE/dx.
  ///
  PropXYZPlane::PropXYZPlane(double tcut) :
    Propagator(tcut, boost::shared_ptr<const Interactor>(new InteractPlane(tcut)))
  {}

  /// Destructor.
  PropXYZPlane::~PropXYZPlane()
  {}

  /// Propagate without error.
  /// Optionally return propagation matrix and noise matrix.
  ///
  /// Arguments:
  ///
  /// trk   - Track to propagate.
  /// psurf - Destination surface.
  /// dir   - Propagation direction (FORWARD, BACKWARD, or UNKNOWN).
  /// doDedx - dE/dx enable/disable flag.
  /// prop_matrix - Pointer to optional propagation matrix.
  /// noise_matrix - Pointer to optional noise matrix.
  ///
  /// Returned value: propagation distance + success flag.
  ///
  boost::optional<double>
  PropXYZPlane::short_vec_prop(KTrack& trk,
			       const boost::shared_ptr<const Surface>& psurf, 
			       Propagator::PropDirection dir,
			       bool doDedx,
			       TrackMatrix* prop_matrix,
			       TrackError* noise_matrix) const
  {
    // Set the default return value to be unitialized with value 0.

    boost::optional<double> result(false, 0.);

    // Cast the initial and destination surfaces to SurfXYZPlane or
    // SurfYZPlane.  If either surface is the wrong type, return failure.

    // Get initial surface and surface parameters.

    double x01 = 0.;
    double y01 = 0.;
    double z01 = 0.;
    double theta1 = 0.;
    double phi1 = 0.;
    const SurfXYZPlane* from = dynamic_cast<const SurfXYZPlane*>(&*trk.getSurface());
    if(from != 0) {
      x01 = from->x0();
      y01 = from->y0();
      z01 = from->z0();
      theta1 = from->theta();
      phi1 = from->phi();
    }
    else {
      const SurfYZPlane* from = dynamic_cast<const SurfYZPlane*>(&*trk.getSurface());
      if(from != 0) {

	// x01=0 and theta1=0 for YZ plane.

	y01 = from->y0();
	z01 = from->z0();
	phi1 = from->phi();
      }
      else
	return result;
    }

    // Get destination surface and surface parameters.

    double x02 = 0.;
    double y02 = 0.;
    double z02 = 0.;
    double theta2 = 0.;
    double phi2 = 0.;
    const SurfXYZPlane* to = dynamic_cast<const SurfXYZPlane*>(&*psurf);
    if(to != 0) {
      x02 = to->x0();
      y02 = to->y0();
      z02 = to->z0();
      theta2 = to->theta();
      phi2 = to->phi();
    }
    else {
      const SurfYZPlane* to = dynamic_cast<const SurfYZPlane*>(&*psurf);
      if(to != 0) {

	// x02=0 and theta2=0 for YZ plane.

	y02 = from->y0();
	z02 = from->z0();
	phi2 = from->phi();
      }
      else
	return result;
    }

    // Calculate transcendental functions.

    double sinth1 = std::sin(theta1);
    double costh1 = std::cos(theta1);
    double sinth2 = std::sin(theta2);
    double costh2 = std::cos(theta2);

    double sinphi2 = std::sin(phi2);
    double cosphi2 = std::cos(phi2);
    double sindphi = std::sin(phi2 - phi1);
    double cosdphi = std::cos(phi2 - phi1);

    // Get the initial track state vector and track parameters.

    const TrackVector& vec = trk.getVector();
    if(vec.size() != 5)
      throw cet::exception("PropXYZPlane") 
	<< "Track state vector has wrong size" << vec.size() << "\n";
    double u1 = vec(0);
    double v1 = vec(1);
    double dudw1 = vec(2);
    double dvdw1 = vec(3);
    double pinv = vec(4);
    Surface::TrackDirection dir1 = trk.getDirection();

    // Make sure initial track has a valid direction.

    if(dir1 == Surface::UNKNOWN)
      return result;
    assert(dir1 == Surface::FORWARD || dir1 == Surface::BACKWARD);

    // Calculate elements of rotation matrix from initial coordinate
    // system to destination coordinte system.

    double ruu = costh1*costh2 + sinth1*sinth2*cosdphi;
    double ruv = sinth2*sindphi;
    double ruw = sinth1*costh2 - costh1*sinth2*cosdphi;

    double rvu = -sinth1*sindphi;
    double rvv = cosdphi;
    double rvw = costh1*sindphi;

    double rwu = costh1*sinth2 - sinth1*costh2*cosdphi;
    double rwv = -costh2*sindphi;
    double rww = sinth1*sinth2 + costh1*costh2*cosdphi;

    // Calculate elements of rotation matrix from global coordinate
    // system to destination coordinate system.

    double rux = costh2;
    double ruy = sinth2*sinphi2;
    double ruz = -sinth2*cosphi2;

    double rvx = 0.;
    double rvy = cosphi2;
    double rvz = sinphi2;

    double rwx = sinth2;
    double rwy = -costh2*sinphi2;
    double rwz = costh2*cosphi2;

    // Calculate the initial position in the destination coordinate
    // system.

    double u2 = (x01-x02)*rux + (y01-y02)*ruy + (z01-z02)*ruz + u1*ruu + v1*ruv;
    double v2 = (x01-x02)*rvx + (y01-y02)*rvy + (z01-z02)*rvz + u1*rvu + v1*rvv;
    double w2 = (x01-x02)*rwx + (y01-y02)*rwy + (z01-z02)*rwz + u1*rwu + v1*rwv;

    // Calculate the derivative dw2/dw1;
    // If dw2/dw1 == 0., that means the track is moving parallel
    // to destination plane.
    // In this case return propagation failure.

    double dw2dw1 = dudw1*rwu + dvdw1*rwv + rww;
    if(dw2dw1 == 0.)
      return result;

    // Calculate slope in destination plane coordinates.

    double dudw2 = (dudw1*ruu + dvdw1*ruv + ruw) / dw2dw1;
    double dvdw2 = (dudw1*rvu + dvdw1*rvv + rvw) / dw2dw1;

    // Calculate position at destination surface (propagate distance -w2).

    double u2p = u2 - w2 * dudw2;
    double v2p = v2 - w2 * dvdw2;

    // Calculate direction parameter at destination surface.
    // Direction will flip if dw2dw1 < 0.;

    Surface::TrackDirection dir2;
    if(dir1 == Surface::FORWARD) {
      if(dw2dw1 > 0.)
    	dir2 = Surface::FORWARD;
      else
	dir2 = Surface::BACKWARD;
    }
    else {
      if(dw2dw1 > 0.)
    	dir2 = Surface::BACKWARD;
      else
    	dir2 = Surface::FORWARD;
    }
    assert(dir2 == Surface::FORWARD || dir2 == Surface::BACKWARD);

    // Calculate the signed propagation distance.

    double s = -w2 * std::sqrt(1. + dudw2*dudw2 + dvdw2*dvdw2);
    if(dir2 == Surface::BACKWARD)
      s = -s;

    // Check if propagation was in the right direction.
    // (Compare sign of s with requested direction).

    bool sok = (dir == Propagator::UNKNOWN ||
		(dir == Propagator::FORWARD && s >= 0.) ||
		(dir == Propagator::BACKWARD && s <= 0.));

    // If wrong direction, return failure without updating the track
    // or propagation matrix.

    if(!sok)
      return result;

    // Find final momentum.

    double deriv = 1.;
    boost::optional<double> pinv2(true, pinv);
    if(doDedx && s != 0.) {
      double* pderiv = (prop_matrix != 0 ? &deriv : 0);
      pinv2 = dedx_prop(pinv, trk.Mass(), s, pderiv);
    }

    // Return failure in case of range out.

    if(!pinv2)
      return result;

    // Update result object (success guaranteed).

    result = boost::optional<double>(true, s);		

    // Update propagation matrix (if requested).

    if(prop_matrix != 0) {
      TrackMatrix& pm = *prop_matrix;
      pm.resize(vec.size(), vec.size(), false);

      // Calculate partial derivatives.

      pm(0,0) = ruu - dudw2*rwu;    // du2/du1
      pm(1,0) = rvu - dvdw2*rwu;    // dv2/du1
      pm(2,0) = 0.;                 // d(dudw2)/du1
      pm(3,0) = 0.;                 // d(dvdw2)/du1
      pm(4,0) = 0.;                 // d(pinv2)/du1

      pm(0,1) = ruv - dudw2*rwv;    // du2/dv1
      pm(1,1) = rvv - dvdw2*rwv;    // dv2/dv1
      pm(2,1) = 0.;                 // d(dudw2)/dv1
      pm(3,1) = 0.;                 // d(dvdw2)/dv1
      pm(4,1) = 0.;                 // d(pinv2)/dv1

      pm(0,2) = -w2 * (ruu - dudw2*rwu) / dw2dw1;    // du2/d(dudw1);
      pm(1,2) = -w2 * (rvu - dvdw2*rwu) / dw2dw1;    // dv2/d(dudw1);
      pm(2,2) = (ruu - dudw2*rwu) / dw2dw1;          // d(dudw2)/d(dudw1);
      pm(3,2) = (rvu - dvdw2*rwu) / dw2dw1;          // d(dvdw2)/d(dudw1);
      pm(4,2) = 0.;                                  // d(pinv2)/d(dudw1);

      pm(0,3) = -w2 * (ruv - dudw2*rwv) / dw2dw1;    // du2/d(dvdw1);
      pm(1,3) = -w2 * (rvv - dvdw2*rwv) / dw2dw1;    // dv2/d(dvdw1);
      pm(2,3) = (ruv - dudw2*rwv) / dw2dw1;          // d(dudw2)/d(dvdw1);
      pm(3,3) = (rvv - dvdw2*rwv) / dw2dw1;          // d(dvdw2)/d(dvdw1);
      pm(4,3) = 0.;                                  // d(pinv2)/d(dvdw1);

      pm(0,4) = 0.;      // du2/d(pinv1);
      pm(1,4) = 0.;      // dv2/d(pinv1);
      pm(2,4) = 0.;      // d(dudw2)/d(pinv1);
      pm(3,4) = 0.;      // d(dvdw2)/d(pinv1);
      pm(4,4) = deriv;   // d(pinv2)/d(pinv1);
    }

    // Update noise matrix (if requested).

    if(noise_matrix != 0) {
      noise_matrix->resize(vec.size(), vec.size(), false);
      if(getInteractor().get() != 0)
	getInteractor()->noise(trk, s, *noise_matrix);
      else
	noise_matrix->clear();
    }

    // Construct track vector at destination surface.

    TrackVector vec2(vec.size());
    vec2(0) = u2p;
    vec2(1) = v2p;
    vec2(2) = dudw2;
    vec2(3) = dvdw2;
    vec2(4) = *pinv2;

    // Update track.

    trk.setSurface(psurf);
    trk.setVector(vec2);
    trk.setDirection(dir2);

    // Done.

    return result;
  }

} // end namespace trkf
