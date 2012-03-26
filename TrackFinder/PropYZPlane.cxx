///////////////////////////////////////////////////////////////////////
///
/// \file   PropYZPlane.cxx
///
/// \brief  Propagate between two SurfYZPlanes.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <cmath>
#include "TrackFinder/PropYZPlane.h"
#include "TrackFinder/SurfYZPlane.h"
#include "cetlib/exception.h"

namespace trkf {

  /// Default constructor.
  PropYZPlane::PropYZPlane()
  {}

  /// Destructor.
  PropYZPlane::~PropYZPlane()
  {}

  /// Propagate without error.
  /// Optionally return propagation matrix and noise matrix.
  ///
  /// Arguments:
  ///
  /// trk   - Track to propagate.
  /// psurf - Destination surface.
  /// dir   - Propagation direction (FORWARD, BACKWARD, or UNKNOWN).
  /// prop_matrix - Pointer to optional propagation matrix.
  /// noise_matrix - Pointer to optional noise matrix.
  ///
  /// Returned value: propagation distance + success flag.
  ///
  boost::optional<double> PropYZPlane::vec_prop(KTrack& trk,
						const boost::shared_ptr<const Surface>& psurf, 
						Propagator::PropDirection dir, 
						TrackMatrix* prop_matrix,
						TrackError* noise_matrix) const
  {
    // Set the default return value to be unitialized with value 0.

    boost::optional<double> result(false, 0.);

    // Cast the initial and destination surfaces to SurfYZPlane.
    // If either surface is the wrong type, return failure.

    // Get initial surface and surface parameters.

    const SurfYZPlane* from = dynamic_cast<const SurfYZPlane*>(&*trk.getSurface());
    if(from == 0)
      return result;
    double y01 = from->y0();
    double z01 = from->z0();
    double phi1 = from->phi();

    // Get destination surface and surface parameters.

    const SurfYZPlane* to = dynamic_cast<const SurfYZPlane*>(&*psurf);
    if(to == 0)
      return result;
    double y02 = to->y0();
    double z02 = to->z0();
    double phi2 = to->phi();

    // Calculate transcendental functions.

    double sinphi2 = std::sin(phi2);
    double cosphi2 = std::cos(phi2);
    double sindphi = std::sin(phi2 - phi1);
    double cosdphi = std::cos(phi2 - phi1);

    // Get the initial track state vector and track parameters.

    const TrackVector& vec = trk.getVector();
    if(vec.size() != 5)
      throw cet::exception("PropYZPlane") 
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

    // Calculate initial position and slope in coordinate system
    // of destination plane.

    double u2 = u1;
    double v2 = (y01 - y02) * cosphi2 + (z01 - z02) * sinphi2 + v1 * cosdphi;
    double w2 = -(y01 - y02) * sinphi2 + (z01 - z02) * cosphi2 - v1 * sindphi;

    // Calculate derivative dw2/dw1.
    // If dw2/dw1 == 0., that means the track is moving parallel
    // to destination plane.
    // In this case return propagation failure.

    double dw2dw1 = cosdphi - dvdw1 * sindphi;
    if(dw2dw1 == 0.)
      return result;

    double dudw2 = dudw1 / dw2dw1;
    double dvdw2 = (sindphi + dvdw1 * cosdphi) / dw2dw1;

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

    // Decide if propagation was in the right direction.
    // (Compare sign of s with requested direction).

    bool sok = (dir == Propagator::UNKNOWN ||
		(dir == Propagator::FORWARD && s >= 0.) ||
		(dir == Propagator::BACKWARD && s <= 0.));

    // Update result.

    result = boost::optional<double>(sok, s);		

    // If final result is failure (because of wrong direction),
    // return without updating the track or propagation matrix.

    if(!result)
      return result;

    // Update propagation matrix (if requested).

    if(prop_matrix != 0) {
      TrackMatrix& pm = *prop_matrix;
      pm.resize(vec.size(), vec.size(), false);

      // Calculate partial derivatives.

      pm(0,0) = 1.;   // du2/du1
      pm(1,0) = 0.;   // dv2/du1
      pm(2,0) = 0.;   // d(dudw2)/du1
      pm(3,0) = 0.;   // d(dvdw2)/du1
      pm(4,0) = 0.;   // d(pinv2)/du1

      pm(0,1) = dudw2 * sindphi;             // du2/dv1
      pm(1,1) = cosdphi + dvdw2 * sindphi;   // dv2/dv1
      pm(2,1) = 0.;                          // d(dudw2)/dv1
      pm(3,1) = 0.;                          // d(dvdw2)/dv1
      pm(4,1) = 0.;                          // d(pinv2)/dv1

      pm(0,2) = -w2 / dw2dw1;   // du2/d(dudw1);
      pm(1,2) = 0.;             // dv2/d(dudw1);
      pm(2,2) = 1. / dw2dw1;    // d(dudw2)/d(dudw1);
      pm(3,2) = 0.;             // d(dvdw2)/d(dudw1);
      pm(4,2) = 0.;             // d(pinv2)/d(dudw1);

      pm(0,3) = -w2 * dudw1 * sindphi / (dw2dw1*dw2dw1);   // du2/d(dvdw1);
      pm(1,3) = -w2 / (dw2dw1*dw2dw1);                     // dv2/d(dvdw1);
      pm(2,3) = dudw1 * sindphi / (dw2dw1*dw2dw1);         // d(dudw2)/d(dvdw1);
      pm(3,3) = 1. / (dw2dw1*dw2dw1);                      // d(dvdw2)/d(dvdw1);
      pm(4,3) = 0.;                                        // d(pinv2)/d(dvdw1);

      pm(0,4) = 0.;   // du2/d(pinv1);
      pm(1,4) = 0.;   // dv2/d(pinv1);
      pm(2,4) = 0.;   // d(dudw2)/d(pinv1);
      pm(3,4) = 0.;   // d(dvdw2)/d(pinv1);
      pm(4,4) = 1.;   // d(pinv2)/d(pinv1);
    }

    // Update noise matrix (if requested).

    if(noise_matrix != 0) {
      noise_matrix->resize(vec.size(), vec.size(), false);
      noise_matrix->clear();
    }

    // Construct track vector at destination surface.

    TrackVector vec2(vec.size());
    vec2(0) = u2p;
    vec2(1) = v2p;
    vec2(2) = dudw2;
    vec2(3) = dvdw2;
    vec2(4) = pinv;

    // Update track.

    trk.setSurface(psurf);
    trk.setVector(vec2);
    trk.setDirection(dir2);

    // Done.

    return result;
  }
} // end namespace trkf
