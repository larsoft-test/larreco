//
// File: PropTest.cxx
//
// Purpose: Single source executable with tests for propagator classes.
//

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <vector>
#include "TrackFinder/KETrack.h"
#include "TrackFinder/SurfYZPlane.h"
#include "TrackFinder/PropYZPlane.h"
#include "TMath.h"

int main()
{
  // Make sure assert is enabled.

  bool assert_flag = false;
  assert((assert_flag = true, assert_flag));
  if ( ! assert_flag ) {
    std::cerr << "Assert is disabled" << std::endl;
    return 1;
  }

  // Make a PropYZPlane propagate to test.

  const trkf::Propagator* prop = new trkf::PropYZPlane;

  // Make some random surfaces.
  // Also make initial tracks.

  std::vector<boost::shared_ptr<const trkf::Surface> > surfaces;
  std::vector<trkf::KETrack> tracks;

  int nsurf = 10;
  for(int isurf = 0; isurf < nsurf; ++isurf) {

    // Make random surface.

    double y0 = 100.*double(rand()) / double(RAND_MAX) - 50.;  // (-50,50)
    double z0 = 1000.*double(rand()) / double(RAND_MAX);       // (0,1000)
    double phi = TMath::TwoPi() * double(rand()) / double(RAND_MAX) - TMath::Pi();  // (-pi,pi)
    boost::shared_ptr<const trkf::Surface> psurf(new trkf::SurfYZPlane(y0, z0, phi));
    surfaces.push_back(psurf);

    // Make random track vector.

    double u = 100.*double(rand()) / double(RAND_MAX);  // (0,100)
    double v = 100.*double(rand()) / double(RAND_MAX) - 50.;  // (-50, 50)
    double dudw = 4.*double(rand()) / double(RAND_MAX) - 2.;  // (-2, 2)
    double dvdw = 4.*double(rand()) / double(RAND_MAX) - 2.;  // (-2, 2)
    double pinv = 0.9*double(rand()) / double(RAND_MAX) + 0.1;  // (0.1, 1.0)
    trkf::TrackVector vec(5);
    vec(0) = u;
    vec(1) = v;
    vec(2) = dudw;
    vec(3) = dvdw;
    vec(4) = pinv;

    // Make initial error matrix (diagonal, not random).

    trkf::TrackError err(5);
    err.clear();
    err(0,0) = 0.1;
    err(1,1) = 0.1;
    err(2,2) = 0.01;
    err(3,3) = 0.01;
    err(4,4) = 0.01;

    // Make random track direction.

    trkf::Surface::TrackDirection dir = trkf::Surface::FORWARD;
    if(rand() % 2 == 0)
      dir = trkf::Surface::BACKWARD;

    // Make KETrack.

    tracks.push_back(trkf::KETrack(psurf, vec, err, dir));    
  }

  // Loop over initial surface/track.

  for(int isurf = 1; isurf < nsurf; ++isurf) {

    std::cout << "\nInitial track " << isurf << std::endl;

    // Get initial surface and track object.

    const boost::shared_ptr<const trkf::Surface>& psurf1 = surfaces[isurf];
    const trkf::KETrack& trk1 = tracks[isurf];
    const trkf::TrackVector& vec1 = trk1.getVector();
    const trkf::TrackError& err1 = trk1.getError();

    // Get initial position and momentum vector.

    double xyz1[3];
    double mom1[3];
    trk1.getPosition(xyz1);
    trk1.getMomentum(mom1);
    std::cout << "Initial position = "
	      << xyz1[0] << ", " << xyz1[1] << ", " << xyz1[2] << std::endl;
    std::cout << "Initial momentum = "
	      << mom1[0] << ", " << mom1[1] << ", " << mom1[2] << std::endl;

    // Loop over destination surface.

    for(int jsurf = 0; jsurf < isurf; ++jsurf) {

      std::cout << "Destination " << jsurf << std::endl;

      // Get destination surface.

      const boost::shared_ptr<const trkf::Surface>& psurf2 = surfaces[jsurf];

      // Make a copy of the initial track and propagate it to the 
      // destination surface.

      trkf::KETrack trk2 = trk1;
      boost::optional<double> dist12 = prop->err_prop(trk2, psurf2);
      if(!!dist12)
	std::cout << "Propagation distance = " << *dist12 << std::endl;
      else
	assert(false);

      // Get final position and momentum vector.

      double xyz2[3];
      double mom2[3];
      trk2.getPosition(xyz2);
      trk2.getMomentum(mom2);
      std::cout << "Final position = "
		<< xyz2[0] << ", " << xyz2[1] << ", " << xyz2[2] << std::endl;
      std::cout << "Final momentum = "
		<< mom2[0] << ", " << mom2[1] << ", " << mom2[2] << std::endl;

      // Calculate displacement vector.

      double dx[3];
      dx[0] = xyz2[0] - xyz1[0];
      dx[1] = xyz2[1] - xyz1[1];
      dx[2] = xyz2[2] - xyz1[2];

      // Check the sign of the propagation distance.
      // The sign of the propatation distance should be the same as 
      // the sign of the dot product of the initial momentum and the
      // displacement vector.

      double pdotdx = mom1[0]*dx[0] + mom1[1]*dx[1] + mom1[2]*dx[2];
      assert(pdotdx * (*dist12) > 0.);

      // Check propagation distance against displacement vector.

      double dist = std::sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      std::cout << "Displacement = " << dist << std::endl;
      assert(std::abs(dist - std::abs(*dist12)) <= 1.e-10);

      // Make sure final momentum is parallel to initial momentum.

      double m0 = mom1[1]*mom2[2] - mom1[2]*mom2[1];
      double m1 = mom1[2]*mom2[0] - mom1[0]*mom2[2];
      double m2 = mom1[0]*mom2[1] - mom1[1]*mom2[0];
      std::cout << "mom1 x mom2 = " << m0 << ", " << m1 << ", " << m2 << std::endl;
      assert(std::abs(m0) <= 1.e-8);
      assert(std::abs(m1) <= 1.e-8);
      assert(std::abs(m2) <= 1.e-8);

      // Make sure displacement vector is parallel to initial momentum.

      double d0 = mom1[1]*dx[2] - mom1[2]*dx[1];
      double d1 = mom1[2]*dx[0] - mom1[0]*dx[2];
      double d2 = mom1[0]*dx[1] - mom1[1]*dx[0];
      std::cout << "mom1 x dx = " << d0 << ", " << d1 << ", " << d2 << std::endl;
      assert(std::abs(d0) <= 1.e-8);
      assert(std::abs(d1) <= 1.e-8);
      assert(std::abs(d2) <= 1.e-8);

      // Test propagation matrix by numerical partial derivatives.

      trkf::TrackMatrix pm(vec1.size(), vec1.size());
      trkf::KTrack trk10 = trk1;
      boost::optional<double> stat =
	prop->vec_prop(trk10, psurf2, trkf::Propagator::UNKNOWN, &pm);
      assert(!!stat);

      double small = 1.e-5;

      for(unsigned int i = 0; i < vec1.size(); ++i) {
	for(unsigned int j = 0; j < vec1.size(); ++j) {

	  // Calculate d(vec[i])/d(vec[j])

	  trkf::KTrack trk1a = trk1;
	  trkf::TrackVector vec1a = vec1;
	  vec1a(j) = vec1(j) - small;
	  trk1a.setVector(vec1a);
	  boost::optional<double> stata = prop->vec_prop(trk1a, psurf2);
	  assert(!!stata);

	  trkf::KTrack trk1b = trk1;
	  trkf::TrackVector vec1b = vec1;
	  vec1b(j) = vec1(j) + small;
	  trk1b.setVector(vec1b);
	  boost::optional<double> statb = prop->vec_prop(trk1b, psurf2);
	  assert(!!statb);

	  // Compare numerical and analytic partial derivative.

	  double dij = (trk1b.getVector()(i) - trk1a.getVector()(i)) / (2.*small);
	  //std::cout << "(" << i << "," << j << "): " << dij << ", " << pm(i,j) << std::endl;
	  assert(std::abs(dij - pm(i,j)) <= 1.e-5*std::abs(dij));
	}
      }

      // Now propagate back to the original surface.

      boost::optional<double> dist21 = prop->err_prop(trk2, psurf1);
      assert(!!dist21);
      assert(std::abs(dist - std::abs(*dist21)) <= 1.e-10);

      // Check that state vector and error matrix returned to the original.
      // This will test that the forward and backward propagation matrices
      // are inverses.

      const trkf::TrackVector& vec2 = trk2.getVector();
      const trkf::TrackError& err2 = trk2.getError();
      int n = vec1.size();
      for(int i=0; i<n; ++i) {
	assert(std::abs(vec1(i) - vec2(i)) <= 1.e-8);
	for(int j=0; j<n; ++j)
	  assert(std::abs(err1(i,j) - err2(i,j)) <= 1.e-8);
      }
    }
  }
  
  // Done (success).

  std::cout << "PropTest: All tests passed." << std::endl;

  return 0;
}
