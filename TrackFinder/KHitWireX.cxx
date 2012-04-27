///////////////////////////////////////////////////////////////////////
///
/// \file   KHitWireX.cxx
///
/// \brief  Kalman filter wire-time measurement on a SurfWireX surface.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "TrackFinder/KHitWireX.h"
#include "TrackFinder/SurfWireX.h"
#include "Geometry/geo.h"

namespace trkf {

  /// Constructor from channel.
  ///
  /// Arguments:
  ///
  /// channel - Readout channel number.
  /// x       - X position.
  /// xerr    - X error.
  ///
  KHitWireX::KHitWireX(unsigned int channel, double x, double xerr) :
    KHit(boost::shared_ptr<const Surface>(new SurfWireX(channel)))
  {
    trkf::KVector<1>::type mvec(1, x);
    setMeasVector(mvec);

    trkf::KSymMatrix<1>::type merr(1);
    merr(0,0) = xerr * xerr;
    setMeasError(merr);

    // Get geometry service.

    art::ServiceHandle<geo::Geometry> geom;

    // Get wire geometry.

    unsigned int cstat, tpc, plane, wire;
    geom->ChannelToWire(channel, cstat, tpc, plane, wire);
    setMeasPlane(plane);
  }

  /// Constructor from surface.
  ///
  /// Arguments:
  ///
  /// psurf - Surface.
  /// x     - X position.
  /// xerr  - X error.
  /// plane - Plane index.
  ///
  KHitWireX::KHitWireX(const boost::shared_ptr<const Surface>& psurf,
		       double x, double xerr, int plane) :
    KHit(psurf)
  {
    trkf::KVector<1>::type mvec(1, x);
    setMeasVector(mvec);

    trkf::KSymMatrix<1>::type merr(1);
    merr(0,0) = xerr * xerr;
    setMeasError(merr);
    setMeasPlane(plane);
  }

  /// Destructor.
  KHitWireX::~KHitWireX()
  {}

  bool KHitWireX::subpredict(const KETrack& tre,
			     KVector<1>::type& pvec,
			     KSymMatrix<1>::type& perr,
			     KHMatrix<1>::type& hmatrix) const
  {
    // Make sure that the track surface and the measurement surface are the same.
    // Throw an exception if they are not.

    if(!getMeasSurface()->isEqual(*tre.getSurface()))
      throw cet::exception("KHitWireX") << "Track surface not the same as measurement surface.\n";
 
    // Prediction is just u track perameter and error.

    int size = tre.getVector().size();
    pvec.resize(1, /* preserve */ false);
    pvec.clear();
    pvec(0) = tre.getVector()(0);

    perr.resize(1, /* preserve */ false);
    perr.clear();
    perr(0,0) = tre.getError()(0,0);

    // Hmatrix - du/du = 1., all others are zero.

    hmatrix.resize(1, size, /* preserve */ false);
    hmatrix.clear();
    hmatrix(0,0) = 1.;

    return true;
  }
} // end namespace trkf
