////////////////////////////////////////////////////////////////////////
///
/// \file   KHitWireX.h
///
/// \brief  Kalman filter wire-time measurement on a SurfWireX surface.
///
/// \author H. Greenlee 
///
/// This class is a type of one-dimensional Kalman filter measurement
/// reprsenting a single wire-time hit on a surface parallel to the
/// x-axis (nonmagnetic LAr tpc).
///
/// This class derives from base class KHit<1>, which is the general
/// one-dimensional measurement base class.  It does not have any
/// data attributes of its own.  It adds its own constructor and
/// overrides virtual base class method subpredict.
///
/// The following data are needed to fully specify an object of this
/// class.
///
/// 1.  Channel (defines measurement surface) or surface.
/// 2.  X position.
/// 3.  X error.
///
/// The x position and error are specified in the global coordinate
/// system, which is the same as the local u coordinate of the
/// measurement surface coordinate system.
///
////////////////////////////////////////////////////////////////////////

#ifndef KHITWIREX_H
#define KHITWIREX_H

#include "TrackFinder/KHit.h"

namespace trkf {

  class KHitWireX : public KHit<1>
  {
  public:

    /// Constructor from channel.
    KHitWireX(unsigned int channel, double x, double xerr);

    /// Constructor from surface.
    KHitWireX(const boost::shared_ptr<const Surface>& psurf,
	      double x, double xerr, int plane);

    /// Destructor.
    virtual ~KHitWireX();

    // Overrides.
    virtual bool subpredict(const KETrack& tre,
			    KVector<1>::type& pvec,
			    KSymMatrix<1>::type& perr,
			    KHMatrix<1>::type& hmatrix) const;
  };
}

#endif
