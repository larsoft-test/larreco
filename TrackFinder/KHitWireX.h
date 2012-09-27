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
/// x-axis (approprite for nonmagnetic LAr tpc).
///
/// This class derives from base class KHit<1>, which is the general
/// one-dimensional measurement base class.  It is constructed from a
/// art::Ptr<recob::Hit>, which pointer is saved in this class as a
/// data member.  This class overrides virtual base class method
/// subpredict.
///
/// The following data are extracted from the Hit, and are stored in
/// the base class.
///
/// 1.  Channel (defines measurement surface).
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
#include "RecoBase/Hit.h"
#include "art/Persistency/Common/Ptr.h"

namespace trkf {

  class KHitWireX : public KHit<1>
  {
  public:

    /// Constructor from Hit.
    KHitWireX(const art::Ptr<recob::Hit>& hit,
	      const std::shared_ptr<const Surface>& psurf);

    /// Constructor from channel (mainly for testing).
    KHitWireX(unsigned int channel, double x, double xerr);

    /// Destructor.
    virtual ~KHitWireX();

    // Accessors.

    /// Get original hit.
    const art::Ptr<recob::Hit>& getHit() const {return fHit;}

    // Overrides.

    // Prediction method.
    virtual bool subpredict(const KETrack& tre,
			    KVector<1>::type& pvec,
			    KSymMatrix<1>::type& perr,
			    KHMatrix<1>::type& hmatrix) const;

  private:

    // Attributes.

    art::Ptr<recob::Hit> fHit;
  };
}

#endif
