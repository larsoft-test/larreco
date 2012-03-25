///////////////////////////////////////////////////////////////////////
///
/// \file   KHitBase.cxx
///
/// \brief  Base class for Kalman filter measurement.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "TrackFinder/KHitBase.h"

namespace trkf {

  /// Default Constructor.
  KHitBase::KHitBase()
  {}

  /// Initializing Constructor.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  ///
  KHitBase::KHitBase(const boost::shared_ptr<const Surface>& psurf) :
    fSurf(psurf)
  {}

  /// Destructor.
  KHitBase::~KHitBase()
  {}

} // end namespace trkf
