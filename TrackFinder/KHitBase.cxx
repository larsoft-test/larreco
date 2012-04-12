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

  /// Printout
  std::ostream& KHitBase::Print(std::ostream& out, bool doTitle) const
  {
    if(doTitle)
      out << "KHitBase:\n";
    out << "  Surface: " << *fSurf << "\n";
    return out;
  }

  /// Output operator.
  std::ostream& operator<<(std::ostream& out, const KHitBase& trk)
  {
    return trk.Print(out);
  }

} // end namespace trkf
