///////////////////////////////////////////////////////////////////////
///
/// \file   Surface.cxx
///
/// \brief  Base class for Kalman filter surface.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "TrackFinder/Surface.h"

namespace trkf {

  /// Default constructor.
  Surface::Surface()
  {}

  /// Destructor.
  Surface::~Surface()
  {}

  /// Output operator.
  std::ostream& operator<<(std::ostream& out, const Surface& surf)
  {
    return surf.Print(out);
  }

} // end namespace trkf
