///////////////////////////////////////////////////////////////////////
///
/// \file   Interactor.cxx
///
/// \brief  Base class for Kalman filter track interactor.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "TrackFinder/Interactor.h"

namespace trkf {

  /// Constructor.
  ///
  /// Arguments:
  ///
  /// tcut - Maximum delta ray energy.
  ///
  Interactor::Interactor(double tcut) :
    fTcut(tcut)
  {}

  /// Destructor.
  Interactor::~Interactor()
  {}

} // end namespace trkf
