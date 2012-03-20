///////////////////////////////////////////////////////////////////////
///
/// \file   KFitTrack.cxx
///
/// \brief  Basic Kalman filter track class, with fit information.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "TrackFinder/KFitTrack.h"

namespace trkf {

  /// Default constructor.
  KFitTrack::KFitTrack() :
    fPath(0.),
    fChisq(0.),
    fStat(INVALID)
  {}

  /// Initializing constructor.
  ///
  /// Arguments:
  ///
  /// tre   - KETrack.
  /// s     - Path distance.
  /// chisq - Fit chisquare.
  /// stat  - Fit status.
  ///
  KFitTrack::KFitTrack(const KETrack& tre,
		       double s,
		       double chisq,
		       FitStatus stat) :
    KETrack(tre),
    fPath(s),
    fChisq(chisq),
    fStat(stat)
  {}

  /// Destructor.
  KFitTrack::~KFitTrack()
  {}

} // end namespace trkf
