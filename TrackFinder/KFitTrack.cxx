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

  /// Printout
  std::ostream& KFitTrack::Print(std::ostream& out, bool doTitle) const
  {
    if(doTitle)
      out << "KFitTrack:\n";

    // Print information specific to this class.

    out << "  Distance = " << fPath << "\n";
    out << "  Chisquare = " << fChisq << "\n";
    out << "  Status = " << (fStat == INVALID ? "INVALID" :
			     (fStat == UNKNOWN ? "UNKNOWN" :
			      (fStat == FORWARD ? "FORWARD" :
			       (fStat == FORWARD_PREDICTED ? "FORWARD_PREDICTED" :
				(fStat == BACKWARD ? "BACKWARD" :
				 (fStat == BACKWARD_PREDICTED ? "BACKWARD_PREDICTED" :
				  (fStat == OPTIMAL ? "OPTIMAL" : 
				   "OPTIMAL_PREDICTED"))))))) << "\n";

    // Print base class.

    KETrack::Print(out, false);
    return out;
  }

} // end namespace trkf
