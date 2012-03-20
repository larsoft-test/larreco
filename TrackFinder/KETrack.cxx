///////////////////////////////////////////////////////////////////////
///
/// \file   KETrack.cxx
///
/// \brief  Basic Kalman filter track class, with error.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "TrackFinder/KETrack.h"

namespace trkf {

  /// Default constructor.
  KETrack::KETrack()
  {
    fErr.clear();
  }

  /// Constructor - specify surface only.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  ///
  KETrack::KETrack(boost::shared_ptr<const Surface> psurf) :
    KTrack(psurf)
  {
    fErr.clear();
  }

  /// Constructor - surface + track parameters + error matrix.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  /// vec   - Track state vector.
  /// err   - Track error matrix.
  /// dir   - Track direction.
  ///
  KETrack::KETrack(boost::shared_ptr<const Surface> psurf,
		   const TrackVector& vec,
		   const TrackError& err, 
		   Surface::TrackDirection dir) :
    KTrack(psurf, vec, dir),
    fErr(err)
  {}

  /// Constructor - KTrack + error matrix.
  ///
  /// Arguments:
  ///
  /// trk - KTrack.
  /// err - Track error matrix.
  ///
  KETrack::KETrack(const KTrack& trk, const TrackError& err) :
    KTrack(trk),
    fErr(err)
  {}

  /// Destructor.
  KETrack::~KETrack()
  {}

} // end namespace trkf
