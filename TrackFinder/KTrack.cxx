///////////////////////////////////////////////////////////////////////
///
/// \file   KTrack.cxx
///
/// \brief  Basic Kalman filter track class, without error.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "TrackFinder/KTrack.h"
#include "cetlib/exception.h"

namespace trkf {

  /// Default constructor.
  KTrack::KTrack() :
    fVec(0.),
    fDir(Surface::UNKNOWN)
  {}

  /// Constructor - specify surface only.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  ///
  KTrack::KTrack(boost::shared_ptr<const Surface> psurf) :
    fSurf(psurf),
    fVec(0.),
    fDir(Surface::UNKNOWN)
  {}

  /// Constructor - surface + track parameters.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  /// vec   - Track state vector.
  /// dir   - Track direction.
  ///
  KTrack::KTrack(boost::shared_ptr<const Surface> psurf,
		 const TrackVector& vec,
		 Surface::TrackDirection dir) :
    fSurf(psurf),
    fVec(vec),
    fDir(dir)
  {}

  /// Destructor.
  KTrack::~KTrack()
  {}

  /// Track direction accessor.
  /// Track direction implied by track parameters has precedence
  /// over track direction attribute.
  /// If the surface pointer is null, return UNKNOWN.
  Surface::TrackDirection KTrack::getDirection() const
  {
    Surface::TrackDirection result = Surface::UNKNOWN;
    if(fSurf.get() != 0)
      result = fSurf->getDirection(fVec, fDir);
    return result;
  }

  /// Get position of track.
  /// Throw an exception if track is not valid.
  ///
  /// Arguments:
  ///
  /// xyz - Position vector.
  ///
  void KTrack::getPosition(double xyz[3]) const
  {
    if(!isValid())
      throw cet::exception("SurfYZPlane") << "Position requested for invalid track.\n";
    fSurf->getPosition(fVec, xyz);
  }

  /// Get momentum vector of track.
  /// Throw an exception if track is not valid.
  ///
  /// Arguments:
  ///
  /// mom - Momentum vector of track.
  ///
  void KTrack::getMomentum(double mom[3]) const
  {
    if(!isValid())
      throw cet::exception("SurfYZPlane") << "Momentum vector requested for invalid track.\n";
    Surface::TrackDirection dir = getDirection();
    fSurf->getMomentum(fVec, mom, dir);
  }

} // end namespace trkf
