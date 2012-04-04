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
    fDir(Surface::UNKNOWN),
    fPdgCode(0)
  {}

  /// Constructor - specify surface only.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  ///
  KTrack::KTrack(const boost::shared_ptr<const Surface>& psurf) :
    fSurf(psurf),
    fDir(Surface::UNKNOWN),
    fPdgCode(0)
  {}

  /// Constructor - surface + track parameters.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  /// vec   - Track state vector.
  /// dir   - Track direction.
  /// pdg   - Pdg code.
  ///
  KTrack::KTrack(boost::shared_ptr<const Surface> psurf,
		 const TrackVector& vec,
		 Surface::TrackDirection dir,
		 int pdg) :
    fSurf(psurf),
    fVec(vec),
    fDir(dir),
    fPdgCode(pdg)
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

  /// Printout
  std::ostream& KTrack::Print(std::ostream& out, bool doTitle) const
  {
    if(doTitle)
      out << "KTrack:\n";
    out << "  Surface direction = " << (fDir == Surface::FORWARD ?
					"FORWARD" : 
					( fDir == Surface::BACKWARD ?
					  "BACKWARD" : "UNKNOWN" )) << "\n"
	<< "  Pdg = " << fPdgCode << "\n"
	<< "  Surface: " << *fSurf << "\n"
	<< "  Track parameters:\n"
	<< "  [";
    for(unsigned int i = 0; i < fVec.size(); ++i) {
      if(i != 0)
	out << ", ";
      out << fVec(i);
    }
    out << "]\n";
    return out;
  }

  /// Output operator.
  std::ostream& operator<<(std::ostream& out, const KTrack& trk)
  {
    return trk.Print(out);
  }

} // end namespace trkf
