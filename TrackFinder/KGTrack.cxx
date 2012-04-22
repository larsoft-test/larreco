///////////////////////////////////////////////////////////////////////
///
/// \file   KGTrack.cxx
///
/// \brief  A collection of KHitTracks.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "TrackFinder/KGTrack.h"
#include "cetlib/exception.h"

namespace trkf {

  /// Default constructor.
  KGTrack::KGTrack()
  {}

  /// Destructor.
  KGTrack::~KGTrack()
  {}

  /// Track at start point.
  const KHitTrack& KGTrack::startTrack() const
  {
    /// Throw exception if track is not valid.

    if(!isValid())
      throw cet::exception("KGTrack") << "Starting track is not valid.\n";

    // Return track.

    return (*fTrackMap.begin()).second;
  }

  /// Track at end point.
  const KHitTrack& KGTrack::endTrack() const
  {
    /// Throw exception if track is not valid.

    if(!isValid())
      throw cet::exception("KGTrack") << "Ending track is not valid.\n";

    // Return track.

    return (*fTrackMap.rbegin()).second;
  }

  /// Modifiable track at start point.
  KHitTrack& KGTrack::startTrack()
  {
    /// Throw exception if track is not valid.

    if(!isValid())
      throw cet::exception("KGTrack") << "Starting track is not valid.\n";

    // Return track.

    return (*fTrackMap.begin()).second;
  }

  /// Modifiable track at end point.
  KHitTrack& KGTrack::endTrack()
  {
    /// Throw exception if track is not valid.

    if(!isValid())
      throw cet::exception("KGTrack") << "Ending track is not valid.\n";

    // Return track.

    return (*fTrackMap.rbegin()).second;
  }

  /// Add track.
  void KGTrack::addTrack(const KHitTrack& trh) {
    fTrackMap.insert(std::make_pair(trh.getPath(), trh));
  }

} // end namespace trkf
