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
    fTrackMap.insert(std::make_pair(trh.getPath() + trh.getHit()->getPredDistance(), trh));
  }

  /// Recalibrate track map.
  ///
  /// Loop over contents of track map.  Copy each KHitTrack into a new multimap track map.
  /// Offset the distance stored in the KHitTracks such that the minimum distance is zero.
  /// Also update multimap keys to agree with distance stored in track.
  ///
  void KGTrack::recalibrate()
  {
    std::multimap<double, KHitTrack> newmap;

    // Loop over old track map.

    bool first = true;
    double s0 = 0.;
    for(std::multimap<double, KHitTrack>::iterator i = fTrackMap.begin();
	i != fTrackMap.end(); ++i) {
      KHitTrack& trh = (*i).second;
      if(first) {
	first = false;
	s0 = trh.getPath();
      }
      double s = trh.getPath()  - s0;
      trh.setPath(s);
      newmap.insert(std::make_pair(s, trh));
    }

    // Update data member track map.

    fTrackMap.swap(newmap);
  }

} // end namespace trkf
