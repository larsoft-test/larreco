///////////////////////////////////////////////////////////////////////
///
/// \file   KGTrack.cxx
///
/// \brief  A collection of KHitTracks.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "TrackFinder/KGTrack.h"
#include "TrackFinder/KHitWireX.h"
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

  /// Fill a recob::Track.
  ///
  /// Arguments:
  ///
  /// track - Track to fill.
  ///
  void KGTrack::fillTrack(recob::Track& track) const
  {
    // Fill collections of trajectory points and direction vectors.

    std::vector<TVector3> xyz;
    std::vector<TVector3> dxdydz;
    std::vector<double> momentum;
    std::vector<std::vector<double> > dQdx;

    xyz.reserve(fTrackMap.size());
    dxdydz.reserve(fTrackMap.size());
    momentum.reserve(fTrackMap.size());

    // Loop over KHitTracks.

    for(std::multimap<double, KHitTrack>::const_iterator itr = fTrackMap.begin();
	itr != fTrackMap.end(); ++itr) {
      const KHitTrack& trh = (*itr).second;

      // Get position.

      double pos[3];
      trh.getPosition(pos);
      xyz.push_back(TVector3(pos[0], pos[1], pos[2]));

      // Get momentum vector.

      double mom[3];
      trh.getMomentum(mom);
      double p = std::sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
      assert(p != 0.);
      dxdydz.push_back(TVector3(mom[0]/p, mom[1]/p, mom[2]/p));
      momentum.push_back(p);
    }

    // Fill track.

    track = recob::Track(xyz, dxdydz, dQdx, momentum);
  }

  /// Fill a PtrVector of Hits.
  ///
  /// Arguments:
  ///
  /// hits - Hit vector to fill.
  ///
  void KGTrack::fillHits(art::PtrVector<recob::Hit>& hits) const
  {
    hits.reserve(hits.size() + fTrackMap.size());

    // Loop over KHitTracks and fill hits belonging to this track.

    for(std::multimap<double, KHitTrack>::const_iterator it = fTrackMap.begin();
	it != fTrackMap.end(); ++it) {
      const KHitTrack& track = (*it).second;

      // Extrack Hit from track.

      const boost::shared_ptr<const KHitBase>& hit = track.getHit();
      const KHitWireX* phit = dynamic_cast<const KHitWireX*>(&*hit);
      if(phit != 0) {
	const art::Ptr<recob::Hit> prhit = phit->getHit();
	if(!prhit.isNull())
	  hits.push_back(prhit);
      }
    }
  }

  /// Fill a collection of space points.
  ///
  /// Arguments:
  ///
  /// spts   - Collection of space points to fill.
  /// sptalg - Space point algorithm object.
  ///
  /// This method uses the hits contained in this track to construct
  /// space points.
  ///
  /// This method does not have any knowledge of what constitutes a
  /// good space point, except that Hits are required to be
  /// consecutive when sorted by path distance, and space points are
  /// required to pass compatibility tests used by the space point
  /// algorithm object.  This method will make space points from
  /// either two or three Hits (even for three-plane detectors), if
  /// the space point algorithm is configured to allow it.
  ///
  void KGTrack::fillSpacePoints(std::vector<recob::SpacePoint>& spts,
				const SpacePointAlg& sptalg) const
  {
    // Loop over KHitTracks.

    art::PtrVector<recob::Hit> hits;
    art::PtrVector<recob::Hit> compatible_hits;
    for(std::multimap<double, KHitTrack>::const_iterator it = fTrackMap.begin();
	it != fTrackMap.end(); ++it) {
      const KHitTrack& track = (*it).second;

      // Extrack Hit from track.

      const boost::shared_ptr<const KHitBase>& hit = track.getHit();
      const KHitWireX* phit = dynamic_cast<const KHitWireX*>(&*hit);
      if(phit != 0) {
	const art::Ptr<recob::Hit> prhit = phit->getHit();

	// Test this hit for compatibility.

	hits.push_back(prhit);
	bool ok = sptalg.compatible(hits);
	if(!ok) {

	  // The new hit is not compatible.  Make a space point out of
	  // the last known compatible hits, provided there are at least
	  // two.

	  if(compatible_hits.size() >= 2) {
	    spts.push_back(recob::SpacePoint());
	    spts.back().SetID(sptalg.numHitMap());
	    sptalg.fillSpacePoint(compatible_hits, spts.back());
	    compatible_hits.clear();
	  }

	  // Forget about any previous hits.

	  hits.clear();
	  hits.push_back(prhit);
	}

	// Update the list of known compatible hits.

	compatible_hits = hits;
      }
    }

    // Maybe make one final space point.

    if(compatible_hits.size() >= 2) {
      spts.push_back(recob::SpacePoint());
      spts.back().SetID(sptalg.numHitMap());
      sptalg.fillSpacePoint(compatible_hits, spts.back());
    }
  }
} // end namespace trkf
