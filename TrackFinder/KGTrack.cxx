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
#include "TrackFinder/SurfXYZPlane.h"
#include "TrackFinder/PropXYZPlane.h"
#include "Geometry/Geometry.h"
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
  void KGTrack::fillTrack(recob::Track& track, 
			  int id) const
  {
    // Get geometry service.

    art::ServiceHandle<geo::Geometry> geom;
    int nview = geom->Nviews();

    // Make propagator for propating to standard track surface.

    PropXYZPlane prop(0.);

    // Fill collections of trajectory points and direction vectors.

    std::vector<TVector3> xyz;
    std::vector<TVector3> dxdydz;
    std::vector<TMatrixD> cov;
    std::vector<double> momentum;
    std::vector<std::vector<double> > dqdx(nview);

    xyz.reserve(fTrackMap.size());
    dxdydz.reserve(fTrackMap.size());
    momentum.reserve(fTrackMap.size());
    for(int iview=0; iview<nview; ++iview)
      dqdx[iview].reserve(fTrackMap.size());

    // Loop over KHitTracks.

    unsigned int n = 0;
    for(std::multimap<double, KHitTrack>::const_iterator itr = fTrackMap.begin();
	itr != fTrackMap.end(); ++itr, ++n) {
      const KHitTrack& trh = (*itr).second;

      // Get position.

      double pos[3];
      trh.getPosition(pos);
      xyz.push_back(TVector3(pos[0], pos[1], pos[2]));

      // Get momentum vector.
      // Fill direction unit vector and momentum.

      double mom[3];
      trh.getMomentum(mom);
      double p = std::sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
      assert(p != 0.);
      dxdydz.push_back(TVector3(mom[0]/p, mom[1]/p, mom[2]/p));
      momentum.push_back(p);

      // Fill error matrix.

      if(n == 0 || n == fTrackMap.size() - 1) {

	TMatrixD covar(5,5);

	// Construct surface perpendicular to track momentun, and
	// propagate track to that surface (zero distance).

	const std::shared_ptr<const Surface> psurf(new SurfXYZPlane(pos[0], pos[1], pos[2],
								    mom[0], mom[1], mom[2]));
	KETrack tre(trh);
	boost::optional<double> dist = prop.err_prop(tre, psurf, Propagator::UNKNOWN, false);
	assert(!!dist);
	for(int i=0; i<5; ++i) {
	  for(int j=0; j<5; ++j)
	    covar(i,j) = tre.getError()(i,j);
	}
	cov.push_back(covar);
      }

      // Get charge.
      // Only implemented for KHitWireX type measurements.

      for(int iview=0; iview<nview; ++iview)
	dqdx[iview].push_back(0.);
      const std::shared_ptr<const KHitBase>& phit = trh.getHit();
      if(phit.get() != 0) {
	const KHitWireX* phitx = dynamic_cast<const KHitWireX*>(&*phit);
	if(phitx != 0) {
	  const art::Ptr<recob::Hit>& parthit = phitx->getHit();
	  if(parthit.get() != 0) {
	    const recob::Hit& arthit = *parthit;
	    geo::View_t view = arthit.View();
	    double pitch = geom->WirePitch(view);
	    double charge = arthit.Charge();
	    double dudw = trh.getVector()[2];
	    double dvdw = trh.getVector()[3];
	    double dist = pitch * std::sqrt(1. + dudw * dudw + dvdw * dvdw);
	    double qdist = charge / dist;
	    assert(view >= 0 && view < nview);
	    dqdx[view].back() = qdist;
	  }
	}
      }
    }

    // Fill track.

    track = recob::Track(xyz, dxdydz, cov, dqdx, momentum, id);
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

      const std::shared_ptr<const KHitBase>& hit = track.getHit();
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

      const std::shared_ptr<const KHitBase>& hit = track.getHit();
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
	    sptalg.fillSpacePoint(compatible_hits, spts, sptalg.numHitMap());
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
      sptalg.fillSpacePoint(compatible_hits, spts, sptalg.numHitMap());
    }
  }
} // end namespace trkf
