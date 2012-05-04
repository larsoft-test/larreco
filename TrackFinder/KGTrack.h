////////////////////////////////////////////////////////////////////////
///
/// \file   KGTrack.h
///
/// \brief  A collection of KHitTracks.
///
/// \author H. Greenlee 
///
/// This class is a colleciton of KHitTracks.  In other words, it has
/// full information every measurement, and every fit on every
/// measurement surface.  This is the maximum amount of information
/// that it is possible to have.
///
/// KHitTrack collection is stored as a map, indexed by path distance.
/// This organization makes it easy to find the one or two nearest 
/// KHitTrack objects to any path distance.
///
/// Note that by combining information from forward and backward fit
/// tracks (Kalman smoothing), it is possible to obtain optimal fit
/// information (based on all past and future hits) for each measurement
/// surface, or any other surface, including interior surfaces.
///
////////////////////////////////////////////////////////////////////////

#ifndef KGTRACK_H
#define KGTRACK_H

#include <map>
#include "TrackFinder/KHitTrack.h"
#include "TrackFinder/SpacePointAlg.h"
#include "RecoBase/Track.h"

namespace trkf {

  class KGTrack
  {
  public:

    /// Default constructor.
    KGTrack();

    /// Destructor.
    ~KGTrack();

    // Accessors.

    /// KHitTrack collection, indexed by path distance.
    const std::multimap<double, KHitTrack>& getTrackMap() const {return fTrackMap;}

    /// Track at start point.
    const KHitTrack& startTrack() const;

    /// Track at end point.
    const KHitTrack& endTrack() const;

    /// Validity flag.
    bool isValid() const {return fTrackMap.size() > 0;}

    // Modifiers.

    /// Modifiable KHitTrack collection, indexed by path distance.
    std::multimap<double, KHitTrack>& getTrackMap() {return fTrackMap;}

    /// Modifiable track at start point.
    KHitTrack& startTrack();

    /// Modifiable track at end point.
    KHitTrack& endTrack();

    /// Add track.
    void addTrack(const KHitTrack& trh);

    /// Recalibrate track map.
    void recalibrate();

    /// Clear track collection.
    void clear() {fTrackMap.clear();}

    // Methods.

    /// Fill a recob::Track.
    void fillTrack(recob::Track& track) const;

    /// Fill a PtrVector of Hits.
    void fillHits(art::PtrVector<recob::Hit>& hits) const;

    /// Fill a collection of space points.
    void fillSpacePoints(std::vector<recob::SpacePoint>& spts,
			 const SpacePointAlg& sptalg) const;

  private:

    // Attributes.

    /// KHitTrack collection, indexed by path distance.
    std::multimap<double, KHitTrack> fTrackMap;
  };
}

#endif
