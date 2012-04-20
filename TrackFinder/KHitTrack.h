////////////////////////////////////////////////////////////////////////
///
/// \file   KHitTrack.h
///
/// \brief  Basic Kalman filter track class, plus one measurement on same surface.
///
/// \author H. Greenlee 
///
/// This class inherits the following attributes from KFitTrack.
///
/// 1. Surface.
/// 2. Track state vector.
/// 3. Track direction parameter.
/// 4. Track error matrix.
/// 5. Propagation distance.
/// 6. Fit chisquare.
/// 7. Fit status.
///
/// This class adds the following attributes of its own.
///
/// 8. A single measurement.
///
////////////////////////////////////////////////////////////////////////

#ifndef KHITSTRACK_H
#define KHITSTRACK_H

#include <vector>
#include "TrackFinder/KFitTrack.h"
#include "TrackFinder/KHitBase.h"

namespace trkf {

  class KHitTrack : public KFitTrack
  {
  public:

    /// Default constructor.
    KHitTrack();

    /// Initializing constructor - KFitTrack + measurement.
    KHitTrack(const KFitTrack& trf, const boost::shared_ptr<const KHitBase>& hit);

    /// Initializing constructor - KETrack.
    KHitTrack(const KETrack& tre);

    /// Destructor.
    virtual ~KHitTrack();

    // Accessor.

    /// Measurement.
    const boost::shared_ptr<const KHitBase>& getHit() {return fHit;}

    // Modifiers.

    /// Set measurement.
    void setHit(const boost::shared_ptr<const KHitBase>& hit) {fHit = hit;}

    /// Printout
    virtual std::ostream& Print(std::ostream& out, bool doTitle = true) const;

  private:

    // Attributes.

    /// Measurement.
    boost::shared_ptr<const KHitBase> fHit;
  };
}

#endif
