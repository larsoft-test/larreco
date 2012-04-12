////////////////////////////////////////////////////////////////////////
///
/// \file   KHitBase.h
///
/// \brief  Base class for Kalman filter measurement.
///
/// \author H. Greenlee 
///
/// This class represents a general measurement on a surface.
///
/// This class includes the following attribute.
///
/// 1.  Measurement surface.
///
/// This class has the following pure virtual methods.
///
/// 1.  Prediction method, in which the predicted track object (state
///     vector + error matrix) is passed into the measurement to
///     generate a prediction in the measurement coordinate system.
///
/// 2.  Accessor for incremental chisquare.
///
/// 3.  Update method, in which the track object is passed in and is
///     updated according to the Kalman updating formula.
///
/// This class does not include in its interface anything having to do
/// with concrete measurements, predictions, or residuals, or anything
/// with a variable dimension.
///
/// As with KTrack, the surface attribute is polymorphic, and is held
/// via boost::shared_ptr type of smart pointer, which handles memory
/// management using reference-counted shared ownership.
///
////////////////////////////////////////////////////////////////////////

#ifndef KHITBASE_H
#define KHITBASE_H

#include <ostream>
#include "TrackFinder/Surface.h"
#include "TrackFinder/KETrack.h"
#include "boost/shared_ptr.hpp"

namespace trkf {

  class KHitBase
  {
  public:

    /// Default constructor.
    KHitBase();

    /// Initializing Constructor.
    KHitBase(const boost::shared_ptr<const Surface>& psurf);

    /// Destructor.
    virtual ~KHitBase();

    // Accessors.

    const boost::shared_ptr<const Surface>& getSurface() const {return fSurf;} ///< Surface.

    // Pure virtual methods.

    /// Prediction method (return false if fail).
    virtual bool predict(const KETrack& tre) const = 0;

    /// Return incremental chisquare.
    virtual double getChisq() const = 0;

    /// Update track method.
    virtual void update(KETrack& tre) const = 0;

    /// Printout
    virtual std::ostream& Print(std::ostream& out, bool doTitle = true) const;

  private:

    // Attributes.

    boost::shared_ptr<const Surface>   fSurf;   ///< Track surface.
  };

  /// Output operator.
  std::ostream& operator<<(std::ostream& out, const KHitBase& trk);
}

#endif
