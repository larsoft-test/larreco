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
/// 2.  Measurement plane index.
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
/// This class has two surface data members.  Member fMeasSurf is the
/// idealized measurement surface, which is set on construction and
/// never changes.  Member fPredSurf is used to remember the track
/// surface used to make a prediction.  It is updated (by derived
/// class) every time method predict is called.
///
/// As with KTrack, the surface attributes are polymorphic, and is held
/// via boost::shared_ptr type of smart pointer, which handles memory
/// management using reference-counted shared ownership.
///
////////////////////////////////////////////////////////////////////////

#ifndef KHITBASE_H
#define KHITBASE_H

#include <ostream>
#include "TrackFinder/Surface.h"
#include "TrackFinder/KETrack.h"
#include "TrackFinder/Propagator.h"
#include "boost/shared_ptr.hpp"

namespace trkf {

  class KHitBase
  {
  public:

    /// Default constructor.
    KHitBase();

    /// Initializing Constructor.
    KHitBase(const boost::shared_ptr<const Surface>& psurf, int plane = -1);

    /// Destructor.
    virtual ~KHitBase();

    // Accessors.

    /// Predition surface.
    const boost::shared_ptr<const Surface>& getPredSurface() const {return fPredSurf;}

    /// Measurement surface.
    const boost::shared_ptr<const Surface>& getMeasSurface() const {return fMeasSurf;}

    /// Measurement plane index.
    int getMeasPlane() const {return fMeasPlane;}

    // Modifiers.
    void setMeasPlane(int plane) {fMeasPlane = plane;}

    // Pure virtual methods.

    /// Prediction method (return false if fail).
    virtual bool predict(const KETrack& tre, const Propagator* prop = 0) const = 0;

    /// Return incremental chisquare.
    virtual double getChisq() const = 0;

    /// Update track method.
    virtual void update(KETrack& tre) const = 0;

    /// Printout
    virtual std::ostream& Print(std::ostream& out, bool doTitle = true) const;

    // Attributes.

  protected:

    mutable boost::shared_ptr<const Surface> fPredSurf;   ///< Prediction surface.

  private:

    boost::shared_ptr<const Surface> fMeasSurf;           ///< Measurement surface.
    int fMeasPlane;                                       ///< Measurement plane index.
  };

  /// Output operator.
  std::ostream& operator<<(std::ostream& out, const KHitBase& trk);
}

#endif
