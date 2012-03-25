////////////////////////////////////////////////////////////////////////
///
/// \file   KTrack.h
///
/// \brief  Basic Kalman filter track class, without error.
///
/// \author H. Greenlee 
///
/// This class include the following attributes.
///
/// 1. Surface.
/// 2. Track state vector.
/// 3. Track direction parameter.
///
/// The surface attribute is polymorphic, and is held via
/// boost::shared_ptr type of smart pointer, which handles memory
/// management using reference-counted shared ownership.
///
////////////////////////////////////////////////////////////////////////

#ifndef KTRACK_H
#define KTRACK_H

#include "TrackFinder/KalmanLinearAlgebra.h"
#include "TrackFinder/Surface.h"
#include "boost/shared_ptr.hpp"

namespace trkf {

  class KTrack
  {
  public:

    /// Default constructor.
    KTrack();

    /// Constructor - specify surface only.
    KTrack(const boost::shared_ptr<const Surface>& psurf);

    /// Constructor - surface + track parameters.
    KTrack(boost::shared_ptr<const Surface> psurf,
	   const TrackVector& vec,
	   Surface::TrackDirection dir = Surface::UNKNOWN);

    /// Destructor.
    virtual ~KTrack();

    // Accessors.

    const boost::shared_ptr<const Surface>& getSurface() const {return fSurf;} ///< Surface.
    const TrackVector& getVector() const {return fVec;}                 ///< Track state vector.
    Surface::TrackDirection getDirection() const;                       ///< Track direction.

    // Modifiers.

    /// Set surface.
    void setSurface(const boost::shared_ptr<const Surface>& psurf) {fSurf = psurf;}
    void setVector(const TrackVector& vec) {fVec = vec;}                ///< Set state vector.
    void setDirection(Surface::TrackDirection dir) {fDir = dir;}        ///< Set direction.

    /// Test if track is valid.
    bool isValid() const {return getDirection() != Surface::UNKNOWN;}

    /// Get position of track.
    void getPosition(double xyz[3]) const;

    /// Get momentum vector of track.
    void getMomentum(double mom[3]) const;

  private:

    // Attributes.

    boost::shared_ptr<const Surface>   fSurf;   ///< Track surface.
    TrackVector fVec;                           ///< Track state vector.
    Surface::TrackDirection fDir;               ///< Track direction.
  };
}

#endif
