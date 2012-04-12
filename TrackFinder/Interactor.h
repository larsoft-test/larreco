////////////////////////////////////////////////////////////////////////
///
/// \file   Interactor.h
///
/// \brief  Base class for Kalman filter track interactor.
///
/// \author H. Greenlee 
///
/// This class defined the general interface for calculating
/// propagation noise.
///
/// This class defined a single pure virtual method called "noise"
/// whose purpose is to calculate the propagation noise matrix
/// associated with the propagation of a track over a specified
/// distance.
///
////////////////////////////////////////////////////////////////////////

#ifndef INTERACTOR_H
#define INTERACTOR_H

#include "TrackFinder/KalmanLinearAlgebra.h"
#include "TrackFinder/KTrack.h"

namespace trkf {

  class Interactor
  {
  public:

    /// Constructor.
    Interactor(double tcut);

    /// Destructor.
    virtual ~Interactor();

    // Accessors.

    double getTcut() const {return fTcut;}

    // Virtual methods.

    /// Clone method.
    virtual Interactor* clone() const = 0;

    /// Calculate noise matrix.
    virtual void noise(const KTrack& trk, double s, TrackError& noise_matrix) const = 0;

  private:

    // Attributes.

    double fTcut;  ///< Maximum delta ray energy for dE/dx.
  };
}

#endif
