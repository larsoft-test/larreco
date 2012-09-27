////////////////////////////////////////////////////////////////////////
///
/// \file   PropYZPlane.h
///
/// \brief  Propagate between two SurfYZPlanes.
///
/// \author H. Greenlee 
///
/// Class for propagating between two SurfYZPlane surfaces.  If the
/// initial or destination surface is not a SurfYZPlane, return
/// propation failure.
///
////////////////////////////////////////////////////////////////////////

#ifndef PROPYZPLANE_H
#define PROPYZPLANE_H

#include "TrackFinder/Propagator.h"

namespace trkf {

  class PropYZPlane : public trkf::Propagator
  {
  public:

    /// Constructor.
    PropYZPlane(double tcut);

    /// Destructor.
    virtual ~PropYZPlane();

    // Overrides.

    /// Clone method.
    Propagator* clone() const {return new PropYZPlane(*this);}

    /// Propagate without error.
    boost::optional<double> short_vec_prop(KTrack& trk,
					   const std::shared_ptr<const Surface>& surf, 
					   Propagator::PropDirection dir, 
					   bool doDedx,
					   TrackMatrix* prop_matrix = 0,
					   TrackError* noise_matrix = 0) const;
  };
}

#endif
