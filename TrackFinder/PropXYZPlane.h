////////////////////////////////////////////////////////////////////////
///
/// \file   PropXYZPlane.h
///
/// \brief  Propagate between two SurfXYZPlanes.
///
/// \author H. Greenlee 
///
/// Class for propagating between two SurfXYZPlane or SurfYZPlane 
/// surfaces.  If the initial or destination surface is not a
/// SurfXYZPlane or SurfYZPlane, return propation failure.
///
////////////////////////////////////////////////////////////////////////

#ifndef PROPXYZPLANE_H
#define PROPXYZPLANE_H

#include "TrackFinder/Propagator.h"

namespace trkf {

  class PropXYZPlane : public trkf::Propagator
  {
  public:

    /// Constructor.
    PropXYZPlane(double tcut);

    /// Destructor.
    virtual ~PropXYZPlane();

    // Overrides.

    /// Clone method.
    Propagator* clone() const {return new PropXYZPlane(*this);}

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
