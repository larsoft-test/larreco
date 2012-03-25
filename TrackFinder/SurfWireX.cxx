///////////////////////////////////////////////////////////////////////
///
/// \file   SurfWireX.cxx
///
/// \brief  Planar surface defined by readout wire and x-axis.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "TrackFinder/SurfWireX.h"
#include "Geometry/geo.h"
#include "TMath.h"

namespace trkf {

  /// Constructor.
  ///
  /// Arguments:
  ///
  /// channel - Readout channel number.
  ///
  SurfWireX::SurfWireX(unsigned int channel)
  {
    // Get geometry service.

    art::ServiceHandle<geo::Geometry> geom;

    // Get wire geometry.

    unsigned int cstat, tpc, plane, wire;
    const geo::WireGeo& wgeom = geom->ChannelToWire(channel, cstat, tpc, plane, wire);

    // Get wire center and angle from the wire geometry.
    // Put local origin at center of wire.

    double xyz[3];
    wgeom.GetCenter(xyz);
    double phi = TMath::PiOver2() - wgeom.ThetaZ();

    // Update base class.

    *static_cast<SurfYZPlane*>(this) = SurfYZPlane(xyz[1], xyz[2], phi);
  }

  /// Destructor.
  SurfWireX::~SurfWireX()
  {}

} // end namespace trkf
