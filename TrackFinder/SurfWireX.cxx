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
#include "Geometry/Geometry.h"
#include "Geometry/WireGeo.h"
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

    //unsigned int cstat, tpc, plane, wire;
	
	std::vector<geo::WireID> vecChannelWireID = geom->ChannelToWire(channel);
	  
	if (vecChannelWireID.size() != 0) {
	  
	  geo::WireGeo const& wgeom = geom->WireIDToWireGeo(vecChannelWireID[0]);
	  
	  // Get wire center and angle from the wire geometry.
	  // Put local origin at center of wire.
	  
	  double xyz[3] = {0.};
	  wgeom.GetCenter(xyz);
	  double phi = TMath::PiOver2() - wgeom.ThetaZ();
	  
	  // Update base class.
	  
	  *static_cast<SurfYZPlane*>(this) = SurfYZPlane(xyz[1], xyz[2], phi);
	} else {
		throw cet::exception("Can't Find Wire for that channel") << channel; 

	}

  }

  /// Destructor.
  SurfWireX::~SurfWireX()
  {}

} // end namespace trkf
