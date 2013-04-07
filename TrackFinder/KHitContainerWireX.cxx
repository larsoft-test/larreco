///////////////////////////////////////////////////////////////////////
///
/// \file   KHitContainerWireX.cxx
///
/// \brief  A KHitContainer for KHitWireX type measurements.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include <map>
#include <stdint.h>

#include "TrackFinder/KHitContainerWireX.h"
#include "TrackFinder/KHitWireX.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Geometry/Geometry.h"

namespace trkf {

  /// Default Constructor.
  KHitContainerWireX::KHitContainerWireX()
  {}

  /// Destructor.
  KHitContainerWireX::~KHitContainerWireX()
  {}

  /// Fill container.
  ///
  /// Arguments:
  ///
  /// hits       - RecoBase/Hit collection.
  /// only_plane - Choose hits from this plane if >= 0.
  ///
  /// This method converts the hits in the input collection into
  /// KHitWireX objects and inserts them into the base class.  Hits
  /// corresponding to the same readout wire are grouped together as
  /// KHitGroup objects.
  ///
  void KHitContainerWireX::fill(const art::PtrVector<recob::Hit>& hits,
				int only_plane)
  {
    // Get services.

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;

    // Make a temporary map from channel number to KHitGroup objects.
    // The KHitGroup pointers are borrowed references to KHitGroup
    // objects stored by value in the base class.

    std::map<unsigned int, KHitGroup*> group_map;

    // Loop over hits.

    for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
	ihit != hits.end(); ++ihit) {
      const recob::Hit& hit = **ihit;

      // Extract the channel number from the Hit.
      geo::WireID hitWireID = hit.WireID();
		
      uint32_t channel = hit.Channel();

      // Extract the cryostat, tpc, and plane.
      // We only care about the plane.

	  //unsigned int cstat, tpc, plane, wire;
      //geom->ChannelToWire(channel, cstat, tpc, plane, wire);

		// Choose plane.

      if(only_plane >= 0 && hitWireID.Plane != (unsigned int)(only_plane))
	continue;

      // See if we need to make a new KHitGroup.

      KHitGroup* pgr = 0;
      if(group_map.count(channel) == 0) {
	getUnsorted().push_back(KHitGroup());
	pgr = &(getUnsorted().back());
	group_map[channel] = pgr;
      }
      else
	pgr = group_map[channel];
      assert(pgr != 0);

      // Get surface from KHitGroup (might be null pointer).

      const std::shared_ptr<const Surface>& psurf = pgr->getSurface();

      // Construct KHitWireX object.

      std::shared_ptr<const KHitBase> phit(new KHitWireX(*ihit, psurf));

      // Insert hit into KHitGroup.

      pgr->addHit(phit);
    }
  }

} // end namespace trkf
