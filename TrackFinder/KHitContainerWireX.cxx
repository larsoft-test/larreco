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
#include "TrackFinder/KHitContainerWireX.h"
#include "TrackFinder/KHitWireX.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Geometry/geo.h"

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

      // Extract the information that we need from the Hit.
      // That's channel, time, and time error.

      unsigned int channel = hit.Channel();
      double t = hit.PeakTime();
      double terr = hit.SigmaPeakTime();

      // Don't let the time error be less than 1./sqrt(12.).
      // This should be removed when hit errors are fixed.

      if(terr < 1./std::sqrt(12.))
	terr = 1./std::sqrt(12.);

      // Extract the cryostat, tpc, and plane.
      // We need these to get the time offset.

      unsigned int cstat, tpc, plane, wire;
      geom->ChannelToWire(channel, cstat, tpc, plane, wire);

      // Choose plane.

      if(only_plane >= 0 && plane != (unsigned int)(only_plane))
	continue;

      // Calculate position and error.

      double x = detprop->ConvertTicksToX(t, plane, tpc, cstat);
      double xerr = terr * detprop->GetXTicksCoefficient();

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

      const boost::shared_ptr<const Surface>& psurf = pgr->getSurface();

      // Construct KHitWireX object.
      // If the surface pointer is valid, use that surface.
      // Otherwise make a new surface.

      boost::shared_ptr<const KHitBase> phit;
      if(psurf.get() == 0)
	phit = boost::shared_ptr<const KHitBase>(new KHitWireX(channel, x, xerr));
      else
	phit = boost::shared_ptr<const KHitBase>(new KHitWireX(psurf, x, xerr, pgr->getPlane()));

      // Insert hit into KHitGroup.

      pgr->addHit(phit);
    }
  }

} // end namespace trkf
