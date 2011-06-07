////////////////////////////////////////////////////////////////////////
///
/// \file   SpacePointService.h
///
/// \brief  Service for generating space points from hits.
///
/// \author H. Greenlee 
///
/// This service calculates space points (recob::SpacePoint) from an 
/// unsorted collection of hits (recob::Hit).  The resulting space
/// points will contain one hit from from two or three views.
///
/// FCL parameters:
///
/// Debug - Enable debugging messages (0=none, 1=some, 2=more).
/// Hist - Enable histograms.
/// MCHist - Fill histograms involving MC information.
/// UseMC - Use MC truth information to find only true space points.
/// MaxDT - The maximum time difference (ticks) between any pair of hits.
/// MaxS  - The maximum 3-view wire separation parameter S (cm).
/// TimeOffsetU - Plane U time offset (ticks).
/// TimeOffsetV - Plane V time offset (ticks).
/// TimeOffsetW - Plane W time offset (ticks).
/// MinViews - Minimum number of views to make a space point (2 or 3).
/// EnableU - Use U view hits.
/// EnableV - Use V view hits.
/// EnableW - Use W view hits.
///
/// The first two parameters are used to implement a notion of whether
/// the input hits are compatible with being a space point.  Parameter
/// MaxS is a cut on the 3-plane wire separation parameter S, which is
/// defined as follows:
///
/// S = sin(theta_vw)*u + sin(theta_wu)*v + sin(theta_uv)*w
///
/// where wire coordinates (u,v,w) are measured in cm with respect to a 
/// common origin.
///
/// The time offsets are subtracted from times embedded in hits before
/// comparing times in different planes, and before converting time
/// to distance. 
///
/// There should eventually be a better way to specify time offsets.
///
////////////////////////////////////////////////////////////////////////

#ifndef SPACEPOINTSERVICE_H
#define SPACEPOINTSERVICE_H

#include <vector>
#include <string>
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "Geometry/geo.h"
#include "RecoBase/SpacePoint.h"

class TH1F;
namespace art{class Event;}
namespace sim{
  class SimChannel;
  class Electrons;
}

namespace trkf {

  class SpacePointService {
  public:

    // Constructor.
    SpacePointService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);

    // Destructor.
    ~SpacePointService();

    // Update configuration parameters.
    void reconfigure(const fhicl::ParameterSet& pset);

    // Update constants obtained from geometry or LAr properties services.
    void postBeginRun(const art::Run& run);   ///< Same as maybe_update.
    void maybe_update();    ///< Update constants if geometry has changed..
    void update();          ///< Update constants unconditionally.

    // Book histograms.
    void bookHistograms();

    // Test whether the specified hits are compatible with a space point.
    bool compatible(const art::PtrVector<recob::Hit>& hits) const;

    // Fill a single space point using the specified hits.
    // Hits are assumed to be compatible.
    void fillSpacePoint(const art::PtrVector<recob::Hit>& hits,
			recob::SpacePoint& spt) const;

    // Fill a vector of space points from an unsorted collection of hits.
    // Space points are generted for all compatible combinations of hits.
    // The event pointer is only needed if fcl parameter MCHist = true or
    // UseMC = true.
    void makeSpacePoints(const art::Handle< std::vector<recob::Hit> >& hith,
			 std::vector<recob::SpacePoint>& spts,
			 const art::Event* pevt = 0) const;
    void makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
			 std::vector<recob::SpacePoint>& spts,
			 const art::Event* pevt = 0) const;

  private:

    // Extract Electrons associated with Hit.
    // (Similar as BackTracker::HitToXYZ in MCCheater.)
    void HitToElectrons(const recob::Hit& hit,
			std::vector<const sim::Electrons*>& electrons) const;

    // Extract average position of vector of electrons.
    // (Similar as BackTracker::HitToXYZ in MCCheater.)
    TVector3 ElectronsToXYZ(const std::vector<const sim::Electrons*>& electrons) const;

    // Configuration paremeters.

    int fDebug;             ///< Debugging level.
    bool fHist;             ///< Enable histograms.
    bool fMCHist;           ///< Enable histograms involving MC information.
    bool fUseMC;            ///< Use MC truth information.
    std::string fMClabel;   ///< MC label for SimChannel.
    double fMaxDT;          ///< Maximum time difference between planes.
    double fMaxS;           ///< Maximum space separation between wires.
    double fTimeOffset[3];  ///< Time offset corrections (indexed by view-1).
    int fMinViews;          ///< Mininum number of views per space point.
    bool fEnable[3];        ///< Enable flag for three views (indexed by view-1).

    // Geometry and LAr constants.

    const geo::Geometry* fGeom;  ///< Geometry service.
    std::string fGDMLPath;   ///< Used to determine if geometry has changed.
    std::string fROOTPath;   ///< Used to determine if geometry has changed.
    double fWirePitch[3];    ///< Wire pitch (cm) [view-1].
    double fWireOffset[3];   ///< Distance of wire 0 from origin (cm) [view-1].
    double fTheta[3];        ///< theta[view-1].
    double fSinTheta[3];     ///< Sin(theta[view-1]).
    double fCosTheta[3];     ///< Cos(theta[view-1]).
    double fSin[3];          ///< Sin(delta theta) [view-1].

    // From Detector properties.

    double fSamplingRate;    ///< ns/tick.
    int fTriggerOffset;      ///< tick offset (aka presamplings).

    // From LArProperties.

    double fEfield;          ///< Electric field in drift volume (kV/cm).
    double fTemperature;     ///< LAr temperature.
    double fDriftVelocity;   ///< Drift velocity (cm/us).
    double fTimePitch;       ///< Time pitch (cm/tick).

    // Histograms.

    bool fBooked;   // Have histograms been booked yet?
    TH1F* fHDTUV;   // U-V time difference.
    TH1F* fHDTVW;   // V-W time difference.
    TH1F* fHDTWU;   // W-U time difference.
    TH1F* fHS;      // Spatial separation parameter.
    TH1F* fHDTUVC;  // U-V time difference (spatially compatible hits).
    TH1F* fHDTVWC;  // V-W time difference (spatially compatible hits).
    TH1F* fHDTWUC;  // W-U time difference (spatially compatible hits).
    TH1F* fHx;      // X position.
    TH1F* fHy;      // Y position.
    TH1F* fHz;      // Z position.
    TH1F* fHMCdx;   // X residual (reco vs. mc truth).
    TH1F* fHMCdy;   // Y residual (reco vs. mc truth).
    TH1F* fHMCdz;   // Z residual (reco vs. mc truth).

    // Temporary variables.

    mutable art::Handle<std::vector<sim::SimChannel> > fSCHandle;
    mutable std::map<const recob::Hit*, std::vector<const sim::Electrons*> > fElecMap;
    mutable std::map<const recob::Hit*, TVector3> fMCPosMap;
  };
}

#endif
