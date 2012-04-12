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
/// MaxDT - The maximum time difference (ticks) between any pair of hits.
/// MaxS  - The maximum 3-view wire separation parameter S (cm).
/// MinViews - Minimum number of views to make a space point (2 or 3).
/// EnableU - Use U view hits.
/// EnableV - Use V view hits.
/// EnableW - Use W view hits.
/// Filter - Filter space points flag.
/// Merge - Merge space points flag.
///
/// The parameters fMaxDT and fMaxS are used to implement a notion of whether
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
/// Actual corrected times for hits use time offsets calculated using
/// method fillTimeOffset.  These time offsets include all known geometric
/// and electronic time offsets from Geometry, DetectorProperties, and
/// LArProperties services, in addition to the configuration time offsets.
///
/// If enabled, filtering eliminates multiple space points with similar
/// times on the same wire of the most popluated plane.
///
/// If enabled, merging combines multiple space points with similar 
/// times on the same wire of the most populated plane (potentially
/// producing space points with more hits than the number of planes).
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
#include "RecoBase/SpacePoint.h"
#include "Simulation/SimChannel.h"

class TH1F;
namespace sim {
  class IDE;
}

namespace trkf {

  class SpacePointService {
  public:

    // Constructor.
    SpacePointService(const fhicl::ParameterSet& pset, art::ActivityRegistry& reg);

    // Destructor.
    ~SpacePointService();

    // Configuration Accessors.

    bool filter() const {return fFilter;}
    bool merge() const {return fMerge;}
    double maxDT() const {return fMaxDT;}
    double maxS() const {return fMaxS;}
    int minViews() const {return fMinViews;}
    bool enableU() const {return fEnableU;}
    bool enableV() const {return fEnableV;}
    bool enableW() const {return fEnableW;}

    // Update configuration parameters.
    void reconfigure(const fhicl::ParameterSet& pset);

    // Print constants obtained from geometry and properties services.
    void update() const;

    // Corrected time accessors.
    double correctedTime(const recob::Hit& hit) const;

    // Spatial separation of hits (zero if two or fewer).
    double separation(const art::PtrVector<recob::Hit>& hits) const;

    // Fill a single simple space point using the specified hits.
    // Hits are assumed to be compatible.
    void fillSpacePoint(const art::PtrVector<recob::Hit>& hits,
			recob::SpacePoint& spt) const;

    // Fill a single complex space point using the specified hits.
    // Complex space points allow multiple hits in one plane.
    // Hits are assumed to be compatible.
    void fillComplexSpacePoint(const art::PtrVector<recob::Hit>& hits,
			       recob::SpacePoint& spt) const;

    // Fill a vector of space points from an unsorted collection of hits.
    // Space points are generated for all compatible combinations of hits.
    // Second version of this method allows to override various
    // configuration parameters.
    void makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
			 std::vector<recob::SpacePoint>& spts) const;
    void makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
			 std::vector<recob::SpacePoint>& spts,
			 bool filter,
			 bool merge,
			 double maxDT,
			 double maxS) const;

    // Fill a vector of space points compatible with mc truth information 
    // contained in SimChannels (with or without overriding configuration
    // parameters).
    void makeMCTruthSpacePoints(const art::PtrVector<recob::Hit>& hits,
				std::vector<recob::SpacePoint>& spts,
				const std::vector<const sim::SimChannel*>& simchans) const;
    void makeMCTruthSpacePoints(const art::PtrVector<recob::Hit>& hits,
				std::vector<recob::SpacePoint>& spts,
				const std::vector<const sim::SimChannel*>& simchans,
				bool filter,
				bool merge,
				double maxDT,
				double maxS) const;

  private:

    // This is the real method for calculating space points (each of
    // the public make*SpacePoints methods comes here).
    void makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
			 std::vector<recob::SpacePoint>& spts,
			 const std::vector<const sim::SimChannel*>& simchans,
			 bool useMC,
			 bool filter,
			 bool merge,
			 double maxDT,
			 double maxS) const;

    // Test whether the specified hits are compatible with a space point.
    // The last two arguments can be used to override the default cuts.
    bool compatible(const art::PtrVector<recob::Hit>& hits,
		    bool useMC,
		    double maxDT,
		    double maxS) const;

    // Extended SpacePoint struct.
    // Contains extra information used for filtering (private type).

    //struct SpacePointX : public recob::SpacePoint {
    //  SpacePointX() : goodness(0) {}
    //  double goodness;
    //};

    // Configuration paremeters.

    double fMaxDT;          ///< Maximum time difference between planes.
    double fMaxS;           ///< Maximum space separation between wires.
    int fMinViews;          ///< Mininum number of views per space point.
    bool fEnableU;          ///< Enable flag (U).
    bool fEnableV;          ///< Enable flag (V).
    bool fEnableW;          ///< Enable flag (W).
    bool fFilter;           ///< Filter flag.
    bool fMerge;            ///< Merge flag.

    // Temporary variables.

    struct HitMCInfo
    {
      std::vector<int> trackIDs;       ///< Parent trackIDs.
      std::vector<double> xyz;         ///< Location of ionization (all tracks).
      std::vector<const recob::Hit*> pchit;   ///< Pointer to nearest neighbor hit (indexed by plane).
      std::vector<double> dist2;              ///< Distance to nearest neighbor hit (indexed by plane).
    };
    mutable std::map<const recob::Hit*, HitMCInfo> fHitMCMap;
  };
}

#endif
