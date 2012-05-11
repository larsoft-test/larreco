#ifndef Track3DKalmanHit_h
#define Track3DKalmanHit_h
////////////////////////////////////////////////////////////////////////
// Class:       Track3DKalmanHit
// Module Type: producer
// File:        Track3DKalmanHit.h
//
// This class produces RecoBase/Track objects using KalmanFilterAlgorithm.
//
// Configuration parameters:
//
// Hist               - Histogram flag (generate histograms if true).
// UseClusterHits     - Use clustered hits if true, use all hits if false.
// HitModuleLabel     - Module label for unclustered Hits.
// ClusterModuleLabel - Module label for Clusters.
// G4ModuleLabel      - Module label for MC truth particles.
// MaxTcut            - Maximum delta ray energy in Mev for dE/dx.
// KalmanFilterAlg    - Parameter set for KalmanFilterAlg.
// SpacePointAlg      - Parmaeter set for space points.
// SeedFinder         - Parameter set for seed finder.
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "TrackFinder/KalmanFilterAlg.h"
#include "TrackFinder/SeedFinder.h"
#include "TH1F.h"

namespace trkf {

  class Propagator;

  class Track3DKalmanHit : public art::EDProducer {
  public:

    // Copnstructors, destructor.

    explicit Track3DKalmanHit(fhicl::ParameterSet const & pset);
    virtual ~Track3DKalmanHit();

    // Overrides.

    virtual void reconfigure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e);
    virtual void beginJob();
    virtual void endJob();

  private:

    // Fcl parameters.

    bool fHist;                        ///< Make histograms.
    bool fUseClusterHits;              ///< Use cluster hits or all hits?
    std::string fHitModuleLabel;       ///< Unclustered Hits.
    std::string fClusterModuleLabel;   ///< Clustered Hits.
    double fMaxTcut;                   ///< Maximum delta ray energy in MeV for restricted dE/dx.
    double fMinSeedHits;               ///< Minimum number of hits per track seed.
    double fMaxSeedChiDF;              ///< Maximum seed track chisquare/dof.
    double fMinSeedSlope;              ///< Minimum seed slope (dx/dz).

    // Algorithm objects.

    KalmanFilterAlg fKFAlg;            ///< Kalman filter algorithm.
    SeedFinder fSeedFinder;            ///< Seed finder.
    SpacePointAlg fSpacePointAlg;      ///< Space point algorithm.

    /// Propagator.
    const Propagator* fProp;

    // Histograms.

    TH1F* fHIncChisq;   ///< Incremental chisquare.
    TH1F* fHPull;       ///< Hit pull.

    // Statistics.

    int fNumEvent;    ///< Number of events seen.
    int fNumTrack;    ///< Number of tracks produced.

  };

} // namespace trkf

#endif
