#ifndef TrackKalmanCheater_h
#define TrackKalmanCheater_h
////////////////////////////////////////////////////////////////////////
// Class:       TrackKalmanCheater
// Module Type: producer
// File:        TrackKalmanCheater.h
//
// This class produces RecoBase/Track objects using KalmanFilterService.
// MC truth information is used to associate Hits used as input to
// the Kalman filter.
//
// Configuration parameters:
//
// Hist               - Histogram flag (generate histograms if true).
// UseClusterHits     - Use clustered hits if true, use all hits if false.
// HitModuleLabel     - Module label for unclustered Hits.
// ClusterModuleLabel - Module label for Clusters.
// G4ModuleLabel      - Module label for MC truth particles.
// MaxIncChisq        - Maximum incremental chisquare cut.
// MaxTcut            - Maximum delta ray energy in Mev for dE/dx.
// KalmanFilterAlg    - Parameter set for KalmanFilterAlg.
// SpacePointAlg      - Parmaeter set for space points.
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "TrackFinder/KalmanFilterAlg.h"
#include "TH1F.h"

namespace trkf {

  class Propagator;

  class TrackKalmanCheater : public art::EDProducer {
  public:

    // Copnstructors, destructor.

    explicit TrackKalmanCheater(fhicl::ParameterSet const & pset);
    virtual ~TrackKalmanCheater();

    // Overrides.

    virtual void reconfigure(fhicl::ParameterSet const & pset);
    virtual void produce(art::Event & e);
    virtual void beginJob();
    virtual void endJob();

  private:

    // Fcl parameters.

    bool fHist;                        ///< Make histograms.
    KalmanFilterAlg fKFAlg;            ///< Kalman filter algorithm.
    SpacePointAlg fSpacePointAlg;      ///< Space point algorithm.
    bool fUseClusterHits;              ///< Use cluster hits or all hits?
    std::string fHitModuleLabel;       ///< Unclustered Hits.
    std::string fClusterModuleLabel;   ///< Clustered Hits.
    std::string fG4ModuleLabel;        ///< For SimChannel.
    double fMaxIncChisq;               ///< Maximum incremental chisquare.
    double fMaxTcut;                   ///< Maximum delta ray energy in MeV for restricted dE/dx.

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
