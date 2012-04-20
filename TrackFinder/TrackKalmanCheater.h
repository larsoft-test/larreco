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
// Generated at Wed Mar 28 13:43:51 2012 by Herbert Greenlee using artmod
// from art v1_00_11.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"

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

    bool fUseClusterHits;              ///< Use cluster hits or all hits?
    std::string fGenModuleLabel;       ///< MC truth.
    std::string fHitModuleLabel;       ///< Unclustered Hits.
    std::string fClusterModuleLabel;   ///< Clustered Hits.
    std::string fG4ModuleLabel;        ///< For SimChannel.

    /// Propagator.
    const Propagator* prop;

    // Statistics.

    int fNumEvent;    ///< Number of events seen.
    int fNumTrack;    ///< Number of tracks produced.

  };

} // namespace trkf

#endif
