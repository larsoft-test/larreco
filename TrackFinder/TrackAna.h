#ifndef TRACKANA_H
#define TRACKANA_H

//
// Name: TrackAna.h
//
// Purpose: Header file for module TrackAna.
//
// Created: 2-Aug-2011  H. Greenlee
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "TH2F.h"

namespace trkf {

  class TrackAna : public art::EDAnalyzer
  {
  public:
 
    // Constructors, destructor

    explicit TrackAna(fhicl::ParameterSet const& pset);
    virtual ~TrackAna();

    // Book histograms.

    void bookHistograms(bool mc);

    // Overrides.

    void analyze(const art::Event& evt);

  private:

    // Fcl Attributes.

    std::string fTrackModuleLabel;
    std::string fG4ModuleLabel;
    double fMinMCKE;           // Minimum MC particle kinetic energy (GeV).
    double fMatchColinearity;  // Minimum matching colinearity.
    double fMatchDisp;         // Maximum matching displacement.

    // Histograms.

    bool fBooked;      // Have histograms been booked yet?
    TH1F* fHstartx;    // Starting x position.
    TH1F* fHstarty;    // Starting y position.
    TH1F* fHstartz;    // Starting z position.
    TH1F* fHstartd;    // Starting distance to boundary.
    TH1F* fHendx;      // Ending x position.
    TH1F* fHendy;      // Ending y position.
    TH1F* fHendz;      // Ending z position.
    TH1F* fHendd;      // Ending distance to boundary.
    TH1F* fHtheta;     // Theta.
    TH1F* fHphi;       // Phi.
    TH1F* fHmom;       // Momentum.
    TH2F* fHduvcosth;  // 2D mc vs. data matching, duv vs. cos(theta).
    TH1F* fHcosth;     // 1D direction matching, cos(theta).
    TH1F* fHmcu;       // 1D endpoint truth u.
    TH1F* fHmcv;       // 1D endpoint truth v.
    TH1F* fHmcw;       // 1D endpoint truth w.
    TH1F* fHupull;     // 1D endpoint u pull.
    TH1F* fHvpull;     // 1D endpoint v pull.
    TH1F* fHmcdudw;    // Truth du/dw.
    TH1F* fHmcdvdw;    // Truth dv/dw.
    TH1F* fHdudwpull;  // du/dw pull.
    TH1F* fHdvdwpull;  // dv/dw pull.
    TH2F* fHpvsp;      // MC vs. reco momentum
    TH2F* fHpvspc;     // MC vs. reco momentum (contained tracks).
    TH1F* fHdp;        // Momentum difference.
    TH1F* fHdpc;       // Momentum difference (contained tracks).
    TH1F* fHppull;     // Momentum pull.
    TH1F* fHppullc;    // Momentum pull (contained tracks).

    // Statistics.

    int fNumEvent;
  };
}

#endif // TRACKANA_H
