#ifndef TRACKANA_H
#define TRACKANA_H

//
// Name: TrackAna.h
//
// Purpose: Header file for module TrackAna.
//
// Created: 2-Aug-2011  H. Greenlee
//

#include <map>
#include "art/Framework/Core/EDAnalyzer.h"
#include "TH2F.h"

namespace trkf {

  class TrackAna : public art::EDAnalyzer
  {
  public:

    // Embedded structs.

    // Struct for histograms that depend on reco track only.

    struct RecoHists
    {
      // Constructors.

      RecoHists();
      RecoHists(const std::string& subdir);

      // Pure reco track histograms.

      TH1F* fHstartx;      // Starting x position.
      TH1F* fHstarty;      // Starting y position.
      TH1F* fHstartz;      // Starting z position.
      TH1F* fHstartd;      // Starting distance to boundary.
      TH1F* fHendx;        // Ending x position.
      TH1F* fHendy;        // Ending y position.
      TH1F* fHendz;        // Ending z position.
      TH1F* fHendd;        // Ending distance to boundary.
      TH1F* fHtheta;       // Theta.
      TH1F* fHphi;         // Phi.
      TH1F* fHtheta_xz;    // Theta_xz.
      TH1F* fHtheta_yz;    // Theta_yz.
      TH1F* fHmom;         // Momentum.
      TH1F* fHlen;         // Length.
    };

    // Struct for mc particles and mc-matched tracks.

    struct MCHists
    {
      // Constructors.

      MCHists();
      MCHists(const std::string& subdir);

      // Reco-MC matching.

      TH2F* fHduvcosth;    // 2D mc vs. data matching, duv vs. cos(theta).
      TH1F* fHcosth;       // 1D direction matching, cos(theta).
      TH1F* fHmcu;         // 1D endpoint truth u.
      TH1F* fHmcv;         // 1D endpoint truth v.
      TH1F* fHmcw;         // 1D endpoint truth w.
      TH1F* fHupull;       // 1D endpoint u pull.
      TH1F* fHvpull;       // 1D endpoint v pull.
      TH1F* fHmcdudw;      // Truth du/dw.
      TH1F* fHmcdvdw;      // Truth dv/dw.
      TH1F* fHdudwpull;    // du/dw pull.
      TH1F* fHdvdwpull;    // dv/dw pull.

      // Histograms for matched tracks.

      TH1F* fHstartdx;     // Start dx.
      TH1F* fHstartdy;     // Start dy.
      TH1F* fHstartdz;     // Start dz.
      TH1F* fHenddx;       // End dx.
      TH1F* fHenddy;       // End dy.
      TH1F* fHenddz;       // End dz.
      TH2F* fHlvsl;        // MC vs. reco length.
      TH1F* fHdl;          // Delta(length).
      TH2F* fHpvsp;        // MC vs. reco momentum.
      TH2F* fHpvspc;       // MC vs. reco momentum (contained tracks).
      TH1F* fHdp;          // Momentum difference.
      TH1F* fHdpc;         // Momentum difference (contained tracks).
      TH1F* fHppull;       // Momentum pull.
      TH1F* fHppullc;      // Momentum pull (contained tracks).

      // Pure MC particle histograms (efficiency denominator).

      TH1F* fHmcstartx;    // Starting x position.
      TH1F* fHmcstarty;    // Starting y position.
      TH1F* fHmcstartz;    // Starting z position.
      TH1F* fHmcendx;      // Ending x position.
      TH1F* fHmcendy;      // Ending y position.
      TH1F* fHmcendz;      // Ending z position.
      TH1F* fHmctheta;     // Theta.
      TH1F* fHmcphi;       // Phi.
      TH1F* fHmctheta_xz;  // Theta_xz.
      TH1F* fHmctheta_yz;  // Theta_yz.
      TH1F* fHmcmom;       // Momentum.
      TH1F* fHmclen;       // Length.

      // Histograms for well-reconstructed matched tracks (efficiency numerator).

      TH1F* fHgstartx;     // Starting x position.
      TH1F* fHgstarty;     // Starting y position.
      TH1F* fHgstartz;     // Starting z position.
      TH1F* fHgendx;       // Ending x position.
      TH1F* fHgendy;       // Ending y position.
      TH1F* fHgendz;       // Ending z position.
      TH1F* fHgtheta;      // Theta.
      TH1F* fHgphi;        // Phi.
      TH1F* fHgtheta_xz;   // Theta_xz.
      TH1F* fHgtheta_yz;   // Theta_yz.
      TH1F* fHgmom;        // Momentum.
      TH1F* fHglen;        // Length.

      // Efficiency histograms.

      TH1F* fHestartx;     // Starting x position.
      TH1F* fHestarty;     // Starting y position.
      TH1F* fHestartz;     // Starting z position.
      TH1F* fHeendx;       // Ending x position.
      TH1F* fHeendy;       // Ending y position.
      TH1F* fHeendz;       // Ending z position.
      TH1F* fHetheta;      // Theta.
      TH1F* fHephi;        // Phi.
      TH1F* fHetheta_xz;   // Theta_xz.
      TH1F* fHetheta_yz;   // Theta_yz.
      TH1F* fHemom;        // Momentum.
      TH1F* fHelen;        // Length.
    };

    // Constructors, destructor

    explicit TrackAna(fhicl::ParameterSet const& pset);
    virtual ~TrackAna();

    // Overrides.

    void analyze(const art::Event& evt);
    void endJob();

  private:

    // Fcl Attributes.

    std::string fTrackModuleLabel;
    std::string fG4ModuleLabel;
    double fMinMCKE;           // Minimum MC particle kinetic energy (GeV).
    double fMatchColinearity;  // Minimum matching colinearity.
    double fMatchDisp;         // Maximum matching displacement.

    // Histograms.

    std::map<int, MCHists> fMCHistMap;       // Indexed by pdg id.
    std::map<int, RecoHists> fRecoHistMap;   // Indexed by pdg id.

    // Statistics.

    int fNumEvent;
  };
}

#endif // TRACKANA_H
