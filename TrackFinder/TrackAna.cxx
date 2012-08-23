//
// Name: TrackAna.cxx
//
// Purpose: Implementation file for module TrackAna.
//
// Created: 26-Jun-2012  H. Greenlee
//

#include <iostream>
#include <sstream>
#include <cmath>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Geometry/geo.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "TrackFinder/TrackAna.h"
#include "art/Framework/Principal/Event.h"
#include "RecoBase/Track.h"
#include "MCCheater/BackTracker.h"
#include "TFile.h"

namespace {

  // Local functions.

  // Calculate distance to boundary.

  double bdist(const TVector3& pos, unsigned int tpc = 0, unsigned int cstat = 0)
  {
    // Get geometry.

    art::ServiceHandle<geo::Geometry> geom;

    double d1 = pos.X();                             // Distance to right side (wires).
    double d2 = 2.*geom->DetHalfWidth() - pos.X();   // Distance to left side (cathode).
    double d3 = pos.Y() + geom->DetHalfHeight();     // Distance to bottom.
    double d4 = geom->DetHalfHeight() - pos.Y();     // Distance to top.
    double d5 = pos.Z();                             // Distance to front.
    double d6 = geom->DetLength() - pos.Z();         // Distance to back.

    double result = std::min(std::min(std::min(std::min(std::min(d1, d2), d3), d4), d5), d6);
    return result;
  }

  // Length of reconstructed track.

  double length(const recob::Track& track)
  {
    double result = 0.;
    TVector3 disp = track.LocationAtPoint(0);
    int n = track.NumberTrajectoryPoints();

    for(int i = 1; i < n; ++i) {
      const TVector3& pos = track.LocationAtPoint(i);
      disp -= pos;
      result += disp.Mag();
      disp = pos;
    }

    return result;
  }

  // Length of MC particle.

  double length(const sim::Particle& part, 
		TVector3& start, TVector3& end,
		unsigned int tpc = 0, unsigned int cstat = 0)
  {
    // Get geometry.

    art::ServiceHandle<geo::Geometry> geom;

    // Get fiducial volume boundary.

    double xmin = 0.;
    double xmax = 2.*geom->DetHalfWidth();
    double ymin = -geom->DetHalfHeight();
    double ymax = geom->DetHalfHeight();
    double zmin = 0.;
    double zmax = geom->DetLength();

    double result = 0.;
    TVector3 disp;
    int n = part.NumberTrajectoryPoints();
    bool first = true;

    for(int i = 0; i < n; ++i) {
      const TVector3& pos = part.Position(i).Vect();
      if(pos.X() >= xmin &&
	 pos.X() <= xmax &&
	 pos.Y() >= ymin &&
	 pos.Y() <= ymax &&
	 pos.Z() >= zmin &&
	 pos.Z() <= zmax) {
	if(first)
	  start = pos;
	else {
	  disp -= pos;
	  result += disp.Mag();
	}
	first = false;
	disp = pos;
	end = pos;
      }
    }

    return result;
  }

  // Fill efficiency histogram assuming binomial errors.

  void effcalc(const TH1* hnum, const TH1* hden, TH1* heff)
  {
    int nbins = hnum->GetNbinsX();
    assert(nbins == hden->GetNbinsX());
    assert(nbins == heff->GetNbinsX());

    // Loop over bins, including underflow and overflow.

    for(int ibin = 0; ibin <= nbins+1; ++ibin) {
      double num = hnum->GetBinContent(ibin);
      double den = hden->GetBinContent(ibin);
      if(den == 0.) {
	heff->SetBinContent(ibin, 0.);
	heff->SetBinError(ibin, 0.);
      }
      else {
	double eff = num / den;
	double err = std::sqrt(eff * (1.-eff) / den);
	heff->SetBinContent(ibin, eff);
	heff->SetBinError(ibin, err);
	assert(eff >= 0. && eff <= 1.);
      }
    }

    heff->SetMinimum(0.);
    heff->SetMaximum(1.05);
    heff->SetMarkerStyle(20);
  }
}

namespace trkf {

  // RecoHists methods.

  TrackAna::RecoHists::RecoHists() :
    //
    // Purpose: Default constructor.
    //
    fHstartx(0),
    fHstarty(0),
    fHstartz(0),
    fHstartd(0),
    fHendx(0),
    fHendy(0),
    fHendz(0),
    fHendd(0),
    fHtheta(0),
    fHphi(0),
    fHtheta_xz(0),
    fHtheta_yz(0),
    fHmom(0),
    fHlen(0)
  {}

  TrackAna::RecoHists::RecoHists(const std::string& subdir)
  //
  // Purpose: Initializing constructor.
  //
  {
    // Make sure all histogram pointers are initially zero.

    *this = RecoHists();

    // Get services.

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;

    // Make histogram directory.

    art::TFileDirectory topdir = tfs->mkdir("trkana", "TrackAna histograms");
    art::TFileDirectory dir = topdir.mkdir(subdir);

    // Book histograms.

    fHstartx = dir.make<TH1F>("xstart", "X Start Position",
			      100, 0., 2.*geom->DetHalfWidth());
    fHstarty = dir.make<TH1F>("ystart", "Y Start Position",
			      100, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHstartz = dir.make<TH1F>("zstart", "Z Start Position",
			      100, 0., geom->DetLength());
    fHstartd = dir.make<TH1F>("dstart", "Start Position Distance to Boundary",
			      100, -10., geom->DetHalfWidth());
    fHendx = dir.make<TH1F>("xend", "X End Position",
			    100, 0., 2.*geom->DetHalfWidth());
    fHendy = dir.make<TH1F>("yend", "Y End Position",
			    100, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHendz = dir.make<TH1F>("zend", "Z End Position",
			    100, 0., geom->DetLength());
    fHendd = dir.make<TH1F>("dend", "End Position Distance to Boundary",
			    100, -10., geom->DetHalfWidth());
    fHtheta = dir.make<TH1F>("theta", "Theta", 100, 0., 3.142);
    fHphi = dir.make<TH1F>("phi", "Phi", 100, -3.142, 3.142);
    fHtheta_xz = dir.make<TH1F>("theta_xz", "Theta_xz", 100, -3.142, 3.142);
    fHtheta_yz = dir.make<TH1F>("theta_yz", "Theta_yz", 100, -3.142, 3.142);
    fHmom = dir.make<TH1F>("mom", "Momentum", 100, 0., 10.);
    fHlen = dir.make<TH1F>("len", "Track Length", 100, 0., 1.1 * geom->DetLength());
  }

  // MCHists methods.

  TrackAna::MCHists::MCHists() :
    //
    // Purpose: Default constructor.
    //
    fHduvcosth(0),
    fHcosth(0),
    fHmcu(0),
    fHmcv(0),
    fHmcw(0),
    fHupull(0),
    fHvpull(0),
    fHmcdudw(0),
    fHmcdvdw(0),
    fHdudwpull(0),
    fHdvdwpull(0),
    fHstartdx(0),
    fHstartdy(0),
    fHstartdz(0),
    fHenddx(0),
    fHenddy(0),
    fHenddz(0),
    fHlvsl(0),
    fHdl(0),
    fHpvsp(0),
    fHpvspc(0),
    fHdp(0),
    fHdpc(0),
    fHppull(0),
    fHppullc(0),
    fHmcstartx(0),
    fHmcstarty(0),
    fHmcstartz(0),
    fHmcendx(0),
    fHmcendy(0),
    fHmcendz(0),
    fHmctheta(0),
    fHmcphi(0),
    fHmctheta_xz(0),
    fHmctheta_yz(0),
    fHmcmom(0),
    fHmclen(0),
    fHgstartx(0),
    fHgstarty(0),
    fHgstartz(0),
    fHgendx(0),
    fHgendy(0),
    fHgendz(0),
    fHgtheta(0),
    fHgphi(0),
    fHgtheta_xz(0),
    fHgtheta_yz(0),
    fHgmom(0),
    fHglen(0),
    fHestartx(0),
    fHestarty(0),
    fHestartz(0),
    fHeendx(0),
    fHeendy(0),
    fHeendz(0),
    fHetheta(0),
    fHephi(0),
    fHetheta_xz(0),
    fHetheta_yz(0),
    fHemom(0),
    fHelen(0)
  {}

  TrackAna::MCHists::MCHists(const std::string& subdir)
  //
  // Purpose: Initializing constructor.
  //
  {
    // Make sure all histogram pointers are initially zero.

    *this = MCHists();

    // Get services.

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<art::TFileService> tfs;

    // Make histogram directory.

    art::TFileDirectory topdir = tfs->mkdir("trkana", "TrackAna histograms");
    art::TFileDirectory dir = topdir.mkdir(subdir);

    // Book histograms.

    fHduvcosth = dir.make<TH2F>("duvcosth", "Delta(uv) vs. Colinearity", 
				100, 0.95, 1., 100, 0., 1.);
    fHcosth = dir.make<TH1F>("colin", "Colinearity", 100, 0.95, 1.);
    fHmcu = dir.make<TH1F>("mcu", "MC Truth U", 100, -1., 1.);
    fHmcv = dir.make<TH1F>("mcv", "MC Truth V", 100, -1., 1.);
    fHmcw = dir.make<TH1F>("mcw", "MC Truth W", 100, -1., 1.);
    fHupull = dir.make<TH1F>("dupull", "U Pull", 100, -10., 10.);
    fHvpull = dir.make<TH1F>("dvpull", "V Pull", 100, -10., 10.);
    fHmcdudw = dir.make<TH1F>("mcdudw", "MC Truth U Slope", 100, -0.2, 0.2);
    fHmcdvdw = dir.make<TH1F>("mcdvdw", "MV Truth V Slope", 100, -0.2, 0.2);
    fHdudwpull = dir.make<TH1F>("dudwpull", "U Slope Pull", 100, -10., 10.);
    fHdvdwpull = dir.make<TH1F>("dvdwpull", "V Slope Pull", 100, -10., 10.);
    fHstartdx = dir.make<TH1F>("startdx", "Start Delta x", 100, -1., 1.);
    fHstartdy = dir.make<TH1F>("startdy", "Start Delta y", 100, -1., 1.);
    fHstartdz = dir.make<TH1F>("startdz", "Start Delta z", 100, -1., 1.);
    fHenddx = dir.make<TH1F>("enddx", "End Delta x", 100, -10., 10.);
    fHenddy = dir.make<TH1F>("enddy", "End Delta y", 100, -10., 10.);
    fHenddz = dir.make<TH1F>("enddz", "End Delta z", 100, -20., 20.);
    fHlvsl = dir.make<TH2F>("lvsl", "Reco Length vs. MC Truth Length",
			    100, 0., 1.1 * geom->DetLength(), 100, 0., 1.1 * geom->DetLength());
    fHdl = dir.make<TH1F>("dl", "Track Length Minus MC Particle Length", 100, -50., 50.);
    fHpvsp = dir.make<TH2F>("pvsp", "Reco Momentum vs. MC Truth Momentum",
			    100, 0., 5., 100, 0., 5.);
    fHpvspc = dir.make<TH2F>("pvspc", "Reco Momentum vs. MC Truth Momentum (Contained Tracks)",
			     100, 0., 5., 100, 0., 5.);
    fHdp = dir.make<TH1F>("dp", "Reco-MC Momentum Difference", 100, -5., 5.);
    fHdpc = dir.make<TH1F>("dpc", "Reco-MC Momentum Difference (Contained Tracks)",
			   100, -5., 5.);
    fHppull = dir.make<TH1F>("ppull", "Momentum Pull", 100, -10., 10.);
    fHppullc = dir.make<TH1F>("ppullc", "Momentum Pull (Contained Tracks)", 100, -10., 10.);

    fHmcstartx = dir.make<TH1F>("mcxstart", "MC X Start Position",
				10, 0., 2.*geom->DetHalfWidth());
    fHmcstarty = dir.make<TH1F>("mcystart", "MC Y Start Position",
				10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHmcstartz = dir.make<TH1F>("mczstart", "MC Z Start Position",
				10, 0., geom->DetLength());
    fHmcendx = dir.make<TH1F>("mcxend", "MC X End Position",
			      10, 0., 2.*geom->DetHalfWidth());
    fHmcendy = dir.make<TH1F>("mcyend", "MC Y End Position",
			      10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHmcendz = dir.make<TH1F>("mczend", "MC Z End Position",
			      10, 0., geom->DetLength());
    fHmctheta = dir.make<TH1F>("mctheta", "MC Theta", 20, 0., 3.142);
    fHmcphi = dir.make<TH1F>("mcphi", "MC Phi", 10, -3.142, 3.142);
    fHmctheta_xz = dir.make<TH1F>("mctheta_xz", "MC Theta_xz", 40, -3.142, 3.142);
    fHmctheta_yz = dir.make<TH1F>("mctheta_yz", "MC Theta_yz", 40, -3.142, 3.142);
    fHmcmom = dir.make<TH1F>("mcmom", "MC Momentum", 10, 0., 10.);
    fHmclen = dir.make<TH1F>("mclen", "MC Particle Length", 10, 0., 1.1 * geom->DetLength());

    fHgstartx = dir.make<TH1F>("gxstart", "Good X Start Position",
			       10, 0., 2.*geom->DetHalfWidth());
    fHgstarty = dir.make<TH1F>("gystart", "Good Y Start Position",
			       10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHgstartz = dir.make<TH1F>("gzstart", "Good Z Start Position",
			       10, 0., geom->DetLength());
    fHgendx = dir.make<TH1F>("gxend", "Good X End Position",
			     10, 0., 2.*geom->DetHalfWidth());
    fHgendy = dir.make<TH1F>("gyend", "Good Y End Position",
			     10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHgendz = dir.make<TH1F>("gzend", "Good Z End Position",
			     10, 0., geom->DetLength());
    fHgtheta = dir.make<TH1F>("gtheta", "Good Theta", 20, 0., 3.142);
    fHgphi = dir.make<TH1F>("gphi", "Good Phi", 10, -3.142, 3.142);
    fHgtheta_xz = dir.make<TH1F>("gtheta_xz", "Good Theta_xz", 40, -3.142, 3.142);
    fHgtheta_yz = dir.make<TH1F>("gtheta_yz", "Good Theta_yz", 40, -3.142, 3.142);
    fHgmom = dir.make<TH1F>("gmom", "Good Momentum", 10, 0., 10.);
    fHglen = dir.make<TH1F>("glen", "Good Particle Length", 10, 0., 1.1 * geom->DetLength());

    fHestartx = dir.make<TH1F>("exstart", "Efficiency vs. X Start Position",
			       10, 0., 2.*geom->DetHalfWidth());
    fHestarty = dir.make<TH1F>("eystart", "Efficiency vs. Y Start Position",
			       10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHestartz = dir.make<TH1F>("ezstart", "Efficiency vs. Z Start Position",
			       10, 0., geom->DetLength());
    fHeendx = dir.make<TH1F>("exend", "Efficiency vs. X End Position",
			     10, 0., 2.*geom->DetHalfWidth());
    fHeendy = dir.make<TH1F>("eyend", "Efficiency vs. Y End Position",
			     10, -geom->DetHalfHeight(), geom->DetHalfHeight());
    fHeendz = dir.make<TH1F>("ezend", "Efficiency vs. Z End Position",
			     10, 0., geom->DetLength());
    fHetheta = dir.make<TH1F>("etheta", "Efficiency vs. Theta", 20, 0., 3.142);
    fHephi = dir.make<TH1F>("ephi", "Efficiency vs. Phi", 10, -3.142, 3.142);
    fHetheta_xz = dir.make<TH1F>("etheta_xz", "Efficiency vs. Theta_xz", 40, -3.142, 3.142);
    fHetheta_yz = dir.make<TH1F>("etheta_yz", "Efficiency vs. Theta_yz", 40, -3.142, 3.142);
    fHemom = dir.make<TH1F>("emom", "Efficiency vs. Momentum", 10, 0., 10.);
    fHelen = dir.make<TH1F>("elen", "Efficiency vs. Particle Length",
			    10, 0., 1.1 * geom->DetLength());
  }

  TrackAna::TrackAna(const fhicl::ParameterSet& pset) :
    //
    // Purpose: Constructor.
    //
    // Arguments: pset - Module parameters.
    //
    fTrackModuleLabel(pset.get<std::string>("TrackModuleLabel")),
    fG4ModuleLabel(pset.get<std::string>("G4ModuleLabel")),
    fMinMCKE(pset.get<double>("MinMCKE")),
    fMatchColinearity(pset.get<double>("MatchColinearity")),
    fMatchDisp(pset.get<double>("MatchDisp")),
    fNumEvent(0)
  {

    // Report.

    mf::LogInfo("TrackAna") 
      << "TrackAna configured with the following parameters:\n"
      << "  TrackModuleLabel = " << fTrackModuleLabel << "\n"
      << "  G4ModuleLabel = " << fG4ModuleLabel << "\n"
      << "  MinMCKE = " << fMinMCKE;
  }

  TrackAna::~TrackAna()
  //
  // Purpose: Destructor.
  //
  {}

  void TrackAna::analyze(const art::Event& evt)
  //
  // Purpose: Analyze method.
  //
  // Arguments: event - Art event.
  //
  {
    ++fNumEvent;

    // Make sure histograms are booked.

    bool mc = !evt.isRealData();

    // Get mc particles.

    sim::ParticleList plist;
    std::vector<const sim::Particle*> plist2;
    plist2.reserve(plist.size());

    if(mc) {
      art::ServiceHandle<cheat::BackTracker> bt;
      plist = bt->ParticleList();

      // Loop over mc particles, and fill histograms that depend only
      // on mc particles.  Also, fill a secondary list of mc particles
      // that pass various selection criteria.

      for(sim::ParticleList::const_iterator ipart = plist.begin();
	  ipart != plist.end(); ++ipart) {
	const sim::Particle* part = (*ipart).second;
	assert(part != 0);
	int pdg = part->PdgCode();

	// Ignore everything except stable charged nonshowering particles.

	int apdg = std::abs(pdg);
	if(apdg == 13 ||     // Muon
	   apdg == 211 ||    // Charged pion
	   apdg == 321 ||    // Charged kaon
	   apdg == 2212) {   // (Anti)proton

	  // Apply minimum energy cut.

	  if(part->E() >= 0.001*part->Mass() + fMinMCKE) {

	    // This is a good mc particle (capable of making a track).

	    plist2.push_back(part);

	    // Fill histograms.

	    if(fMCHistMap.count(pdg) == 0) {
	      std::ostringstream ostr;
	      ostr << "MC" << (pdg > 0 ? "Pos" : "Neg") << std::abs(pdg);
	      fMCHistMap[pdg] = MCHists(ostr.str());
	    }
	    const MCHists& mchists = fMCHistMap[pdg];

	    TVector3 mcstart;
	    TVector3 mcend;
	    double plen = length(*part, mcstart, mcend);
	    double mctheta_xz = std::atan2(part->Px(), part->Pz());
	    double mctheta_yz = std::atan2(part->Py(), part->Pz());

	    mchists.fHmcstartx->Fill(mcstart.X());
	    mchists.fHmcstarty->Fill(mcstart.Y());
	    mchists.fHmcstartz->Fill(mcstart.Z());
	    mchists.fHmcendx->Fill(mcend.X());
	    mchists.fHmcendy->Fill(mcend.Y());
	    mchists.fHmcendz->Fill(mcend.Z());
	    mchists.fHmctheta->Fill(part->Momentum().Theta());
	    mchists.fHmcphi->Fill(part->Momentum().Phi());
	    mchists.fHmctheta_xz->Fill(mctheta_xz);
	    mchists.fHmctheta_yz->Fill(mctheta_yz);
	    mchists.fHmcmom->Fill(part->Momentum().Vect().Mag());
	    mchists.fHmclen->Fill(plen);
	  }
	}
      }
    }

    // Get tracks.

    art::Handle< std::vector<recob::Track> > trackh;
    evt.getByLabel(fTrackModuleLabel, trackh);
    if(trackh.isValid()) {

      // Loop over tracks.

      int ntrack = trackh->size();
      for(int i = 0; i < ntrack; ++i) {
	art::Ptr<recob::Track> ptrack(trackh, i);
	const recob::Track& track = *ptrack;

	// Fill histograms involving reco tracks only.

	int ntraj = track.NumberTrajectoryPoints();
	if(ntraj > 0) {
	  const TVector3& pos = track.Vertex();
	  const TVector3& dir = track.VertexDirection();
	  const TVector3& end = track.End();

	  double dpos = bdist(pos);
	  double dend = bdist(end);
	  double tlen = length(track);
	  double theta_xz = std::atan2(dir.X(), dir.Z());
	  double theta_yz = std::atan2(dir.Y(), dir.Z());

	  if(fRecoHistMap.count(0) == 0)
	    fRecoHistMap[0] = RecoHists("Reco");
	  const RecoHists& rhists = fRecoHistMap[0];

	  rhists.fHstartx->Fill(pos.X());
	  rhists.fHstarty->Fill(pos.Y());
	  rhists.fHstartz->Fill(pos.Z());
	  rhists.fHstartd->Fill(dpos);
	  rhists.fHendx->Fill(end.X());
	  rhists.fHendy->Fill(end.Y());
	  rhists.fHendz->Fill(end.Z());
	  rhists.fHendd->Fill(dend);
	  rhists.fHtheta->Fill(dir.Theta());
	  rhists.fHphi->Fill(dir.Phi());
	  rhists.fHtheta_xz->Fill(theta_xz);
	  rhists.fHtheta_yz->Fill(theta_yz);

	  double mom = 0.;
	  if(track.NumberFitMomentum() > 0)
	    mom = track.VertexMomentum();
	  rhists.fHmom->Fill(mom);
	  rhists.fHlen->Fill(tlen);

	  // Calculate the global-to-local rotation matrix.

	  TMatrixD rot(3,3);
	  track.GlobalToLocalRotationAtPoint(0, rot);

	  // Get covariance matrix.

	  const TMatrixD& cov = track.VertexCovariance();

	  // Loop over track-like mc particles.

	  for(std::vector<const sim::Particle*>::const_iterator ipart = plist2.begin();
	      ipart != plist2.end(); ++ipart) {
	    const sim::Particle* part = *ipart;
	    assert(part != 0);
	    int pdg = part->PdgCode();
	    assert(fMCHistMap.count(pdg) > 0);
	    const MCHists& mchists = fMCHistMap[pdg];

	    // Get the momentum and displacement of this mc
	    // particle in the global coordinate system.

	    TVector3 mcmom = part->Momentum().Vect();
	    TVector3 mcpos = part->Position().Vect() - pos;

	    // Rotate the momentum and position to the
	    // track-local coordinate system.

	    TVector3 mcmoml = rot * mcmom;
	    TVector3 mcposl = rot * mcpos;

	    double colinearity = mcmoml.Z() / mcmoml.Mag();

	    double u = mcposl.X();
	    double v = mcposl.Y();
	    double w = mcposl.Z();

	    double pu = mcmoml.X();
	    double pv = mcmoml.Y();
	    double pw = mcmoml.Z();

	    double dudw = pu / pw;
	    double dvdw = pv / pw;

	    double u0 = u - w * dudw;
	    double v0 = v - w * dvdw;
	    double uv0 = std::sqrt(u0*u0 + v0*v0);

	    mchists.fHduvcosth->Fill(colinearity, uv0);
	    if(std::abs(uv0) < fMatchDisp) {

	      // Fill slope matching histograms.

	      mchists.fHmcdudw->Fill(dudw);
	      mchists.fHmcdvdw->Fill(dvdw);
	      mchists.fHdudwpull->Fill(dudw / std::sqrt(cov(2,2)));
	      mchists.fHdvdwpull->Fill(dvdw / std::sqrt(cov(3,3)));
	    }
	    mchists.fHcosth->Fill(colinearity);
	    if(colinearity > fMatchColinearity) {

	      // Fill displacement matching histograms.

	      mchists.fHmcu->Fill(u0);
	      mchists.fHmcv->Fill(v0);
	      mchists.fHmcw->Fill(w);
	      mchists.fHupull->Fill(u0 / std::sqrt(cov(0,0)));
	      mchists.fHvpull->Fill(v0 / std::sqrt(cov(1,1)));

	      if(std::abs(uv0) < fMatchDisp) {

		// Fill matching histograms.

		TVector3 mcstart;
		TVector3 mcend;
		double plen = length(*part, mcstart, mcend);
		double mctheta_xz = std::atan2(part->Px(), part->Pz());
		double mctheta_yz = std::atan2(part->Py(), part->Pz());

		mchists.fHstartdx->Fill(pos.X() - mcstart.X());
		mchists.fHstartdy->Fill(pos.Y() - mcstart.Y());
		mchists.fHstartdz->Fill(pos.Z() - mcstart.Z());
		mchists.fHenddx->Fill(end.X() - mcend.X());
		mchists.fHenddy->Fill(end.Y() - mcend.Y());
		mchists.fHenddz->Fill(end.Z() - mcend.Z());
		mchists.fHlvsl->Fill(plen, tlen);
		mchists.fHdl->Fill(tlen - plen);
		mchists.fHpvsp->Fill(mcmom.Mag(), mom);
		double dp = mom - mcmom.Mag();
		mchists.fHdp->Fill(dp);
		mchists.fHppull->Fill(dp / std::sqrt(cov(4,4)));
		if(std::abs(dpos) >= 5. && std::abs(dend) >= 5.) {
		  mchists.fHpvspc->Fill(mcmom.Mag(), mom);
		  mchists.fHdpc->Fill(dp);
		  mchists.fHppullc->Fill(dp / std::sqrt(cov(4,4)));
		}

		// Count this track as well-reconstructed if it is matched to an
		// mc particle (which it is if get here), and if in addition the
		// starting position (w) matches and the reconstructed track length
		// is more than 0.5 of the mc particle trajectory length.

		bool good = std::abs(w) <= fMatchDisp &&
		  tlen > 0.5 * plen;
		if(good) {
		  mchists.fHgstartx->Fill(mcstart.X());
		  mchists.fHgstarty->Fill(mcstart.Y());
		  mchists.fHgstartz->Fill(mcstart.Z());
		  mchists.fHgendx->Fill(mcend.X());
		  mchists.fHgendy->Fill(mcend.Y());
		  mchists.fHgendz->Fill(mcend.Z());
		  mchists.fHgtheta->Fill(part->Momentum().Theta());
		  mchists.fHgphi->Fill(part->Momentum().Phi());
		  mchists.fHgtheta_xz->Fill(mctheta_xz);
		  mchists.fHgtheta_yz->Fill(mctheta_yz);
		  mchists.fHgmom->Fill(part->Momentum().Vect().Mag());
		  mchists.fHglen->Fill(plen);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  void TrackAna::endJob()
  //
  // Purpose: End of job.
  //
  {
    // Print summary.

    mf::LogInfo("TrackAna") 
      << "TrackAna statistics:\n"
      << "  Number of events = " << fNumEvent;

    // Fill efficiency histograms.

    for(std::map<int, MCHists>::const_iterator i = fMCHistMap.begin();
	i != fMCHistMap.end(); ++i) {
      const MCHists& mchists = i->second;
      effcalc(mchists.fHgstartx, mchists.fHmcstartx, mchists.fHestartx);
      effcalc(mchists.fHgstarty, mchists.fHmcstarty, mchists.fHestarty);
      effcalc(mchists.fHgstartz, mchists.fHmcstartz, mchists.fHestartz);
      effcalc(mchists.fHgendx, mchists.fHmcendx, mchists.fHeendx);
      effcalc(mchists.fHgendy, mchists.fHmcendy, mchists.fHeendy);
      effcalc(mchists.fHgendz, mchists.fHmcendz, mchists.fHeendz);
      effcalc(mchists.fHgtheta, mchists.fHmctheta, mchists.fHetheta);
      effcalc(mchists.fHgphi, mchists.fHmcphi, mchists.fHephi);
      effcalc(mchists.fHgtheta_xz, mchists.fHmctheta_xz, mchists.fHetheta_xz);
      effcalc(mchists.fHgtheta_yz, mchists.fHmctheta_yz, mchists.fHetheta_yz);
      effcalc(mchists.fHgmom, mchists.fHmcmom, mchists.fHemom);
      effcalc(mchists.fHglen, mchists.fHmclen, mchists.fHelen);
    }
  }
}
