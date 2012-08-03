//
// Name: TrackAna.cxx
//
// Purpose: Implementation file for module TrackAna.
//
// Created: 26-Jun-2012  H. Greenlee
//

#include <iostream>
#include <cmath>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Geometry/geo.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "TrackFinder/TrackAna.h"
#include "art/Framework/Principal/Event.h"
#include "RecoBase/Track.h"
#include "Simulation/SimListUtils.h"

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
}

namespace trkf {

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
    fBooked(false),
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
    fHmom(0),
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
    fHpvsp(0),
    fHpvspc(0),
    fHdp(0),
    fHppull(0),
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

  void TrackAna::bookHistograms(bool mc)
  //
  // Purpose: Book histograms.
  //
  {
    if(!fBooked) {
      fBooked = true;

      // Get services.

      art::ServiceHandle<geo::Geometry> geom;
      art::ServiceHandle<art::TFileService> tfs;

      // Make histogram directory.

      art::TFileDirectory dir = tfs->mkdir("trkana", "TrackAna histograms");

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
      fHmom = dir.make<TH1F>("mom", "Momentum", 100, 0., 10.);
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
      fHpvsp = dir.make<TH2F>("pvsp", "Reco Momentum vs. MC Truth Momentum",
			      100, 0., 5., 100, 0., 5.);
      fHpvspc = dir.make<TH2F>("pvspc", "Reco Momentum vs. MC Truth Momentum (Contained Tracks)",
			       100, 0., 5., 100, 0., 5.);
      fHdp = dir.make<TH1F>("dp", "Reco-MC Momentum Difference", 100, -5., 5.);
      fHdpc = dir.make<TH1F>("dpc", "Reco-MC Momentum Difference (Contained Tracks)",
			     100, -5., 5.);
      fHppull = dir.make<TH1F>("ppull", "Momentum Pull", 100, -10., 10.);
      fHppullc = dir.make<TH1F>("ppullc", "Momentum Pull (Contained Tracks)", 100, -10., 10.);
    }
  }

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
    bookHistograms(mc);

    // Get mc particles.

    sim::ParticleList plist;
    if(mc)
      plist = sim::SimListUtils::GetParticleList(evt, fG4ModuleLabel);

    // Get tracks.

    art::Handle< std::vector<recob::Track> > trackh;
    evt.getByLabel(fTrackModuleLabel, trackh);
    if(trackh.isValid()) {

      // Loop over tracks.

      int ntrack = trackh->size();
      for(int i = 0; i < ntrack; ++i) {
	art::Ptr<recob::Track> ptrack(trackh, i);
	const recob::Track& track = *ptrack;

	// Fill histograms involving tracks only.

	int ntraj = track.NumberTrajectoryPoints();
	if(ntraj > 0) {
	  const TVector3& pos = track.Vertex();
	  const TVector3& dir = track.VertexDirection();
	  const TVector3& end = track.End();

	  double dpos = bdist(pos);
	  double dend = bdist(end);

	  fHstartx->Fill(pos.X());
	  fHstarty->Fill(pos.Y());
	  fHstartz->Fill(pos.Z());
	  fHstartd->Fill(dpos);
	  fHendx->Fill(end.X());
	  fHendy->Fill(end.Y());
	  fHendz->Fill(end.Z());
	  fHendd->Fill(dend);
	  fHtheta->Fill(dir.Theta());
	  fHphi->Fill(dir.Phi());

	  double mom = 0.;
	  if(track.NumberFitMomentum() > 0)
	    mom = track.VertexMomentum();
	  fHmom->Fill(mom);

	  // Calculate the global-to-local rotation matrix.

	  TMatrixD rot(3,3);
	  track.GlobalToLocalRotationAtPoint(0, rot);

	  // Get covariance matrix.

	  const TMatrixD& cov = track.VertexCovariance();

	  // Loop over mc particls.

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

		fHduvcosth->Fill(colinearity, uv0);
		if(std::abs(uv0) < fMatchDisp) {

		  // Fill slope matching histograms.

		  fHmcdudw->Fill(dudw);
		  fHmcdvdw->Fill(dvdw);
		  fHdudwpull->Fill(dudw / std::sqrt(cov(2,2)));
		  fHdvdwpull->Fill(dvdw / std::sqrt(cov(3,3)));
		}
		fHcosth->Fill(colinearity);
		if(colinearity > fMatchColinearity) {

		  // Fill displacement matching histograms.

		  fHmcu->Fill(u0);
		  fHmcv->Fill(v0);
		  fHmcw->Fill(w);
		  fHupull->Fill(u0 / std::sqrt(cov(0,0)));
		  fHvpull->Fill(v0 / std::sqrt(cov(1,1)));

		  if(std::abs(uv0) < fMatchDisp) {

		    // Fill matching histograms.

		    fHpvsp->Fill(mcmom.Mag(), mom);
		    double dp = mom - mcmom.Mag();
		    fHdp->Fill(dp);
		    fHppull->Fill(dp / std::sqrt(cov(4,4)));
		    if(std::abs(dpos) >= 5. && std::abs(dend) >= 5.) {
		      fHpvspc->Fill(mcmom.Mag(), mom);
		      fHdpc->Fill(dp);
		      fHppullc->Fill(dp / std::sqrt(cov(4,4)));
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}
