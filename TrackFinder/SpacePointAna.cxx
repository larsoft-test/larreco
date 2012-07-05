//
// Name: SpacePointAna.cxx
//
// Purpose: Implementation file for module SpacePointAna.
//
// Created: 2-Aug-2011  H. Greenlee
//

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

#include "TrackFinder/SpacePointAna.h"
#include "Geometry/geo.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RecoBase/recobase.h"
#include "MCCheater/BackTracker.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 

namespace trkf {

  SpacePointAna::SpacePointAna(const fhicl::ParameterSet& pset) :
    //
    // Purpose: Constructor.
    //
    // Arguments: pset - Module parameters.
    //
    fSptalg1(pset.get<fhicl::ParameterSet>("SpacePointAlg")),
    fSptalg2(pset.get<fhicl::ParameterSet>("SpacePointAlg")),
    fSptalg3(pset.get<fhicl::ParameterSet>("SpacePointAlg")),
    fHitModuleLabel(pset.get<std::string>("HitModuleLabel")),
    fUseClusterHits(pset.get<bool>("UseClusterHits")),
    fClusterModuleLabel(pset.get<std::string>("ClusterModuleLabel")),
    fUseMC(pset.get<bool>("UseMC")),
    fG4ModuleLabel(pset.get<std::string>("G4ModuleLabel")),
    fMaxDT(pset.get<double>("MaxDT")),
    fMaxS(pset.get<double>("MaxS")),
    fBooked(false),
    fHDTUE(0),
    fHDTVE(0),
    fHDTWE(0),
    fHDTUPull(0),
    fHDTVPull(0),
    fHDTWPull(0),
    fHDTUV(0),
    fHDTVW(0),
    fHDTWU(0),
    fHS(0),
    fHchisq(0),
    fHx(0),
    fHy(0),
    fHz(0),
    fHMCdx(0),
    fHMCdy(0),
    fHMCdz(0),
    fHMCxpull(0),
    fHMCypull(0),
    fHMCzpull(0),
    fNumEvent(0)
  {

    // Report.

    mf::LogInfo("SpacePointAna") 
      << "SpacePointAna configured with the following parameters:\n"
      << "  HitModuleLabel = " << fHitModuleLabel << "\n"
      << "  UseClusterHits = " << fUseClusterHits << "\n"
      << "  ClusterModuleLabel = " << fClusterModuleLabel << "\n"
      << "  UseMC = " << fUseMC << "\n"
      << "  G4ModuleLabel = " << fG4ModuleLabel << "\n"
      << "  MaxDT = " << fMaxDT << "\n"
      << "  MaxS = " << fMaxS;
  }

  SpacePointAna::~SpacePointAna()
  //
  // Purpose: Destructor.
  //
  {}

  void SpacePointAna::bookHistograms(bool mc)
  //
  // Purpose: Book histograms.
  //
  {
    if(!fBooked) {
      fBooked = true;

      art::ServiceHandle<geo::Geometry> geom;
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory dir = tfs->mkdir("sptana", "SpacePointAna histograms");

      if(mc) {
	fHDTUE = dir.make<TH1F>("MCDTUE", "U-Drift Electrons Time Difference", 100, -2., 2.);
	fHDTVE = dir.make<TH1F>("MCDTVE", "V-Drift Electrons Time Difference", 100, -2., 2.);
	fHDTWE = dir.make<TH1F>("MCDTWE", "W-Drift Electrons Time Difference", 100, -2., 2.);
	fHDTUPull = dir.make<TH1F>("MCDTUPull", "U-Drift Electrons Time Pull", 100, -50., 50.);
	fHDTVPull = dir.make<TH1F>("MCDTVPull", "V-Drift Electrons Time Pull", 100, -50., 50.);
	fHDTWPull = dir.make<TH1F>("MCDTWPull", "W-Drift Electrons Time Pull", 100, -50., 50.);
      }
      if(!fSptalg1.merge()) {
	fHDTUV = dir.make<TH1F>("DTUV", "U-V time difference", 100, -20., 20.);
	fHDTVW = dir.make<TH1F>("DTVW", "V-W time difference", 100, -20., 20.);
	fHDTWU = dir.make<TH1F>("DTWU", "W-U time difference", 100, -20., 20.);
	fHS = dir.make<TH1F>("DS", "Spatial Separatoin", 100, -2., 2.);
      }
      fHchisq = dir.make<TH1F>("chisq", "Chisquare", 100, 0., 20.);

      fHx = dir.make<TH1F>("xpos", "X Position",
			   100, 0., 2.*geom->DetHalfWidth());
      fHy = dir.make<TH1F>("ypos", "Y Position",
			   100, -geom->DetHalfHeight(), geom->DetHalfHeight());
      fHz = dir.make<TH1F>("zpos", "Z Position",
			   100, 0., geom->DetLength());
      if(mc) {
	fHMCdx = dir.make<TH1F>("MCdx", "X MC Residual", 100, -1., 1.);
	fHMCdy = dir.make<TH1F>("MCdy", "Y MC Residual", 100, -1., 1.);
	fHMCdz = dir.make<TH1F>("MCdz", "Z MC Residual", 100, -1., 1.);
	fHMCxpull = dir.make<TH1F>("MCxpull", "X MC Pull", 100, -50., 50.);
	fHMCypull = dir.make<TH1F>("MCypull", "Y MC Pull", 100, -50., 50.);
	fHMCzpull = dir.make<TH1F>("MCzpull", "Z MC Pull", 100, -50., 50.);
      }
    }
  }

  void SpacePointAna::analyze(const art::Event& evt)
  //
  // Purpose: Analyze method.
  //
  // Arguments: event - Art event.
  //
  {
    art::ServiceHandle<cheat::BackTracker> bt;

    ++fNumEvent;

    // Make sure histograms are booked.

    bool mc = !evt.isRealData();
    bookHistograms(mc);

    // Get Services.

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detprop;
    art::ServiceHandle<util::LArProperties> larprop;

    // Get Hits.

    art::PtrVector<recob::Hit> hits;

    if(fUseClusterHits) {

      // Get clusters.

      art::Handle< std::vector<recob::Cluster> > clusterh;
      evt.getByLabel(fClusterModuleLabel, clusterh);

      // Get hits from all clusters.
      art::FindManyP<recob::Hit> fm(clusterh, evt, fClusterModuleLabel);

      if(clusterh.isValid()) {
	int nclus = clusterh->size();

	for(int i = 0; i < nclus; ++i) {
	  art::Ptr<recob::Cluster> pclus(clusterh, i);
	  std::vector< art::Ptr<recob::Hit> > clushits = fm.at(i);
	  int nhits = clushits.size();
	  hits.reserve(hits.size() + nhits);

	  for(std::vector< art::Ptr<recob::Hit> >::const_iterator ihit = clushits.begin();
	      ihit != clushits.end(); ++ihit) {
	    hits.push_back(*ihit);
	  }
	}
      }
    }
    else {

      // Get unclustered hits.

      art::Handle< std::vector<recob::Hit> > hith;
      evt.getByLabel(fHitModuleLabel, hith);
      if(hith.isValid()) {
	int nhits = hith->size();
	hits.reserve(nhits);

	for(int i = 0; i < nhits; ++i)
	  hits.push_back(art::Ptr<recob::Hit>(hith, i));
      }
    }

    // Fill histograms that don't depend on space points.

    if(mc) {

      // Loop over hits and fill hit-electron time difference histogram.

      for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
	  ihit != hits.end(); ++ihit) {
	const recob::Hit& hit = **ihit;

	unsigned short channel = hit.Channel();
	unsigned int tpc, plane, wire, cstat;
	geom->ChannelToWire(channel, cstat, tpc, plane, wire);
	geo::View_t view = hit.View();
	assert(geom->TPC(tpc).Plane(plane).View() == view);
	double tpeak = hit.PeakTime();
	//double tstart = hit.StartTime();
	//double tend = hit.EndTime();
	double terr = hit.SigmaPeakTime();

	assert(channel == hit.Channel());

	// Loop over electrons associated with this hit/channel and fill
	// hit-electron time difference histograms.

	// loop over the map of TDC to sim::IDE to get the TDC for each energy dep
	// Find the average time in ticks for this hit.

	//double sumw = 0.;
	//double sumt = 0.;

	//const std::map<unsigned short, std::vector<sim::IDE> > &idemap = simchan.TDCIDEMap();
	//std::map<unsigned short, std::vector<sim::IDE> >::const_iterator mitr = idemap.begin();
	//for(mitr = idemap.begin(); mitr != idemap.end(); mitr++) {
	//  double tdc = double((*mitr).first);
	//  if(tdc >= tstart && tdc <= tend) {
	//    const std::vector<sim::IDE>& idevec = (*mitr).second;
	//    for(std::vector<sim::IDE>::const_iterator iide=idevec.begin();
	//	iide != idevec.end(); ++iide) {
	//      const sim::IDE& ide = *iide;
	//      double w = ide.numElectrons;
	//      sumw += w;
	//      sumt += w*tdc;
	//    }
	//  }
	//}
	//double tav = 0.;
	//if(sumw != 0.)
	//  tav = sumt / sumw;

	std::vector<double> hitxyz = bt->HitToXYZ(*ihit);
	double tav = detprop->ConvertXToTicks(hitxyz[0], plane, tpc, cstat);

	if(view == geo::kU) {
	  fHDTUE->Fill(tpeak - tav);
	  fHDTUPull->Fill((tpeak - tav) / terr);
	}
	else if(view == geo::kV) {
	  fHDTVE->Fill(tpeak - tav);
	  fHDTVPull->Fill((tpeak - tav) / terr);
	}
	else if(view == geo::kW) {
	  fHDTWE->Fill(tpeak - tav);
	  fHDTWPull->Fill((tpeak - tav) / terr);
	}
	else
	  throw cet::exception("SpacePointAna") << "Bad view = " << view << "\n";
      }
    }
    
    std::vector<recob::SpacePoint> spts1;  // For time histograms.
    std::vector<recob::SpacePoint> spts2;  // For separation histogram.
    std::vector<recob::SpacePoint> spts3;  // Default cuts.

    // If nonzero time cut is specified, make space points using that
    // time cut (for time histograms).

    if(fMaxDT != 0 && !fSptalg1.merge()) {
      if(mc && fUseMC)
	fSptalg1.makeMCTruthSpacePoints(hits, spts1, 
					fSptalg1.filter(), fSptalg1.merge(),
					fMaxDT, 0.);
      else
	fSptalg1.makeSpacePoints(hits, spts1,
				 fSptalg1.filter(), fSptalg1.merge(),
				 fMaxDT, 0.);

      // Report number of space points.

      mf::LogDebug("SpacePointAna") << "Found " << spts1.size() 
				    << " space points using special time cut.";
    }

    // If nonzero separation cut is specified, make space points using that 
    // separation cut (for separation histogram).

    if(fMaxS != 0. && !fSptalg2.merge()) {
      if(mc && fUseMC)
	fSptalg2.makeMCTruthSpacePoints(hits, spts2, 
					fSptalg2.filter(), fSptalg2.merge(),
					0., fMaxS);
      else
	fSptalg2.makeSpacePoints(hits, spts2,
				 fSptalg2.filter(), fSptalg2.merge(),
				 0., fMaxS);

      // Report number of space points.

      mf::LogDebug("SpacePointAna") << "Found " << spts2.size() 
				    << " space points using special seperation cut.";
    }

    // Make space points using default cuts.

    if(mc && fUseMC)
      fSptalg3.makeMCTruthSpacePoints(hits, spts3);
    else
      fSptalg3.makeSpacePoints(hits, spts3);

    // Report number of space points.

    mf::LogDebug("SpacePointAna") << "Found " << spts3.size() 
				  << " space points using default cuts.";

    std::vector<recob::SpacePoint>::const_iterator ibegin;
    std::vector<recob::SpacePoint>::const_iterator iend;

    if(!fSptalg1.merge()) {

      // Loop over space points and fill time histograms.

      ibegin = (fMaxDT != 0. ? spts1.begin() : spts3.begin());
      iend = (fMaxDT != 0. ? spts1.end() : spts3.end());
      const SpacePointAlg* palg = (fMaxDT != 0. ? &fSptalg1 : &fSptalg3);

      for(std::vector<recob::SpacePoint>::const_iterator i = ibegin; 
	  i != iend; ++i) {
	const recob::SpacePoint& spt = *i;

	// Get hits associated with this SpacePoint.

	const art::PtrVector<recob::Hit>& spthits = palg->getAssociatedHits(spt);

	// Make a double loop over hits and fill hit time difference histograms.

	for(art::PtrVector<recob::Hit>::const_iterator ihit = spthits.begin();
	    ihit != spthits.end(); ++ihit) {
	  const recob::Hit& hit1 = **ihit;

	  unsigned short channel1 = hit1.Channel();
	  unsigned int tpc1, plane1, wire1, cs1;
	  geom->ChannelToWire(channel1, cs1, tpc1, plane1, wire1);
	  geo::View_t view1 = hit1.View();
	  double t1 = palg->correctedTime(hit1);

	  for(art::PtrVector<recob::Hit>::const_iterator jhit = spthits.begin();
	      jhit != spthits.end(); ++jhit) {
	    const recob::Hit& hit2 = **jhit;

	    unsigned short channel2 = hit2.Channel();
	    unsigned int tpc2, plane2, wire2, cs2;
	    geom->ChannelToWire(channel2, cs2, tpc2, plane2, wire2);

	    // Require same tpc, different view.

	    if(tpc1 == tpc2 && plane1 != plane2) {

	      geo::View_t view2 = hit2.View();
	      double t2 = palg->correctedTime(hit2);

	      if(view1 == geo::kU) {
		if(view2 == geo::kV)
		  fHDTUV->Fill(t1-t2);
		if(view2 == geo::kW)
		  fHDTWU->Fill(t2-t1);
	      }
	      if(view1 == geo::kV) {
		if(view2 == geo::kW)
		  fHDTVW->Fill(t1-t2);
		if(view2 == geo::kU)
		  fHDTUV->Fill(t2-t1);
	      }
	      if(view1 == geo::kW) {
		if(view2 == geo::kU)
		  fHDTWU->Fill(t1-t2);
		if(view2 == geo::kV)
		  fHDTVW->Fill(t2-t1);
	      }
	    }
	  }
	}
      }

      // Loop over space points and fill seperation histograms.

      ibegin = (fMaxS != 0. ? spts2.begin() : spts3.begin());
      iend = (fMaxS != 0. ? spts2.end() : spts3.end());
      palg = (fMaxDT != 0. ? &fSptalg2 : &fSptalg3);

      for(std::vector<recob::SpacePoint>::const_iterator i = ibegin; 
	  i != iend; ++i) {
	const recob::SpacePoint& spt = *i;

	// Get hits associated with this SpacePoint.

	const art::PtrVector<recob::Hit>& spthits = palg->getAssociatedHits(spt);

	// Fill separation histogram.

	double sep = palg->separation(spthits);
	fHS->Fill(sep);
      }
    }

    // Loop over default space points and fill histograms.

    ibegin = spts3.begin();
    iend = spts3.end();

    for(std::vector<recob::SpacePoint>::const_iterator i = ibegin;
	i != iend; ++i) {
      const recob::SpacePoint& spt = *i;

      fHchisq->Fill(spt.Chisq());
      fHx->Fill(spt.XYZ()[0]);
      fHy->Fill(spt.XYZ()[1]);
      fHz->Fill(spt.XYZ()[2]);
      if(mc) {
	std::vector<double> mcxyz = bt->SpacePointHitsToXYZ(fSptalg3.getAssociatedHits(spt));
	fHMCdx->Fill(spt.XYZ()[0] - mcxyz[0]);
	fHMCdy->Fill(spt.XYZ()[1] - mcxyz[1]);
	fHMCdz->Fill(spt.XYZ()[2] - mcxyz[2]);
	if(spt.ErrXYZ()[0] > 0.)
	  fHMCxpull->Fill((spt.XYZ()[0] - mcxyz[0]) / std::sqrt(spt.ErrXYZ()[0]));
	if(spt.ErrXYZ()[2] > 0.)
	  fHMCypull->Fill((spt.XYZ()[1] - mcxyz[1]) / std::sqrt(spt.ErrXYZ()[2]));
	if(spt.ErrXYZ()[5] > 0.)
	  fHMCzpull->Fill((spt.XYZ()[2] - mcxyz[2]) / std::sqrt(spt.ErrXYZ()[5]));
      }
    }
  }
}
