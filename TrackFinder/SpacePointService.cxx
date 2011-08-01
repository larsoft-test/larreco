////////////////////////////////////////////////////////////////////////
///
/// \file   SpacePointService.cxx
///
/// \brief  Service for generating space points from hits.
///
/// \author H. Greenlee 
///
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <cmath>
#include <map>
#include <algorithm>
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TrackFinder/SpacePointService.h"
#include "Geometry/geo.h"
#include "Utilities/LArProperties.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "Simulation/SimChannel.h"
#include "Simulation/Electrons.h"
#include "Simulation/LArVoxelData.h"
#include "art/Framework/Core/View.h"
#include "Utilities/DetectorProperties.h"
#include "TH1F.h"

namespace {

  // Function classes for sorting pointer to Electrons according to
  // embedded position TVector3.

  class ElectronsLess {
  public:
    bool operator()(const sim::Electrons* p1, const sim::Electrons* p2) {
      TVector3 v1 = p1->XYZ();
      TVector3 v2 = p2->XYZ();
      bool result = (v1.x() < v2.x() || 
		     (v1.x() == v2.x() && 
		      (v1.y() < v2.y() || 
		       (v1.y() == v2.y() && v1.z() < v2.z()))));
      return result;
    }
  };

  class ElectronsEqual {
  public:
    bool operator()(const sim::Electrons* p1, const sim::Electrons* p2) {
      TVector3 v1 = p1->XYZ();
      TVector3 v2 = p2->XYZ();
      bool result = (v1.x() == v2.x() && v1.y() == v2.y() && v1.z() < v2.z());
      return result;
    }
  };

}

//----------------------------------------------------------------------
// Constructor.
//
trkf::SpacePointService::SpacePointService(const fhicl::ParameterSet& pset,
					   art::ActivityRegistry& reg) :
  fHist(false),
  fMCHist(false),
  fUseMC(false),
  fMaxDT(0.),
  fMaxS(0.),
  fTimeOffsetU(0.),
  fTimeOffsetV(0.),
  fTimeOffsetW(0.),
  fMinViews(1000),
  fEnableU(false),
  fEnableV(false),
  fEnableW(false),
  fGeom(0),
  fSamplingRate(0.),
  fTriggerOffset(0),
  fEfield(0.),
  fTemperature(0.),
  fDriftVelocity(0.),
  fTimePitch(0.),
  fBooked(false),
  fHDTUV(0),
  fHDTVW(0),
  fHDTWU(0),
  fHS(0),
  fHDTUVC(0),
  fHDTVWC(0),
  fHDTWUC(0),
  fHx(0),
  fHy(0),
  fHz(0),
  fHMCdx(0),
  fHMCdy(0),
  fHMCdz(0)
{
  reg.watchPostBeginRun(this, &SpacePointService::postBeginRun);
  reconfigure(pset);
}

//----------------------------------------------------------------------
// Destructor.
//
trkf::SpacePointService::~SpacePointService()
{
}

//----------------------------------------------------------------------
// Update configuration parameters.
//
void trkf::SpacePointService::reconfigure(const fhicl::ParameterSet& pset)
{
  // Get configuration parameters.

  fHist = pset.get<bool>("Hist", false);
  fMCHist = pset.get<bool>("MCHist", false);
  fUseMC = pset.get<bool>("UseMC", false);
  fMaxDT = pset.get<double>("MaxDT", 0.);
  fMaxS = pset.get<double>("MaxS", 0.);

  fTimeOffsetU = pset.get<double>("TimeOffsetU", 0.);
  fTimeOffsetV = pset.get<double>("TimeOffsetV", 0.);
  fTimeOffsetW = pset.get<double>("TimeOffsetW", 0.);

  fMinViews = pset.get<int>("MinViews", 1000);

  fEnableU = pset.get<bool>("EnableU", false);
  fEnableV = pset.get<bool>("EnableV", false);
  fEnableW = pset.get<bool>("EnableW", false);

  fGDMLPath = "";
  fROOTPath = "";

  // Report.

  mf::LogInfo("SpacePointService") 
    << "SpacePointService configured with the following parameters:\n"
    << "  Hist = " << fHist << "\n"
    << "  MC Hist = " << fMCHist << "\n"
    << "  UseMC = " << fUseMC << "\n"
    << "  MaxDT = " << fMaxDT << "\n"
    << "  MaxS = " << fMaxS << "\n"
    << "  TimeOffsetU = " << fTimeOffsetU << "\n"
    << "  TimeOffsetV = " << fTimeOffsetV << "\n"
    << "  TimeOffsetW = " << fTimeOffsetW << "\n" 
    << "  MinViews = " << fMinViews << "\n"
    << "  EnableU = " << fEnableU << "\n"
    << "  EnableV = " << fEnableV << "\n"
    << "  EnableW = " << fEnableW;
}

//----------------------------------------------------------------------
// ART signal.  We check here if geometry needs to be updated.
// The Geometry service also updates in the same signal, so this service
// should come after geometry in service list.  There should be a way
// to register with the geometry service (as opposed to the framework)
// to be notified if the geometry has been updated...
//
void trkf::SpacePointService::postBeginRun(art::Run const& run)
{
  maybe_update();
}

//----------------------------------------------------------------------
// This method updates geometry constants, if needed.
//
void trkf::SpacePointService::maybe_update()
{
  art::ServiceHandle<geo::Geometry> geom;

  std::string gdml_path = geom->GetGDMLPath();
  std::string root_path = geom->GetROOTPath();

  if(fGeom == 0 || gdml_path != fGDMLPath || root_path != fROOTPath) {
    fGeom = &*geom;
    fGDMLPath = gdml_path;
    fROOTPath = root_path;
    update();
  }
}

//----------------------------------------------------------------------
// Update geometry constants.
//
void trkf::SpacePointService::update()
{
  // Reset geometry constants so that we don't
  // accidentally inherit from previous update.

  fEnable.clear();
  fTimeOffset.clear();
  fWirePitch.clear();
  fWireOffset.clear();
  fTheta.clear();
  fSinTheta.clear();
  fCosTheta.clear();
  fSin.clear();

  // Calculate and print geometry information.

  mf::LogInfo log("SpacePointService");
  log << "Updating geometry constants.\n";

  // Loop over TPCs.

  int ntpc = fGeom->NTPC();

  fEnable.resize(ntpc);
  fTimeOffset.resize(ntpc);
  fWirePitch.resize(ntpc);
  fWireOffset.resize(ntpc);
  fTheta.resize(ntpc);
  fSinTheta.resize(ntpc);
  fCosTheta.resize(ntpc);
  fSin.resize(ntpc);

  for(int tpc = 0; tpc < ntpc; ++tpc) {
    const geo::TPCGeo& tpcgeom = fGeom->TPC(tpc);

    // Loop over planes.

    int nplane = tpcgeom.Nplanes();

    fEnable[tpc].resize(nplane, false);
    fTimeOffset[tpc].resize(nplane, 0.);
    fWirePitch[tpc].resize(nplane, 0.);
    fWireOffset[tpc].resize(nplane, 0.);
    fTheta[tpc].resize(nplane, 0.);
    fSinTheta[tpc].resize(nplane, 0.);
    fCosTheta[tpc].resize(nplane, 0.);
    fSin[tpc].resize(nplane, 0.);

    for(int plane = 0; plane < nplane; ++plane) {
      const geo::PlaneGeo& pgeom = fGeom->Plane(plane);
      double theta = pgeom.Wire(0).ThetaZ();

      // Calculate the perpendicular distance of the first and last wires
      // from the origin in the (y,z) plane.  Use this information to
      // calculate the wire pitch and offset.

      int wire[2] = {0, 0};
      wire[1] = pgeom.Nwires() - 1;
      double dist[2] = {0., 0.};
      for(int j=0; j<2; ++j) {
	double xyz1[3], xyz2[3];
	double z1 = -100;
	double z2 = 100.;
	const geo::WireGeo& wgeom = pgeom.Wire(wire[j]);
	wgeom.GetCenter(xyz1, z1);
	wgeom.GetCenter(xyz2, z2);
	dist[j] = (xyz2[1] * xyz1[2] - xyz1[1] * xyz2[2]) / (z2 - z1);
      }
      double pitch = (dist[1] - dist[0]) / (wire[1] - wire[0]);
      double offset = dist[0] - pitch * wire[0];

      geo::View_t view = pgeom.View();
      std::string viewname = "?";
      if(view == geo::kU) {
	fEnable[tpc][plane] = fEnableU;
	fTimeOffset[tpc][plane] = fTimeOffsetU;
	viewname = "U";
      }
      else if(view == geo::kV) {
	fEnable[tpc][plane] = fEnableV;
	fTimeOffset[tpc][plane] = fTimeOffsetV;
	viewname = "V";
      }
      else if(view == geo::kW) {
	fEnable[tpc][plane] = fEnableW;
	fTimeOffset[tpc][plane] = fTimeOffsetW;
	viewname = "W";
      }
      else
	throw cet::exception("SpacePointService") << "Bad view = " << view << "\n";

      fWirePitch[tpc][plane] = pitch;
      fWireOffset[tpc][plane] = offset;
      fTheta[tpc][plane] = theta;
      fSinTheta[tpc][plane] = std::sin(theta);
      fCosTheta[tpc][plane] = std::cos(theta);

      log << "\nTPC, Plane: " << tpc << ", " << plane << "\n"
	  << "  View " << view << " (" << viewname << ")\n"
	  << "  SignalType " << pgeom.SignalType() << "\n"
	  << "  Orientation " << pgeom.Orientation() << "\n"
	  << "  Theta = " << theta << "\n"
	  << "  Wire pitch = " << pitch << "cm\n"
	  << "  Wire offset = " << offset << "cm\n";
    }

    if(fTheta[tpc].size() == 3) {
      fSin[tpc][0] = std::sin(fTheta[tpc][1] - fTheta[tpc][2]);
      fSin[tpc][1] = std::sin(fTheta[tpc][2] - fTheta[tpc][0]);
      fSin[tpc][2] = std::sin(fTheta[tpc][0] - fTheta[tpc][1]);
    }

    log << "  sin(V-W) = " << fSin[tpc][0] << "\n"
	<< "  sin(W-U) = " << fSin[tpc][1] << "\n"
	<< "  sin(U-V) = " << fSin[tpc][2] << "\n";
  }


  // Update detector properties.
  art::ServiceHandle<util::DetectorProperties> detprop;

  fSamplingRate = detprop->SamplingRate();
  fTriggerOffset = detprop->TriggerOffset();
  log << "\nDetector properties:\n"
      << "  Sampling Rate = " << fSamplingRate << " ns/tick\n"
      << "  Trigger offset = " << fTriggerOffset << " ticks\n";

  // Update LArProperties.
  art::ServiceHandle<util::LArProperties> larprop;
  //  art::ServiceHandle<util::LArProperties> larprop;
  fEfield = larprop->Efield();
  fTemperature = larprop->Temperature();
  fDriftVelocity = larprop->DriftVelocity(fEfield, fTemperature);
  fTimePitch = 0.001 * fDriftVelocity * fSamplingRate;
  log << "\nLAr propertoes:\n"
      << "  E field = " << fEfield << " kV/cm\n"
      << "  Temperature = " << fTemperature << " K\n"
      << "  Drift velocity = " << fDriftVelocity << " cm/us\n"
      << "  Time pitch = " << fTimePitch << " cm/tick";

  // Make sure histograms have been booked.

  bookHistograms();
}

//----------------------------------------------------------------------
// Book histograms.
//
void trkf::SpacePointService::bookHistograms()
{
  if(fHist && !fBooked) {
    fBooked = true;

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir("sptsvc", "SpacePointService histograms");
    if(fMCHist) {
      fHDTUE = dir.make<TH1F>("MCDTUE", "U-Drift Electrons Time Difference", 100, -50., 50.);
      fHDTVE = dir.make<TH1F>("MCDTVE", "V-Drift Electrons Time Difference", 100, -50., 50.);
      fHDTWE = dir.make<TH1F>("MCDTWE", "W-Drift Electrons Time Difference", 100, -50., 50.);
    }
    fHDTUV = dir.make<TH1F>("DTUV", "U-V time difference", 100, -50., 50.);
    fHDTVW = dir.make<TH1F>("DTVW", "V-W time difference", 100, -50., 50.);
    fHDTWU = dir.make<TH1F>("DTWU", "W-U time difference", 100, -50., 50.);
    fHS = dir.make<TH1F>("DS", "Wire Coincidence", 100, -10., 10.);
    fHDTUVC = dir.make<TH1F>("DTUVC", "U-V time difference (space compatible hits)",
			      100, -50., 50.);
    fHDTVWC = dir.make<TH1F>("DTVWC", "V-W time difference (space compatible hits)",
			      100, -50., 50.);
    fHDTWUC = dir.make<TH1F>("DTWUC", "W-U time difference (space compatible hits)",
			      100, -50., 50.);
    fHx = dir.make<TH1F>("xpos", "X Position",
			 100, 0., 2.*fGeom->DetHalfWidth());
    fHy = dir.make<TH1F>("ypos", "Y Position",
			 100, -fGeom->DetHalfHeight(), fGeom->DetHalfHeight());
    fHz = dir.make<TH1F>("zpos", "Z Position",
			 100, 0., fGeom->DetLength());
    if(fMCHist) {
      fHMCdx = dir.make<TH1F>("MCdx", "X MC Residual", 100, -5., 5.);
      fHMCdy = dir.make<TH1F>("MCdy", "Y MC Residual", 100, -5., 5.);
      fHMCdz = dir.make<TH1F>("MCdz", "Z MC Residual", 100, -5., 5.);
    }
  }
}

//----------------------------------------------------------------------
// Check hits for compatibility.
// Check hits pairwise for different views and maximum time difference.
// Check three hits for spatial compatibility.
// This method fills histograms if histograms are enabled.
bool trkf::SpacePointService::compatible(const art::PtrVector<recob::Hit>& hits) const
{
  int nhits = hits.size();

  // Fewer than two or more than three hits can never be compatible.

  bool result = nhits >= 2 && nhits <= 3;
  bool mc_ok = true;
  unsigned int tpc = 0;

  if(result) {

    // First do pairwise tests.
    // Do double loop over hits.

    for(int ihit1 = 0; result && ihit1 < nhits-1; ++ihit1) {
      const recob::Hit& hit1 = *(hits[ihit1]);
      unsigned short channel1 = hit1.Channel();
      unsigned int tpc1, plane1, wire1;
      fGeom->ChannelToWire(channel1, tpc1, plane1, wire1);
      geo::View_t view1 = hit1.View();
      double t1 = hit1.PeakTime() - fTimeOffset[tpc1][plane1];

      // If using mc information, make a collection of electrons for hit 1.

      const std::vector<const sim::Electrons*>& electrons1 = fElecMap[&hit1];

      // Loop over second hit.

      for(int ihit2 = ihit1+1; result && ihit2 < nhits; ++ihit2) {
	const recob::Hit& hit2 = *(hits[ihit2]);
	unsigned short channel2 = hit2.Channel();
	unsigned int tpc2, plane2, wire2;
	fGeom->ChannelToWire(channel2, tpc2, plane2, wire2);
	geo::View_t view2 = hit2.View();

	// Test for same tpc and different views.

	result = result && tpc1 == tpc2 && view1 != view2;
	if(result) {

	  // Remember which tpc we are in.

	  tpc = tpc1;

	  double t2 = hit2.PeakTime() - fTimeOffset[tpc2][plane2];
    
	  // Test maximum time difference.

	  result = result && std::abs(t1-t2) <= fMaxDT;

	  // If using mc information, make a collection of electrons for hit 2.

	  if(fUseMC) {

	    // Find electrons that are in time with hit 2.

	    std::vector<const sim::Electrons*> electrons2 = fElecMap[&hit2];
	    std::vector<const sim::Electrons*>::iterator it =
	      std::set_intersection(electrons1.begin(), electrons1.end(),
				    electrons2.begin(), electrons2.end(),
				    electrons2.begin(), ElectronsLess());
	    electrons2.resize(it - electrons2.begin());

	    // Hits are compatible if they have points in common.

	    mc_ok = electrons2.size() > 0;
	    result = result && mc_ok;
	  }

	  // Fill histograms.

	  if(fHist && mc_ok) {
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

    // If there are exactly three hits, and they pass pairwise tests, check
    // for spatial compatibility.

    if(result && nhits == 3) {

      // Loop over hits.

      double dist[3] = {0., 0., 0.};
      double time[3];
      geo::View_t view[3];

      for(int i=0; i<3; ++i) {

	// Get tpc, plane, wire.

	const recob::Hit& hit = *(hits[i]);
	unsigned short channel = hit.Channel();
	unsigned int tpc0, plane, wire;
	fGeom->ChannelToWire(channel, tpc0, plane, wire);
	assert(tpc0 == tpc);
	view[i] = hit.View();

	// Get corrected time.

	time[i] = hit.PeakTime() - fTimeOffset[tpc][plane];

	// Get distance with offset correction.

	dist[plane] = fWirePitch[tpc][plane] * wire + fWireOffset[tpc][plane];
      }

      // Do space cut.

      double S = fSin[tpc][0]*dist[0] + fSin[tpc][1]*dist[1] + fSin[tpc][2]*dist[2];
      result = result && std::abs(S) < fMaxS;

      // Fill histograms.

      if(fHist) {
	fHS->Fill(S);
	if(result) {
	  if(view[0] == geo::kU) {
	    if(view[1] == geo::kV) {
	      assert(view[2] == geo::kW);
	      fHDTUVC->Fill(time[0]-time[1]);
	      fHDTVWC->Fill(time[1]-time[2]);
	      fHDTWUC->Fill(time[2]-time[0]);
	    }
	    if(view[1] == geo::kW) {
	      assert(view[2] == geo::kV);
	      fHDTUVC->Fill(time[0]-time[2]);
	      fHDTVWC->Fill(time[2]-time[1]);
	      fHDTWUC->Fill(time[1]-time[0]);
	    }
	  }
	  if(view[0] == geo::kV) {
	    if(view[1] == geo::kW) {
	      assert(view[2] == geo::kU);
	      fHDTUVC->Fill(time[2]-time[0]);
	      fHDTVWC->Fill(time[0]-time[1]);
	      fHDTWUC->Fill(time[1]-time[2]);
	    }
	    if(view[1] == geo::kU) {
	      assert(view[2] == geo::kW);
	      fHDTUVC->Fill(time[1]-time[0]);
	      fHDTVWC->Fill(time[0]-time[2]);
	      fHDTWUC->Fill(time[2]-time[1]);
	    }
	  }
	  if(view[0] == geo::kW) {
	    if(view[1] == geo::kU) {
	      assert(view[2] == geo::kV);
	      fHDTUVC->Fill(time[1]-time[2]);
	      fHDTVWC->Fill(time[2]-time[0]);
	      fHDTWUC->Fill(time[0]-time[1]);
	    }
	    if(view[1] == geo::kV) {
	      assert(view[2] == geo::kU);
	      fHDTUVC->Fill(time[2]-time[1]);
	      fHDTVWC->Fill(time[1]-time[0]);
	      fHDTWUC->Fill(time[0]-time[2]);
	    }
	  }
	}
      }
    }
  }

  // Done.
    
  return result;
}

//----------------------------------------------------------------------
// Fill one space point using a colleciton of hits.
// Assume points have already been tested for compatibility.
//
void trkf::SpacePointService::fillSpacePoint(const art::PtrVector<recob::Hit>& hits,
					     recob::SpacePoint& spt) const
{
  spt = recob::SpacePoint(recob::SpacePoint(hits));
  int nhits = hits.size();

  // Calculate mc position as the average of the mc position over all hits.

  double mcxyz[3] = {0., 0., 0.};
  if(fMCHist) {
    for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
	ihit != hits.end(); ++ihit) {

      const recob::Hit& hit = **ihit;

      // Each hit should have a precalculated mc position.

      std::map<const recob::Hit*, TVector3>::const_iterator imc =
	fMCPosMap.find(&hit);
      assert(imc != fMCPosMap.end());
      const TVector3& hvec = imc->second;
      mcxyz[0] += hvec.x();
      mcxyz[1] += hvec.y();
      mcxyz[2] += hvec.z();
    }
    if(nhits > 0) {
      mcxyz[0] /= nhits;
      mcxyz[1] /= nhits;
      mcxyz[2] /= nhits;
    }
  }

  // Calculate position.

  double xyz[3] = {0., 0., 0.};

  // Calculate x using drift times.
  // Loop over all hits and calculate the weighted average drift time.

  double sumtw = 0.;
  double sumw = 0.;

  for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
      ihit != hits.end(); ++ihit) {

    const recob::Hit& hit = **ihit;
    unsigned short channel = hit.Channel();
    unsigned int tpc, plane, wire;
    fGeom->ChannelToWire(channel, tpc, plane, wire);

    // Correct time for trigger offset and view-dependent time offsets.
    // Assume time error is proportional to (end time - start time).

    double t = hit.PeakTime() - fTriggerOffset - fTimeOffset[tpc][plane];
    double et = hit.EndTime() - hit.StartTime();
    double w = 1./(et*et);

    sumtw += w*t;
    sumw += w;
  }

  double drift_time = 0.;
  if(sumw != 0.)
    drift_time = sumtw / sumw;
  xyz[0] = drift_time * fTimePitch;

  // Calculate y, z using wires (need at least two hits).

  if(nhits >= 2) {

    // Calculate y and z by chisquare minimization of wire coordinates.

    double sus = 0.;   // sum u_i sin_th_i
    double suc = 0.;   // sum u_i cos_th_i
    double sc2 = 0.;   // sum cos2_th_i
    double ss2 = 0.;   // sum sin2_th_i
    double ssc = 0.;   // sum sin_th_i cos_th_i

    // Loop over points.

    for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
	ihit != hits.end(); ++ihit) {

      const recob::Hit& hit = **ihit;
      unsigned short channel = hit.Channel();
      unsigned int tpc, plane, wire;
      fGeom->ChannelToWire(channel, tpc, plane, wire);

      // Calculate wire coordinate in this view.
    
      double u = wire * fWirePitch[tpc][plane] + fWireOffset[tpc][plane];

      // Summations

      double s = fSinTheta[tpc][plane];
      double c = fCosTheta[tpc][plane];
      sus += u*s;
      suc += u*c;
      sc2 += c*c;
      ss2 += s*s;
      ssc += s*c;
    }

    // Calculate y,z

    double denom = sc2 * ss2 - ssc * ssc;
    if(denom != 0.) {
      xyz[1] = (-suc * ss2 + sus * ssc) / denom;
      xyz[2] = (sus * sc2 - suc * ssc) / denom;
    }

    // Set coordintates in space point.

    spt.SetXYZ(xyz);

    // Fill histograms.

    if(fHist) {
      fHx->Fill(xyz[0]);
      fHy->Fill(xyz[1]);
      fHz->Fill(xyz[2]);
      if(fMCHist) {
	fHMCdx->Fill(xyz[0] - mcxyz[0]);
	fHMCdy->Fill(xyz[1] - mcxyz[1]);
	fHMCdz->Fill(xyz[2] - mcxyz[2]);
      }
    }
  }
}

//----------------------------------------------------------------------
// Fill a vector of space points for all compatible combinations of hits
// from an input vector of hits (non-mc-truth version).
//
void trkf::SpacePointService::makeSpacePoints(const art::Handle< std::vector<recob::Hit> >& hith,
					      std::vector<recob::SpacePoint>& spts) const
{
  art::PtrVector<recob::Hit> hits;
  int nhits = hith->size();
  hits.reserve(nhits);

  for(int i = 0; i < nhits; ++i)
    hits.push_back(art::Ptr<recob::Hit>(hith, i));

  makeSpacePoints(hits, spts);
}

//----------------------------------------------------------------------
// Fill a vector of space points for all compatible combinations of hits
// from an input vector of hits (non-mc-truth version).
//
void trkf::SpacePointService::makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
					      std::vector<recob::SpacePoint>& spts) const
{
  std::vector<sim::SimChannel> empty;
  makeSpacePoints(hits, spts, empty);
}

//----------------------------------------------------------------------
// Fill a vector of space points for all compatible combinations of hits
// from an input vector of hits (mc truth version).
//
void trkf::SpacePointService::makeSpacePoints(const art::Handle< std::vector<recob::Hit> >& hith,
					      std::vector<recob::SpacePoint>& spts,
					      const std::vector<sim::SimChannel>& simchans) const
{
  art::PtrVector<recob::Hit> hits;
  int nhits = hith->size();
  hits.reserve(nhits);

  for(int i = 0; i < nhits; ++i)
    hits.push_back(art::Ptr<recob::Hit>(hith, i));

  makeSpacePoints(hits, spts, simchans);
}

//----------------------------------------------------------------------
// Fill a vector of space points for all compatible combinations of hits
// from an input vector of hits (mc truth version).
//
void trkf::SpacePointService::makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
					      std::vector<recob::SpacePoint>& spts,
					      const std::vector<sim::SimChannel>& simchans) const
{
  // First make result vector is empty.

  spts.erase(spts.begin(), spts.end());

  // Statistics.

  int n2 = 0;  // Number of two-hit space points.
  int n3 = 0;  // Number of three-hit space points.

  // If fUseMC pr fMCHist is true, verify that channels are sorted by
  // channel number.

  if(fUseMC || fMCHist) {

    unsigned int nsc = simchans.size();
    for(unsigned int isc = 0; isc < nsc; ++isc) {
      const sim::SimChannel& sc = simchans[isc];
      if(isc != sc.Channel())
	throw cet::exception("SpacePointService") << "MC channels not sorted.\n";
    }
  }

  // Sort hits into maps indexed by [tpc][plane][time].
  // If using mc information, also generate maps of electrons and mc 
  // position indexed by hit.

  std::vector<std::vector<std::map<double, art::Ptr<recob::Hit> > > > hitmap;
  fElecMap.clear();
  fMCPosMap.clear();

  int ntpc = fEnable.size();
  hitmap.resize(ntpc);
  for(int tpc = 0; tpc < ntpc; ++tpc) {
    int nplane = fEnable[tpc].size();
    hitmap[tpc].resize(nplane);
  }

  /*
  art::View<art::PtrVector<recob::Hit> >  lhits = hits; 
  for(art::PtrVector<recob::Hit>::const_iterator ihit = lhits::begin();
      ihit != lhits::end(); ++ihit) {
  */

  for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
      ihit != hits.end(); ++ihit) {
    const art::Ptr<recob::Hit>& phit = *ihit;
    unsigned short channel = phit->Channel();
    unsigned int tpc, plane, wire;
    fGeom->ChannelToWire(channel, tpc, plane, wire);

    if(fEnable[tpc][plane]) {
      double t = phit->PeakTime() - fTimeOffset[tpc][plane];
      hitmap[tpc][plane][t] = phit;
      const recob::Hit& hit = *phit;

      // Get Electrons.

      if(fUseMC || fMCHist) {
	std::vector<const sim::Electrons*>& electrons = fElecMap[&hit];   // Empty vector.
	HitToElectrons(hit, simchans, electrons);
	fMCPosMap[&hit] = ElectronsToXYZ(electrons);
	std::sort(electrons.begin(), electrons.end(), ElectronsLess());
	std::vector<const sim::Electrons*>::iterator it = 
	  std::unique(electrons.begin(), electrons.end(), ElectronsEqual());
	electrons.resize(it - electrons.begin());
      }
    }
  }

  mf::LogDebug debug("SpacePointService");
  debug << "Total hits = " << hits.size() << "\n\n";

  for(int tpc = 0; tpc < ntpc; ++tpc) {
    int nplane=hitmap[tpc].size();
    for(int plane = 0; plane < nplane; ++plane) {
      debug << "TPC, Plane: " << tpc << ", " << plane 
	    << ", hits = " << hitmap[tpc][plane].size() << "\n";
    }
  }

  // Loop over TPCs.

  for(int tpc = 0; tpc < ntpc; ++tpc) {

    // Sort maps in increasing order of number of hits.
    // This is so that we can do the outer loops over hits 
    // over the views with fewer hits.

    int index[3] = {0, 1, 2};
    for(int i=0; i<2; ++i) {
      for(int j=i+1; j<3; ++j) {
	if(hitmap[tpc][index[i]].size() > hitmap[tpc][index[j]].size()) {
	  int temp = index[i];
	  index[i] = index[j];
	  index[j] = temp;
	}
      }
    }

    // If two-view space points are allowed, make a double loop
    // over hits and produce space points for compatible hit-pairs.

    if(fMinViews <= 2) {

      // Loop over pairs of views.

      for(int i=0; i<2; ++i) {
	int plane1 = index[i];
	for(int j=i+1; j<3; ++j) {
	  int plane2 = index[j];

	  assert(hitmap[tpc][plane1].size() <= hitmap[tpc][plane2].size());

	  // Loop over pairs of hits.

	  art::PtrVector<recob::Hit> hitvec;
	  hitvec.reserve(2);

	  for(std::map<double, art::Ptr<recob::Hit> >::const_iterator
		ihit1 = hitmap[tpc][plane1].begin();
	      ihit1 != hitmap[tpc][plane1].end(); ++ihit1) {

	    const art::Ptr<recob::Hit>& phit1 = ihit1->second;

	    double t1 = phit1->PeakTime() - fTimeOffset[tpc][plane1];
	    double t2min = t1 - fMaxDT;
	    double t2max = t1 + fMaxDT;

	    for(std::map<double, art::Ptr<recob::Hit> >::const_iterator 
		  ihit2 = hitmap[tpc][plane2].lower_bound(t2min);
		ihit2 != hitmap[tpc][plane2].upper_bound(t2max); ++ihit2) {

	      const art::Ptr<recob::Hit>& phit2 = ihit2->second;

	      // Check current pair of hits for compatibility.
	      // By construction, hits should always have compatible views 
	      // and times, but may not have compatible mc information.
	      // Calling method has the side effect of filling
	      // histograms if enabled.

	      hitvec.clear();
	      hitvec.push_back(phit1);
	      hitvec.push_back(phit2);
	      bool ok = compatible(hitvec);
	      if(ok) {

		// Add a space point.

		++n2;
		spts.push_back(recob::SpacePoint());
		fillSpacePoint(hitvec, spts.back());
	      }
	    }
	  }
	}
      }
    }

    // If three-view space points are allowed, make a tripe loop
    // over hits and produce space points for compatible triplets.

    if(fMinViews <= 3) {

      // Loop over triplets of hits.

      art::PtrVector<recob::Hit> hitvec;
      hitvec.reserve(3);

      int plane1 = index[0];
      for(std::map<double, art::Ptr<recob::Hit> >::const_iterator
	    ihit1 = hitmap[tpc][plane1].begin();
	  ihit1 != hitmap[tpc][plane1].end(); ++ihit1) {

	const art::Ptr<recob::Hit>& phit1 = ihit1->second;

	double t1 = phit1->PeakTime() - fTimeOffset[tpc][plane1];
	double t2min = t1 - fMaxDT;
	double t2max = t1 + fMaxDT;

	int plane2 = index[1];
	for(std::map<double, art::Ptr<recob::Hit> >::const_iterator 
	      ihit2 = hitmap[tpc][plane2].lower_bound(t2min);
	    ihit2 != hitmap[tpc][plane2].upper_bound(t2max); ++ihit2) {

	  const art::Ptr<recob::Hit>& phit2 = ihit2->second;

	  // Test first two hits for compatibility before looping 
	  // over third hit.

	  hitvec.clear();
	  hitvec.push_back(phit1);
	  hitvec.push_back(phit2);
	  bool ok = compatible(hitvec);
	  if(ok) {

	    double t2 = phit2->PeakTime() - fTimeOffset[tpc][plane2];
	    double t3min = std::max(t1, t2) - fMaxDT;
	    double t3max = std::min(t1, t2) + fMaxDT;

	    int plane3 = index[2];
	    for(std::map<double, art::Ptr<recob::Hit> >::const_iterator 
		  ihit3 = hitmap[tpc][plane3].lower_bound(t3min);
		ihit3 != hitmap[tpc][plane3].upper_bound(t3max); ++ihit3) {

	      const art::Ptr<recob::Hit>& phit3 = ihit3->second;

	      // Test triplet for compatibility.

	      hitvec.clear();
	      hitvec.push_back(phit1);
	      hitvec.push_back(phit2);
	      hitvec.push_back(phit3);
	      bool ok = compatible(hitvec);

	      if(ok) {

		// Add a space point.

		++n3;
		spts.push_back(recob::SpacePoint());
		fillSpacePoint(hitvec, spts.back());
	      }
	    }
	  }
	}
      }
    }
  }

  debug << "\n2-hit space points = " << n2 << "\n"
	<< "3-hit space points = " << n3;
}

//----------------------------------------------------------------------
// Extract Electrons associated with Hit.
// (Similar as BackTracker::HitToXYZ in MCCheater.)
void trkf::SpacePointService::HitToElectrons(const recob::Hit& hit,
					     const std::vector<sim::SimChannel>& simchans,
					     std::vector<const sim::Electrons*>& electrons) const
{
  // Find SimChannel associated with Hit.

  unsigned short channel = hit.Channel();
  const sim::SimChannel& sc = simchans[channel];
  assert(sc.Channel() == channel);

  // Find electrons in time with Hit.

  double tstart = hit.StartTime();
  double tpeak = hit.PeakTime();
  double tend = hit.EndTime();
  int nelec = sc.NumberOfElectrons();
  for(int ie = 0; ie < nelec; ++ie) {
    const sim::Electrons* elec = sc.GetElectrons(ie);
    assert(elec != 0);
    double tdc = floor(elec->ArrivalT()/fSamplingRate) + fTriggerOffset;
    if(fMCHist) {
      geo::View_t view = hit.View();
      if(view == geo::kU)
	fHDTUE->Fill(tpeak - tdc);
      else if(view == geo::kV)
	fHDTVE->Fill(tpeak - tdc);
      else if(view == geo::kW)
	fHDTWE->Fill(tpeak - tdc);
      else
	throw cet::exception("SpacePointService") << "Bad view = " << view << "\n";
    }
    if(tdc >= tstart && tdc <= tend)
      electrons.push_back(elec);
  }
}

//----------------------------------------------------------------------
// Extract average position of vector of electrons.
// (Similar as BackTracker::HitToXYZ in MCCheater.)
TVector3 trkf::SpacePointService::ElectronsToXYZ(const std::vector<const sim::Electrons*>& electrons) const
{
  double w = 0.;
  double x = 0.;
  double y = 0.;
  double z = 0.;

  // Loop over electrons.

  for(std::vector<const sim::Electrons*>::const_iterator ie = electrons.begin();
      ie != electrons.end(); ++ie) {
    const sim::Electrons* elec = *ie;
    double weight = elec->NumElectrons();
    TVector3 pos = elec->XYZ();
    w += weight;
    x += weight * pos.x();
    y += weight * pos.y();
    z += weight * pos.z();
  }
  if(w != 0.) {
    x /= w;
    y /= w;
    z /= w;
  }
  TVector3 result(x, y, z);
  return result;
}

