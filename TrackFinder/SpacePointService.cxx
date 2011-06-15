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
#include <cmath>
#include <map>
#include <algorithm>
#include "cetlib/exception.h"
#include "TrackFinder/SpacePointService.h"
#include "Geometry/geo.h"
#include "Utilities/LArProperties.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "Simulation/SimChannel.h"
#include "Simulation/Electrons.h"
#include "Simulation/LArVoxelData.h"
#include "art/Framework/Core/Event.h"
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
  fDebug(0),
  fHist(false),
  fMCHist(false),
  fUseMC(false),
  fMaxDT(0.),
  fMaxS(0.),
  fMinViews(1000),
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
  for(int i=0; i<3; ++i) {
    fTimeOffset[i] = 0.;
    fEnable[i] = false;
    fWirePitch[i] = 0.;
    fWireOffset[i] = 0.;
    fTheta[i] = 0.;
    fSinTheta[i] = 0.;
    fCosTheta[i] = 0.;
    fSin[i] = 0.;
  }
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
  // Verify view numbering convention that views are numbered 1-3.

  assert(geo::kU >= 1 && geo::kU <= 3);
  assert(geo::kV >= 1 && geo::kV <= 3);
  assert(geo::kW >= 1 && geo::kW <= 3);
  assert(geo::kU != geo::kV);
  assert(geo::kU != geo::kW);
  assert(geo::kV != geo::kW);

  // Get configuration parameters.

  fDebug = pset.get<int>("Debug", 0);
  fHist = pset.get<bool>("Hist", false);
  fMCHist = pset.get<bool>("MCHist", false);
  fUseMC = pset.get<bool>("UseMC", false);
  fMClabel = pset.get<std::string>("MClabel", "daq");
  fMaxDT = pset.get<double>("MaxDT", 0.);
  fMaxS = pset.get<double>("MaxS", 0.);

  double timeOffsetU = pset.get<double>("TimeOffsetU", 0.);
  double timeOffsetV = pset.get<double>("TimeOffsetV", 0.);
  double timeOffsetW = pset.get<double>("TimeOffsetW", 0.);
  fTimeOffset[geo::kU-1] = timeOffsetU;
  fTimeOffset[geo::kV-1] = timeOffsetV;
  fTimeOffset[geo::kW-1] = timeOffsetW;

  fMinViews = pset.get<int>("MinViews", 1000);

  bool enableU = pset.get<bool>("EnableU", false);
  bool enableV = pset.get<bool>("EnableV", false);
  bool enableW = pset.get<bool>("EnableW", false);
  fEnable[geo::kU-1] = enableU;
  fEnable[geo::kV-1] = enableV;
  fEnable[geo::kW-1] = enableW;

  fGDMLPath = "";
  fROOTPath = "";

  // Report.

  std::cout << "\nSpacePointService configured with the following parameters:\n"
	    << "  Debug = " << fDebug << "\n"
	    << "  Hist = " << fHist << "\n"
	    << "  MC Hist = " << fMCHist << "\n"
	    << "  UseMC = " << fUseMC << "\n"
	    << "  fMClabel = " << fMClabel << "\n"
	    << "  MaxDT = " << fMaxDT << "\n"
	    << "  MaxS = " << fMaxS << "\n"
	    << "  TiemOffsetU = " << timeOffsetU << "\n"
	    << "  TiemOffsetV = " << timeOffsetV << "\n"
	    << "  TiemOffsetW = " << timeOffsetW << "\n" 
	    << "  MinViews = " << fMinViews << "\n"
	    << "  EnableU = " << enableU << "\n"
	    << "  EnableV = " << enableV << "\n"
	    << "  EnableW = " << enableW << "\n"
	    << std::endl;
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
  // First reset all geometry constants to null values.
  // We do this to make sure we don't accidentally inherit geomerty
  // constants from the previous geometry.

  for(int i=0; i<3; ++i) {
    fWirePitch[i] = 0.;
    fWireOffset[i] = 0.;
    fTheta[i] = 0.;
    fSinTheta[i] = 0.;
    fCosTheta[i] = 0.;
    fSin[i] = 0.;
  }

  // Calculate and print geometry information.

  std::cout << "SpacePointService updating geometry constants." << std::endl;

  int n = fGeom->Nplanes();
  for(int i=0; i<n; ++i) {
    const geo::PlaneGeo& pgeom = fGeom->Plane(i);
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

    int view = pgeom.View();
    assert(view >= 1 && view <= 3);
    fTheta[view-1] = theta;
    fSinTheta[view-1] = std::sin(theta);
    fCosTheta[view-1] = std::cos(theta);
    fWirePitch[view-1] = pitch;
    fWireOffset[view-1] = offset;

    std::cout << "Plane " << i << "\n"
	      << "  View " << view << "\n"
	      << "  SignalType " << pgeom.SignalType() << "\n"
	      << "  Orientation " << pgeom.Orientation() << "\n"
	      << "  Theta = " << theta << "\n"
	      << "  Wire pitch = " << pitch << "cm\n"
	      << "  Wire offset = " << offset << "\n";
  }

  fSin[0] = std::sin(fTheta[1] - fTheta[2]);
  fSin[1] = std::sin(fTheta[2] - fTheta[0]);
  fSin[2] = std::sin(fTheta[0] - fTheta[1]);

  std::cout << "  sin(V-W) = " << fSin[0] << "\n"
	    << "  sin(W-U) = " << fSin[1] << "\n"
	    << "  sin(U-V) = " << fSin[2] << "\n";

  // Update detector properties.

  art::ServiceHandle<util::DetectorProperties> detprop;
  fSamplingRate = detprop->SamplingRate();
  fTriggerOffset = detprop->TriggerOffset();
  std::cout << "\nDetector propertoes:\n"
	    << "  Sampling Rate = " << fSamplingRate << " ns/tick\n"
	    << "  Trigger offset = " << fTriggerOffset << " ticks\n" << std::endl;

  // Update LArProperties.

  art::ServiceHandle<util::LArProperties> larprop;
  fEfield = larprop->Efield();
  fTemperature = larprop->Temperature();
  fDriftVelocity = larprop->DriftVelocity(fEfield, fTemperature);
  fTimePitch = 0.001 * fDriftVelocity * fSamplingRate;
  std::cout << "\nLAr propertoes:\n"
	    << "  E field = " << fEfield << " kV/cm\n"
	    << "  Temperature = " << fTemperature << " K\n"
	    << "  Drift velocity = " << fDriftVelocity << " cm/us\n"
	    << "  Time pitch = " << fTimePitch << " cm/tick\n" << std::endl;

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
    fHDTUV = dir.make<TH1F>("DTUV", "U-V time difference", 100, -50, 50);
    fHDTVW = dir.make<TH1F>("DTVW", "V-W time difference", 100, -50, 50);
    fHDTWU = dir.make<TH1F>("DTWU", "W-U time difference", 100, -50, 50);
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

  if(result) {

    // First do pairwise tests.
    // Do double loop over hits.

    for(int ihit1 = 0; ihit1 < nhits-1; ++ihit1) {
      const recob::Hit& hit1 = *(hits[ihit1]);
      geo::View_t view1 = hit1.View();
      if(view1 < 1 || view1 > 3)
	throw cet::exception("SPTError") << "Bad view = " << view1 << "\n";
      double t1 = hit1.PeakTime() - fTimeOffset[view1-1];

      // If using mc information, make a collection of electrons for hit 1.

      const std::vector<const sim::Electrons*>& electrons1 = fElecMap[&hit1];

      // Loop over second hit.

      for(int ihit2 = ihit1+1; ihit2 < nhits; ++ihit2) {
	const recob::Hit& hit2 = *(hits[ihit2]);
	geo::View_t view2 = hit2.View();
	if(view2 < 1 || view2 > 3)
	  throw cet::exception("SPTError") << "Bad view = " << view2 << "\n";
	double t2 = hit2.PeakTime() - fTimeOffset[view2-1];
    
	// Test for different views and maximum time difference.

	result = view1 != view2 && std::abs(t1-t2) <= fMaxDT;

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

    // If there are exactly three hits, and they pass pairwise tests, check
    // for spatial compatibility.

    if(result && nhits == 3) {

      geo::View_t view[3];
      double time[3];
      unsigned int wire[3];
      unsigned int plane[3];
      unsigned int tpc;
      for(int i=0; i<3; ++i) {
	geo::View_t v = hits[i]->View();
	assert(v >= 1 && v <= 3);
	view[i] = v;
	time[i] = hits[i]->PeakTime() - fTimeOffset[v-1];
	unsigned short channel = hits[i]->Channel();
	fGeom->ChannelToWire(channel, tpc, plane[i], wire[i]);
      }

      // Get distance with offset correction.

      double dist[3] = {0., 0., 0.};
      for(int i = 0; i < 3; ++i) {
	geo::View_t v = view[i];
	assert(v >= 1 && v <= 3);
	dist[v-1] = fWirePitch[v-1] * wire[i] + fWireOffset[v-1];
      }

      // Do space cut.

      double S = fSin[0]*dist[0] + fSin[1]*dist[1] + fSin[2]*dist[2];
      result = std::abs(S) < fMaxS;

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

    geo::View_t v = hit.View();
    assert(v >= 1 && v <= 3);

    // Correct time for trigger offset and view-dependent time offsets.
    // Assume time error is proportional to (end time - start time).

    double t = hit.PeakTime() - fTriggerOffset - fTimeOffset[v-1];
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

      geo::View_t v = hit.View();
      assert(v >= 1 && v <= 3);

      // Calculate wire coordinate in this view.
    
      unsigned int plane;
      unsigned int wire;
      unsigned int tpc;
      unsigned short channel = hit.Channel();
      fGeom->ChannelToWire(channel, tpc, plane, wire);
      double u = wire * fWirePitch[v-1] + fWireOffset[v-1];

      // \todo probably need to account for different TPCs somehow here.

      // Summations

      double s = fSinTheta[v-1];
      double c = fCosTheta[v-1];
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

    // Debugging printout.

    if(fDebug >= 2) {
      std::cout << "\nmc   x=" << mcxyz[0] 
		<< ", y=" << mcxyz[1] 
		<< ", z=" << mcxyz[2] << std::endl;
      std::cout << "reco x=" << xyz[0]
		<< ", y=" << xyz[1] 
		<< ", z=" << xyz[2] << std::endl;
    }
  }
}

//----------------------------------------------------------------------
// Fill a vector of space points for all compatible combinations of hits
// from an input vector of hits.
//
void trkf::SpacePointService::makeSpacePoints(const art::Handle< std::vector<recob::Hit> >& hith,
					      std::vector<recob::SpacePoint>& spts,
					      const art::Event* pevt) const
{
  art::PtrVector<recob::Hit> hits;
  int nhits = hith->size();
  hits.reserve(nhits);

  for(int i = 0; i < nhits; ++i)
    hits.push_back(art::Ptr<recob::Hit>(hith, i));

  makeSpacePoints(hits, spts, pevt);
}

//----------------------------------------------------------------------
// Fill a vector of space points for all compatible combinations of hits
// from an input vector of hits.
//
// See remarks for previous method.
//
void trkf::SpacePointService::makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
					      std::vector<recob::SpacePoint>& spts,
					      const art::Event* pevt) const
{
  // First make result vector is empty.

  spts.erase(spts.begin(), spts.end());

  // Statistics.

  int n2 = 0;  // Number of two-hit space points.
  int n3 = 0;  // Number of three-hit space points.

  // If fUseMC pr fMCHist is true, extract SimChannel's from event.

  fSCHandle.clear();
  assert(!fSCHandle.isValid());
  if(fUseMC || fMCHist) {
    if(pevt == 0)
      throw cet::exception("SPTError") << "Event can not be null when requesting to use MC.\n";

    // Get SimChannels.

    pevt->getByLabel(fMClabel, fSCHandle);
    if(!fSCHandle.isValid())
      throw cet::exception("SPTError") << "Did not find any SimChannel's.\n";

    unsigned int nsc = fSCHandle->size();
    if(fDebug) {
      std::cout << "SpacePointService:\n"
		<< "  Found " << nsc << " SimChannel's." << std::endl;
    }
    if(nsc != fGeom->Nchannels()) {
      throw cet::exception("SPTError") 
	<< "Number of channels does not agree with geometry: " 
	<< fGeom->Nchannels() << "\n";
    }

    // Loop over SimChannels.  Verify that channels are sorted by channel
    // number.

    for(unsigned int isc = 0; isc < nsc; ++isc) {
      const sim::SimChannel& sc = (*fSCHandle)[isc];
      if(isc != sc.Channel())
	throw cet::exception("SPTError") << "MC channels not sorted.\n";
    }
  }

  // Sort hits into maps indexed by time for each view.
  // If using mc information, also generate maps of electrons and mc 
  // position indexed by hit.

  std::map<double, art::Ptr<recob::Hit> > hitmap[3];
  fElecMap.clear();
  fMCPosMap.clear();

  for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
      ihit != hits.end(); ++ihit) {
    const art::Ptr<recob::Hit>& phit = *ihit;
    geo::View_t view = phit->View();
    if(view < 1 || view > 3)
      throw cet::exception("SPTError") << "Bad view = " << view << "\n";
    if(fEnable[view-1]) {
      double t = phit->PeakTime() - fTimeOffset[view-1];
      hitmap[view-1][t] = phit;
      const recob::Hit& hit = *phit;

      // Get Electrons.

      if(fUseMC || fMCHist) {
	std::vector<const sim::Electrons*>& electrons = fElecMap[&hit];   // Empty vector.
	HitToElectrons(hit, electrons);
	fMCPosMap[&hit] = ElectronsToXYZ(electrons);
	std::sort(electrons.begin(), electrons.end(), ElectronsLess());
	std::vector<const sim::Electrons*>::iterator it = 
	  std::unique(electrons.begin(), electrons.end(), ElectronsEqual());
	electrons.resize(it - electrons.begin());
      }
    }
  }

  if(fDebug) {
    std::cout << "\nSpacePointService:\n"
	      << "  Total hits = " << hits.size() << "\n"
	      << "  U hits = " << hitmap[0].size() << "\n"
	      << "  V hits = " << hitmap[1].size() << "\n"
	      << "  W hits = " << hitmap[2].size() << "\n"
	      << std::endl;
  }

  // Sort maps in increasing order of number of hits.
  // This is so that we can do the outer loops over hits 
  // over the views with fewer hits.

  int index[3] = {0, 1, 2};
  for(int i=0; i<2; ++i) {
    for(int j=i+1; j<3; ++j) {
      if(hitmap[index[i]].size() > hitmap[index[j]].size()) {
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
      int index1 = index[i];
      for(int j=i+1; j<3; ++j) {
	int index2 = index[j];

	assert(hitmap[index1].size() <= hitmap[index2].size());

	// Loop over pairs of hits.

	art::PtrVector<recob::Hit> hitvec;
	hitvec.reserve(2);

	for(std::map<double, art::Ptr<recob::Hit> >::const_iterator ihit1 = hitmap[index1].begin();
	    ihit1 != hitmap[index1].end(); ++ihit1) {

	  const art::Ptr<recob::Hit>& phit1 = ihit1->second;
	  assert(phit1->View() == index1+1);

	  double t1 = phit1->PeakTime() - fTimeOffset[index1];
	  double t2min = t1 - fMaxDT;
	  double t2max = t1 + fMaxDT;

	  for(std::map<double, art::Ptr<recob::Hit> >::const_iterator 
		ihit2 = hitmap[index2].lower_bound(t2min);
	      ihit2 != hitmap[index2].upper_bound(t2max); ++ihit2) {

	    const art::Ptr<recob::Hit>& phit2 = ihit2->second;
	    assert(phit2->View() == index2+1);

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

    for(std::map<double, art::Ptr<recob::Hit> >::const_iterator ihit1 = hitmap[index[0]].begin();
	ihit1 != hitmap[index[0]].end(); ++ihit1) {

      const art::Ptr<recob::Hit>& phit1 = ihit1->second;
      assert(phit1->View() == index[0]+1);

      double t1 = phit1->PeakTime() - fTimeOffset[index[0]];
      double t2min = t1 - fMaxDT;
      double t2max = t1 + fMaxDT;

      for(std::map<double, art::Ptr<recob::Hit> >::const_iterator 
	    ihit2 = hitmap[index[1]].lower_bound(t2min);
	  ihit2 != hitmap[index[1]].upper_bound(t2max); ++ihit2) {

	const art::Ptr<recob::Hit>& phit2 = ihit2->second;
	assert(phit2->View() == index[1]+1);

	// Test first two hits for compatibility before looping 
	// over third hit.

	hitvec.clear();
	hitvec.push_back(phit1);
	hitvec.push_back(phit2);
	bool ok = compatible(hitvec);
	if(ok) {

	  double t2 = phit2->PeakTime() - fTimeOffset[index[1]];
	  double t3min = std::max(t1, t2) - fMaxDT;
	  double t3max = std::min(t1, t2) + fMaxDT;

	  for(std::map<double, art::Ptr<recob::Hit> >::const_iterator 
		ihit3 = hitmap[index[2]].lower_bound(t3min);
	      ihit3 != hitmap[index[2]].upper_bound(t3max); ++ihit3) {

	    const art::Ptr<recob::Hit>& phit3 = ihit3->second;
	    assert(phit3->View() == index[2]+1);

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

  if(fDebug) {
    std::cout << "SpacePointService:\n"
	      << "  2-hit space points = " << n2 << "\n"
	      << "  3-hit space points = " << n3 << "\n" 
	      << std::endl;
  }
}

//----------------------------------------------------------------------
// Extract Electrons associated with Hit.
// (Similar as BackTracker::HitToXYZ in MCCheater.)
void trkf::SpacePointService::HitToElectrons(const recob::Hit& hit,
					     std::vector<const sim::Electrons*>& electrons) const
{
  // Find SimChannel associated with Hit.

  unsigned short channel = hit.Channel();
  const sim::SimChannel& sc = (*fSCHandle)[channel];
  assert(sc.Channel() == channel);

  // Find electrons in time with Hit.

  double tstart = hit.StartTime();
  double tend = hit.EndTime();
  int nelec = sc.NumberOfElectrons();
  for(int ie = 0; ie < nelec; ++ie) {
    const sim::Electrons* elec = sc.GetElectrons(ie);
    assert(elec != 0);
    double tdc = floor(elec->ArrivalT()/fSamplingRate) + fTriggerOffset;
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

