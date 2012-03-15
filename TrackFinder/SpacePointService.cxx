///////////////////////////////////////////////////////////////////////
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
#include "art/Framework/Principal/View.h"
#include "Utilities/DetectorProperties.h"
#include "MCCheater/BackTracker.h"
#include "TH1F.h"

namespace  {

  // Function classes for sorting sim::IDEs according to track id.

  //class IDELess {
  //public:
  //  bool operator()(const sim::IDE& p1, const sim::IDE& p2) {
  //    bool result = p1.trackID < p2.trackID;
  //    return result;
  //  }
  //};

  //class IDEEqual {
  //public:
  //  bool operator()(const sim::IDE& p1, const sim::IDE& p2) {
  //    bool result = p1.trackID == p2.trackID;
  //    return result;
  //  }
  //};

}

//----------------------------------------------------------------------
// Constructor.
//
namespace  trkf{
  
  SpacePointService::SpacePointService(const fhicl::ParameterSet& pset,
				       art::ActivityRegistry& reg) :
    fMaxDT(0.),
    fMaxS(0.),
    fTimeOffsetU(0.),
    fTimeOffsetV(0.),
    fTimeOffsetW(0.),
    fMinViews(1000),
    fEnableU(false),
    fEnableV(false),
    fEnableW(false),
    fFilter(false),
    fMerge(false)
  {
    reconfigure(pset);
  }
  
  //----------------------------------------------------------------------
  // Destructor.
  //
  SpacePointService::~SpacePointService()
  {
  }

  //----------------------------------------------------------------------
  // Update configuration parameters.
  //
  void SpacePointService::reconfigure(const fhicl::ParameterSet& pset)
  {
    // Get configuration parameters.

    fMaxDT = pset.get<double>("MaxDT", 0.);
    fMaxS = pset.get<double>("MaxS", 0.);

    fTimeOffsetU = pset.get<double>("TimeOffsetU", 0.);
    fTimeOffsetV = pset.get<double>("TimeOffsetV", 0.);
    fTimeOffsetW = pset.get<double>("TimeOffsetW", 0.);

    fMinViews = pset.get<int>("MinViews", 1000);

    fEnableU = pset.get<bool>("EnableU", false);
    fEnableV = pset.get<bool>("EnableV", false);
    fEnableW = pset.get<bool>("EnableW", false);
    fFilter = pset.get<bool>("Filter", false);
    fMerge = pset.get<bool>("Merge", false);

    // Only allow one of fFilter and fMerge to be true.

    if(fFilter && fMerge)
      throw cet::exception("SpacePointService") << "Filter and Merge flags are both true.\n";

    // Report.

    mf::LogInfo("SpacePointService") 
      << "SpacePointService configured with the following parameters:\n"
      << "  MaxDT = " << fMaxDT << "\n"
      << "  MaxS = " << fMaxS << "\n"
      << "  TimeOffsetU = " << fTimeOffsetU << "\n"
      << "  TimeOffsetV = " << fTimeOffsetV << "\n"
      << "  TimeOffsetW = " << fTimeOffsetW << "\n" 
      << "  MinViews = " << fMinViews << "\n"
      << "  EnableU = " << fEnableU << "\n"
      << "  EnableV = " << fEnableV << "\n"
      << "  EnableW = " << fEnableW << "\n"
      << "  Filter = " << fFilter << "\n"
      << "  Merge = " << fMerge;
  }

  //----------------------------------------------------------------------
  // Print geometry and properties constants.
  //
  void SpacePointService::update() const
  {
    // Generate info report on first call only.

    static bool first = true;
    bool report = first;
    first = false;

    // Get services.

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detprop;
    art::ServiceHandle<util::LArProperties> larprop;

    // Calculate and print geometry information.

    mf::LogInfo log("SpacePointService");
    if(report)
      log << "Updating geometry constants.\n";

    // Update detector properties.

    double samplingRate = detprop->SamplingRate();
    double triggerOffset = detprop->TriggerOffset();
    if(report) {
      log << "\nDetector properties:\n"
	  << "  Sampling Rate = " << samplingRate << " ns/tick\n"
	  << "  Trigger offset = " << triggerOffset << " ticks\n";
    }
  
    // Update LArProperties.
  
    double efield = larprop->Efield();
    double temperature = larprop->Temperature();
    double driftVelocity = larprop->DriftVelocity(efield, temperature);
    double timePitch = 0.001 * driftVelocity * samplingRate;
    if(report) {
      log << "\nLAr propertoes:\n"
	  << "  E field = " << efield << " kV/cm\n"
	  << "  Temperature = " << temperature << " K\n"
	  << "  Drift velocity = " << driftVelocity << " cm/us\n"
	  << "  Time pitch = " << timePitch << " cm/tick\n";
    }
  
    // Get time offsets.
  
    std::vector<std::vector<std::vector<double> > > timeOffset;
    fillTimeOffset(timeOffset);
  
    for(unsigned int cstat = 0; cstat < geom->Ncryostats(); ++cstat){
    
      // Loop over TPCs.
    
      int ntpc = geom->Cryostat(cstat).NTPC();
    
      for(int tpc = 0; tpc < ntpc; ++tpc) {
	const geo::TPCGeo& tpcgeom = geom->Cryostat(cstat).TPC(tpc);
      
	// Loop over planes.
      
	int nplane = tpcgeom.Nplanes();
      
	for(int plane = 0; plane < nplane; ++plane) {
	  const geo::PlaneGeo& pgeom = tpcgeom.Plane(plane);
	
	  // Fill view-dependent quantities.
	
	  geo::View_t view = pgeom.View();
	  std::string viewname = "?";
	  if(view == geo::kU) {
	    viewname = "U";
	  }
	  else if(view == geo::kV) {
	    viewname = "V";
	  }
	  else if(view == geo::kW) {
	    viewname = "W";
	  }
	  else
	    throw cet::exception("SpacePointService") << "Bad view = " 
						      << view << "\n";
	
	  std::string sigtypename = "?";
	  geo::SigType_t sigtype = pgeom.SignalType();
	  if(sigtype == geo::kInduction)
	    sigtypename = "Induction";
	  else if(sigtype == geo::kCollection)
	    sigtypename = "Collection";
	  else
	    throw cet::exception("SpacePointService") << "Bad signal type = " 
						      << sigtype << "\n";
	
	  std::string orientname = "?";
	  geo::Orient_t orient = pgeom.Orientation();
	  if(orient == geo::kVertical)
	    orientname = "Vertical";
	  else if(orient == geo::kHorizontal)
	    orientname = "Horizontal";
	  else
	    throw cet::exception("SpacePointService") << "Bad orientation = " 
						      << orient << "\n";
	
	  if(report) {
	    const double* xyz = tpcgeom.PlaneLocation(plane);
	    log << "\nCryostat, TPC, Plane: " << cstat << "," << tpc << ", " << plane << "\n"
		<< "  View: " << viewname << "\n"
		<< "  SignalType: " << sigtypename << "\n"
		<< "  Orientation: " << orientname << "\n"
		<< "  Plane location: " << xyz[0] << "\n"
		<< "  Plane pitch: " << tpcgeom.Plane0Pitch(plane) << "\n"
		<< "  Wire angle: " << tpcgeom.Plane(plane).Wire(0).ThetaZ() << "\n"
		<< "  Wire pitch: " << tpcgeom.WirePitch() << "\n"
		<< "  Time offset: " << timeOffset[cstat][tpc][plane] << "\n";
	  }
	
	  if(orient != geo::kVertical)
	    throw cet::exception("SpacePointService") 
	      << "Horizontal wire geometry not implemented.\n";
	}// end loop over planes
      }// end loop over tpcs
    }// end loop over cryostats
  }
 
  //----------------------------------------------------------------------
  // Calculate time offsets.
  // Results stored in nested vector indexed by [cryostat][tpc][plane]
  void SpacePointService::fillTimeOffset(std::vector< std::vector<std::vector<double> > >& timeOffset) const
  {
    // Get services.

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detprop;
    art::ServiceHandle<util::LArProperties> larprop;

    // Clear result.

    timeOffset.clear();

    // Get properties needed to calculate time offsets.

    double samplingRate = detprop->SamplingRate();
    double triggerOffset = detprop->TriggerOffset();
    double efield = larprop->Efield();
    double temperature = larprop->Temperature();
    double driftVelocity = larprop->DriftVelocity(efield, temperature);
    double timePitch = 0.001 * driftVelocity * samplingRate;

    // Loop over TPCs.
    timeOffset.resize(geom->Ncryostats());

    for(size_t cstat = 0; cstat < geom->Ncryostats(); ++cstat){
      timeOffset[cstat].resize(geom->Cryostat(cstat).NTPC());
    
      for(size_t tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc) {
	const geo::TPCGeo& tpcgeom = geom->Cryostat(cstat).TPC(tpc);

	// Loop over planes.
      
	int nplane = tpcgeom.Nplanes();
	timeOffset[cstat][tpc].resize(nplane, 0.);

	for(int plane = 0; plane < nplane; ++plane) {
	  const geo::PlaneGeo& pgeom = tpcgeom.Plane(plane);

	  // Calculate geometric time offset.
	
	  const double* xyz = tpcgeom.PlaneLocation(0);
	  timeOffset[cstat][tpc][plane] =
	    (-xyz[0] + tpcgeom.Plane0Pitch(plane)) / timePitch + triggerOffset;
	
	  // Add view-dependent time offset.
	
	  geo::View_t view = pgeom.View();
	  if(view == geo::kU)
	    timeOffset[cstat][tpc][plane] += fTimeOffsetU;
	  else if(view == geo::kV)
	    timeOffset[cstat][tpc][plane] += fTimeOffsetV;
	  else if(view == geo::kW)
	    timeOffset[cstat][tpc][plane] += fTimeOffsetW;
	  else
	    throw cet::exception("SpacePointService") << "Bad view = " 
						      << view << "\n";
	}
      }
    }// end loop over cryostats

    return;
  }



  //----------------------------------------------------------------------
  // Get corrected time for the specified hit.
  double SpacePointService::correctedTime(const recob::Hit& hit,
					  const std::vector< std::vector<std::vector<double> > >& timeOffset) const
  {
    // Get services.

    art::ServiceHandle<geo::Geometry> geom;

    // Get tpc, plane.

    unsigned short channel = hit.Channel();
    unsigned int tpc, plane, wire, cstat;
    geom->ChannelToWire(channel, cstat, tpc, plane, wire);

    // Correct time for trigger offset and plane-dependent time offsets.

    double t = hit.PeakTime() - timeOffset[cstat][tpc][plane];

    return t;
  }

  //----------------------------------------------------------------------
  // Spatial separation of hits (zero if two or fewer).
  double SpacePointService::separation(const art::PtrVector<recob::Hit>& hits) const
  {
    // Get geometry service.

    art::ServiceHandle<geo::Geometry> geom;

    // Trivial case - fewer than three hits.

    if(hits.size() < 3)
      return 0.;

    // Error case - more than three hits.

    if(hits.size() > 3) {
      mf::LogError("SpacePointService") << "Method separation called with more than three htis.";
      return 0.;
    }

    // Got exactly three hits.

    assert(hits.size() == 3);

    // Calculate angles and distance of each hit from origin.

    double dist[3] = {0., 0., 0.};
    double sinth[3] = {0., 0., 0.};
    double costh[3] = {0., 0., 0.};
    unsigned int cstats[3];
    unsigned int tpcs[3];
    unsigned int planes[3];

    for(int i=0; i<3; ++i) {

      // Get tpc, plane, wire.

      const recob::Hit& hit = *(hits[i]);
      unsigned short channel = hit.Channel();
      unsigned int tpc, plane, wire, cstat;
      const geo::WireGeo& wgeom = geom->ChannelToWire(channel, cstat, tpc, plane, wire);
      tpcs[i] = tpc;
      planes[i] = plane;

      // Check tpc and plane errors.

      for(int j=0; j<i; ++j) {

	if(cstats[j] != cstat) {
	  mf::LogError("SpacePointService") << "Method separation called with hits from multiple cryostats..";
	  return 0.;
	}

	if(tpcs[j] != tpc) {
	  mf::LogError("SpacePointService") << "Method separation called with hits from multiple tpcs..";
	  return 0.;
	}

	if(planes[j] == plane) {
	  mf::LogError("SpacePointService") << "Method separation called with hits from the same plane..";
	  return 0.;
	}
      }

      // Get angles and distance of wire.

      double hl = wgeom.HalfL();
      double xyz[3];
      double xyz1[3];
      wgeom.GetCenter(xyz);
      wgeom.GetCenter(xyz1, hl);
      double s = (xyz1[1] - xyz[1]) / hl;
      double c = (xyz1[2] - xyz[2]) / hl;
      sinth[plane] = s;
      costh[plane] = c;
      dist[plane] = xyz[2] * s - xyz[1] * c;
    }

    double S = ((sinth[1] * costh[2] - costh[1] * sinth[2]) * dist[0] 
		+(sinth[2] * costh[0] - costh[2] * sinth[0]) * dist[1] 
		+(sinth[0] * costh[1] - costh[0] * sinth[1]) * dist[2]);
    return S;
  }

  //----------------------------------------------------------------------
  // Check hits for compatibility.
  // Check hits pairwise for different views and maximum time difference.
  // Check three hits for spatial compatibility.
  bool SpacePointService::compatible(const art::PtrVector<recob::Hit>& hits,
				     const std::vector<std::vector<std::vector<double> > > & timeOffset,
				     bool useMC,
				     double maxDT, double maxS) const
  {
    // Get geometry service.

    art::ServiceHandle<geo::Geometry> geom;

    // Get cuts.

    if(maxDT == 0.)
      maxDT = fMaxDT;
    if(maxS == 0.)
      maxS = fMaxS;  

    int nhits = hits.size();

    // Fewer than two or more than three hits can never be compatible.

    bool result = nhits >= 2 && nhits <= 3;
    bool mc_ok = true;
    unsigned int tpc = 0;
    unsigned int cstat = 0;

    if(result) {

      // First do pairwise tests.
      // Do double loop over hits.

      for(int ihit1 = 0; result && ihit1 < nhits-1; ++ihit1) {
	const recob::Hit& hit1 = *(hits[ihit1]);
	unsigned short channel1 = hit1.Channel();
	unsigned int tpc1, plane1, wire1, cstat1;
	geom->ChannelToWire(channel1, cstat1, tpc1, plane1, wire1);
	geo::View_t view1 = hit1.View();
	double t1 = hit1.PeakTime() - timeOffset[cstat1][tpc1][plane1];

	// If using mc information, get a collection of track ids for hit 1.
	// If not using mc information, this section of code will trigger the 
	// insertion of a single invalid HitMCInfo object into fHitMCMap.

	const HitMCInfo& mcinfo1 = fHitMCMap[(useMC ? &hit1 : 0)];
	const std::vector<int>& tid1 = mcinfo1.trackIDs;
	bool only_neg1 = tid1.size() > 0 && tid1.back() < 0;

	// Loop over second hit.

	for(int ihit2 = ihit1+1; result && ihit2 < nhits; ++ihit2) {
	  const recob::Hit& hit2 = *(hits[ihit2]);
	  unsigned short channel2 = hit2.Channel();
	  unsigned int tpc2, plane2, wire2, cstat2;
	  geom->ChannelToWire(channel2, cstat2, tpc2, plane2, wire2);
	  geo::View_t view2 = hit2.View();

	  // Test for same tpc and different views.

	  result = result && tpc1 == tpc2 && view1 != view2 && cstat1 == cstat2;
	  if(result) {

	    // Remember which tpc and cryostat we are in.

	    tpc = tpc1;
	    cstat = cstat1;

	    double t2 = hit2.PeakTime() - timeOffset[cstat2][tpc2][plane2];
    
	    // Test maximum time difference.

	    result = result && std::abs(t1-t2) <= maxDT;

	    // Test mc truth.

	    if(result && useMC) {

	      // Test whether hits have a common parent track id.

	      const HitMCInfo& mcinfo2 = fHitMCMap[&hit2];
	      std::vector<int> tid2 = mcinfo2.trackIDs;
	      bool only_neg2 = tid2.size() > 0 && tid2.back() < 0;
	      std::vector<int>::iterator it =
		std::set_intersection(tid1.begin(), tid1.end(),
				      tid2.begin(), tid2.end(),
				      tid2.begin());
	      tid2.resize(it - tid2.begin());

	      // Hits are compatible if they have parents in common.
	      // If the only parent id in common is negative (-999),
	      // then hits are compatible only if both hits have only
	      // negative parent tracks.

	      bool only_neg3 = tid2.size() > 0 && tid2.back() < 0;
	      mc_ok = tid2.size() > 0 && 
		(!only_neg3 || (only_neg1 && only_neg2));
	      result = result && mc_ok;

	      // If we are still OK, check that either hit is
	      // the nearest neighbor of the other.

	      if(result) {
		result = mcinfo1.pchit[plane2] == &hit2 || 
		  mcinfo2.pchit[plane1] == &hit1;
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
	double sinth[3] = {0., 0., 0.};
	double costh[3] = {0., 0., 0.};

	for(int i=0; i<3; ++i) {

	  // Get tpc, plane, wire.

	  const recob::Hit& hit = *(hits[i]);
	  unsigned short channel = hit.Channel();
	  unsigned int tpc0, plane, wire, cstat0;
	  const geo::WireGeo& wgeom = geom->ChannelToWire(channel, cstat0, tpc0, plane, wire);
	  assert(tpc0 == tpc && cstat0 == cstat);

	  // Get angles and distance of wire.

	  double hl = wgeom.HalfL();
	  double xyz[3];
	  double xyz1[3];
	  wgeom.GetCenter(xyz);
	  wgeom.GetCenter(xyz1, hl);
	  double s  = (xyz1[1] - xyz[1]) / hl;
	  double c = (xyz1[2] - xyz[2]) / hl;
	  sinth[plane] = s;
	  costh[plane] = c;
	  dist[plane] = xyz[2] * s - xyz[1] * c;
	}

	// Do space cut.

	double S = ((sinth[1] * costh[2] - costh[1] * sinth[2]) * dist[0] 
		    +(sinth[2] * costh[0] - costh[2] * sinth[0]) * dist[1] 
		    +(sinth[0] * costh[1] - costh[0] * sinth[1]) * dist[2]);

	result = result && std::abs(S) < maxS;
      }
    }

    // Done.
    
    return result;
  }

  //----------------------------------------------------------------------
  // Fill one space point using a colleciton of hits.
  // Assume points have already been tested for compatibility.
  //
  void SpacePointService::fillSpacePoint(const art::PtrVector<recob::Hit>& hits,
					 const std::vector< std::vector<std::vector<double> > >& timeOffset,
					 recob::SpacePoint& spt) const
  {
    // Get services.

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detprop;
    art::ServiceHandle<util::LArProperties> larprop;

    // Calculate time pitch.

    double efield = larprop->Efield();
    double temperature = larprop->Temperature();
    double driftVelocity = larprop->DriftVelocity(efield, temperature); // cm / us
    double samplingRate = detprop->SamplingRate();                      // ns
    double timePitch = 0.001 * driftVelocity * samplingRate;            // cm / tick

    // Store hits in SpacePoint.

    spt = recob::SpacePoint(recob::SpacePoint(hits));
    int nhits = hits.size();

    // Calculate position and error matrix.

    double xyz[3] = {0., 0., 0.};
    double errxyz[6] = {0.,
			0., 0.,
			0., 0., 0.};

    // Calculate x using drift times.
    // Loop over all hits and calculate the weighted average drift time.
    // Also calculate time variance and chisquare.

    double sumt2w = 0.;
    double sumtw = 0.;
    double sumw = 0.;

    for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
	ihit != hits.end(); ++ihit) {

      const recob::Hit& hit = **ihit;
      unsigned short channel = hit.Channel();
      unsigned int tpc, plane, wire, cstat;
      geom->ChannelToWire(channel, cstat, tpc, plane, wire);

      // Correct time for trigger offset and view-dependent time offsets.

      double t0 = timeOffset[cstat][tpc][plane];
      double t = hit.PeakTime() - t0;
      double et = hit.SigmaPeakTime();
      double w = 1./(et*et);

      sumt2w += w*t*t;
      sumtw += w*t;
      sumw += w;
    }

    double drift_time = 0.;
    double var_time = 0.;
    double chisq = 0.;
    if(sumw != 0.) {
      drift_time = sumtw / sumw;
      var_time = sumt2w / sumw - drift_time * drift_time;
      if(var_time < 0.)
	var_time = 0.;
      chisq = sumt2w - sumtw * drift_time;
      if(chisq < 0.)
	chisq = 0.;
    }
    xyz[0] = drift_time * timePitch;
    errxyz[0] = var_time * timePitch * timePitch;

    // Calculate y, z using wires (need at least two hits).

    if(nhits >= 2) {

      // Calculate y and z by chisquare minimization of wire coordinates.

      double sw = 0.;    // sum w_i
      double sus = 0.;   // sum w_i u_i sin_th_i
      double suc = 0.;   // sum w_i u_i cos_th_i
      double sc2 = 0.;   // sum w_i cos2_th_i
      double ss2 = 0.;   // sum w_i sin2_th_i
      double ssc = 0.;   // sum w_i sin_th_i cos_th_i

      // Loop over points.

      for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
	  ihit != hits.end(); ++ihit) {

	const recob::Hit& hit = **ihit;
	unsigned short channel = hit.Channel();
	unsigned int tpc, plane, wire, cstat;
	const geo::WireGeo& wgeom = geom->ChannelToWire(channel, cstat, tpc, plane, wire);

	// Calculate angle and wire coordinate in this view.
    
	double hl = wgeom.HalfL();
	double xyz[3];
	double xyz1[3];
	wgeom.GetCenter(xyz);
	wgeom.GetCenter(xyz1, hl);
	double s  = (xyz1[1] - xyz[1]) / hl;
	double c = (xyz1[2] - xyz[2]) / hl;
	double u = xyz[2] * s - xyz[1] * c;
	double eu = geom->WirePitch(0, 1, plane, tpc) / std::sqrt(12.);
	double w = 1. / (eu * eu);

	// Summations

	sw += w;
	sus += w*u*s;
	suc += w*u*c;
	sc2 += w*c*c;
	ss2 += w*s*s;
	ssc += w*s*c;
      }

      // Calculate y,z

      double denom = sc2 * ss2 - ssc * ssc;
      if(denom != 0.) {
	xyz[1] = (-suc * ss2 + sus * ssc) / denom;
	xyz[2] = (sus * sc2 - suc * ssc) / denom;
	errxyz[2] = ss2 / denom;
	errxyz[4] = ssc / denom;
	errxyz[5] = sc2 / denom;
      }

      // Set coordintates, error matrix, and chisquare in space point.

      spt.SetXYZ(xyz);
      spt.SetErrXYZ(errxyz);
      spt.SetChisq(chisq);
    }
    return;
  }

  //----------------------------------------------------------------------
  // Fill one space point using a colleciton of hits.
  // Assume points have already been tested for compatibility.
  // This version assumes there can be multiple hits per view,
  // and gives unequal weight to different hits.
  //
  void SpacePointService::
  fillComplexSpacePoint(const art::PtrVector<recob::Hit>& hits,
			const std::vector<std::vector<std::vector<double> > >& timeOffset,
			recob::SpacePoint& spt) const
  {
    // Get services.

    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::DetectorProperties> detprop;
    art::ServiceHandle<util::LArProperties> larprop;

    // Calculate time pitch.

    double efield = larprop->Efield();
    double temperature = larprop->Temperature();
    double driftVelocity = larprop->DriftVelocity(efield, temperature); // cm / us
    double samplingRate = detprop->SamplingRate();                      // ns
    double timePitch = 0.001 * driftVelocity * samplingRate;            // cm / tick

    // Figure out which tpc we are in.

    unsigned int tpc0 = 0;
    unsigned int cstat0 = 0;
    int nhits = hits.size();
    if(nhits > 0) {
      unsigned short channel = hits.front()->Channel();
      unsigned int plane, wire;
      geom->ChannelToWire(channel, cstat0, tpc0, plane, wire);
    }

    // Store hits in SpacePoint.

    spt = recob::SpacePoint(recob::SpacePoint(hits));

    // Do a preliminary scan of hits.
    // Determine weight given to hits in each view.

    unsigned int nplanes = geom->Cryostat(cstat0).TPC(tpc0).Nplanes();
    std::vector<int> numhits(nplanes, 0);
    std::vector<double> weight(nplanes, 0.);

    for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
	ihit != hits.end(); ++ihit) {

      const recob::Hit& hit = **ihit;
      unsigned short channel = hit.Channel();
      unsigned int tpc, plane, wire, cstat;
      geom->ChannelToWire(channel, cstat, tpc, plane, wire);
      assert(cstat == cstat0);
      assert(tpc == tpc0);
      assert(plane < nplanes);
      ++numhits[plane];
    }

    for(unsigned int plane = 0; plane < nplanes; ++plane) {
      double np = numhits[plane];
      if(np > 0.)
	weight[plane] = 1. / (np*np*np);
    }

    // Calculate position and error matrix.

    double xyz[3] = {0., 0., 0.};
    double errxyz[6] = {0.,
			0., 0.,
			0., 0., 0.};

    // Calculate x using drift times.
    // Loop over all hits and calculate the weighted average drift time.
    // Also calculate time variance and chisquare.

    double sumt2w = 0.;
    double sumtw = 0.;
    double sumw = 0.;

    for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
	ihit != hits.end(); ++ihit) {

      const recob::Hit& hit = **ihit;
      unsigned short channel = hit.Channel();
      unsigned int tpc, plane, wire, cstat;
      geom->ChannelToWire(channel, cstat, tpc, plane, wire);

      // Correct time for trigger offset and view-dependent time offsets.

      double t0 = timeOffset[cstat][tpc][plane];
      double t = hit.PeakTime() - t0;
      double et = hit.SigmaPeakTime();
      double w = weight[plane]/(et*et);

      sumt2w += w*t*t;
      sumtw += w*t;
      sumw += w;
    }

    double drift_time = 0.;
    double var_time = 0.;
    double chisq = 0.;
    if(sumw != 0.) {
      drift_time = sumtw / sumw;
      var_time = sumt2w / sumw - drift_time * drift_time;
      if(var_time < 0.)
	var_time = 0.;
      chisq = sumt2w - sumtw * drift_time;
      if(chisq < 0.)
	chisq = 0.;
    }
    xyz[0] = drift_time * timePitch;
    errxyz[0] = var_time * timePitch * timePitch;

    // Calculate y, z using wires (need at least two hits).

    if(nhits >= 2) {

      // Calculate y and z by chisquare minimization of wire coordinates.

      double sw = 0.;    // sum w_i
      double sus = 0.;   // sum w_i u_i sin_th_i
      double suc = 0.;   // sum w_i u_i cos_th_i
      double sc2 = 0.;   // sum w_i cos2_th_i
      double ss2 = 0.;   // sum w_i sin2_th_i
      double ssc = 0.;   // sum w_i sin_th_i cos_th_i

      // Loop over points.

      for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
	  ihit != hits.end(); ++ihit) {

	const recob::Hit& hit = **ihit;
	unsigned short channel = hit.Channel();
	unsigned int tpc, plane, wire, cstat;
	const geo::WireGeo& wgeom = geom->ChannelToWire(channel, cstat, tpc, plane, wire);

	// Calculate angle and wire coordinate in this view.
    
	double hl = wgeom.HalfL();
	double xyz[3];
	double xyz1[3];
	wgeom.GetCenter(xyz);
	wgeom.GetCenter(xyz1, hl);
	double s  = (xyz1[1] - xyz[1]) / hl;
	double c = (xyz1[2] - xyz[2]) / hl;
	double u = xyz[2] * s - xyz[1] * c;
	double eu = geom->WirePitch(0, 1, plane, tpc) / std::sqrt(12.);
	double w = weight[plane] / (eu * eu);

	// Summations

	sw += w;
	sus += w*u*s;
	suc += w*u*c;
	sc2 += w*c*c;
	ss2 += w*s*s;
	ssc += w*s*c;
      }

      // Calculate y,z

      double denom = sc2 * ss2 - ssc * ssc;
      if(denom != 0.) {
	xyz[1] = (-suc * ss2 + sus * ssc) / denom;
	xyz[2] = (sus * sc2 - suc * ssc) / denom;
	errxyz[2] = ss2 / denom;
	errxyz[4] = ssc / denom;
	errxyz[5] = sc2 / denom;
      }

      // Set coordintates, error matrix, and chisquare in space point.

      spt.SetXYZ(xyz);
      spt.SetErrXYZ(errxyz);
      spt.SetChisq(chisq);
    }
    return;
  }

  //----------------------------------------------------------------------
  // Fill a vector of space points for all compatible combinations of hits
  // from an input vector of hits (non-config-overriding, non-mc-truth version).
  //
  void SpacePointService::makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
					  std::vector<recob::SpacePoint>& spts) const
  {
    std::vector<const sim::SimChannel*> empty;
    makeSpacePoints(hits, spts, empty, false, fFilter, fMerge, fMaxDT, fMaxS);
  }

  //----------------------------------------------------------------------
  // Fill a vector of space points for all compatible combinations of hits
  // from an input vector of hits (config-overriding, non-mc-truth version).
  //
  void SpacePointService::makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
					  std::vector<recob::SpacePoint>& spts,
					  bool filter, bool merge,
					  double maxDT, double maxS) const
  {
    std::vector<const sim::SimChannel*> empty;
    makeSpacePoints(hits, spts, empty, false, filter, merge, maxDT, maxS);
  }

  //----------------------------------------------------------------------
  // Fill a vector of space points for all compatible combinations of hits
  // from an input vector of hits (non-config-overriding, mc-truth version).
  //
  void SpacePointService::makeMCTruthSpacePoints(const art::PtrVector<recob::Hit>& hits,
						 std::vector<recob::SpacePoint>& spts,
						 const std::vector<const sim::SimChannel*>& simchans) const
  {
    makeSpacePoints(hits, spts, simchans, true, fFilter, fMerge, fMaxDT, fMaxS);
  }

  //----------------------------------------------------------------------
  // Fill a vector of space points for all compatible combinations of hits
  // from an input vector of hits (config-overriding, mc-truth version).
  //
  void SpacePointService::makeMCTruthSpacePoints(const art::PtrVector<recob::Hit>& hits,
						 std::vector<recob::SpacePoint>& spts,
						 const std::vector<const sim::SimChannel*>& simchans,
						 bool filter, bool merge,
						 double maxDT, double maxS) const
  {
    makeSpacePoints(hits, spts, simchans, true, filter, merge, maxDT, maxS);
  }

  //----------------------------------------------------------------------
  // Fill a vector of space points for all compatible combinations of hits
  // from an input vector of hits (general version).
  //
  void SpacePointService::makeSpacePoints(const art::PtrVector<recob::Hit>& hits,
					  std::vector<recob::SpacePoint>& spts,
					  const std::vector<const sim::SimChannel*>& simchans,
					  bool useMC,
					  bool filter, bool merge,
					  double maxDT, double maxS) const
  {
    // Get cuts.

    if(maxDT == 0.)
      maxDT = fMaxDT;
    if(maxS == 0.)
      maxS = fMaxS;  

    // Get geometry service.

    art::ServiceHandle<geo::Geometry> geom;

    // Print diagnostic information.

    update();

    // Get time offsets.

    std::vector< std::vector<std::vector<double> > > timeOffset;
    fillTimeOffset(timeOffset);

    // First make result vector is empty.

    spts.erase(spts.begin(), spts.end());

    // Statistics.

    int n2 = 0;  // Number of two-hit space points.
    int n3 = 0;  // Number of three-hit space points.
    int n2filt = 0;  // Number of two-hit space points after filtering/merging.
    int n3filt = 0;  // Number of three-hit space pointe after filtering/merging.

    // If useMC is true, verify that channels are sorted by channel number.

    if(useMC) {

      unsigned int nsc = simchans.size();
      for(unsigned int isc = 0; isc < nsc; ++isc) {
	const sim::SimChannel* psc = simchans[isc];
	if(psc != 0 && isc != psc->Channel())
	  throw cet::exception("SpacePointService") << "MC channels not sorted.\n";
      }
    }

    // Sort hits into maps indexed by [cryostat][tpc][plane][wire].
    // If using mc information, also generate maps of sim::IDEs and mc 
    // position indexed by hit.

    std::vector< std::vector<std::vector<std::map<unsigned int, art::Ptr<recob::Hit> > > > > hitmap;
    fHitMCMap.clear();

    unsigned int ncstat = geom->Ncryostats();
    hitmap.resize(ncstat);
    for(unsigned int cstat = 0; cstat < ncstat; ++cstat){
      unsigned int ntpc = geom->Cryostat(cstat).NTPC();
      hitmap[cstat].resize(ntpc);
      for(unsigned int tpc = 0; tpc < ntpc; ++tpc) {
	int nplane = geom->Cryostat(cstat).TPC(tpc).Nplanes();
	hitmap[cstat][tpc].resize(nplane);
      }
    }
  
    for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
	ihit != hits.end(); ++ihit) {
      const art::Ptr<recob::Hit>& phit = *ihit;
      geo::View_t view = phit->View();
      if((view == geo::kU && fEnableU) ||
	 (view == geo::kV && fEnableV) ||
	 (view == geo::kW && fEnableW)) {
      
	unsigned short channel = phit->Channel();
	unsigned int tpc, plane, wire, cstat;
	geom->ChannelToWire(channel, cstat, tpc, plane, wire);
	hitmap[cstat][tpc][plane][wire] = phit;
      }
    }

    // Fill mc information, including IDEs and closest neighbors
    // of each hit.
    ///\todo Why are we still checking on whether this is MC or not?  
    ///\todo Such checks should not be in reconstruction code.
    if(useMC) {

      // First loop over hits and fill track ids and mc position.
      for(unsigned int cstat = 0; cstat < ncstat; ++cstat){
	for(unsigned int tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc) {
	  int nplane = geom->Cryostat(cstat).TPC(tpc).Nplanes();
	  for(int plane = 0; plane < nplane; ++plane) {
	    for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator ihit = hitmap[cstat][tpc][plane].begin();
		ihit != hitmap[cstat][tpc][plane].end(); ++ihit) {
	      const art::Ptr<recob::Hit>& phit = ihit->second;
	      const recob::Hit& hit = *phit;
	      HitMCInfo& mcinfo = fHitMCMap[&hit];   // Default HitMCInfo.
	    
	      // Fill default nearest neighbor information (i.e. none).
	    
	      mcinfo.pchit.resize(nplane, 0);
	      mcinfo.dist2.resize(nplane, 1.e20);
	    
	      // Get sim::IDEs for this hit.

	      std::vector<sim::IDE> ides;
	      cheat::BackTracker::HitToSimIDEs(*simchans[hit.Channel()], phit, ides);
	    
	      // Get sorted track ids. for this hit.
	    
	      mcinfo.trackIDs.reserve(ides.size());
	      for(std::vector<sim::IDE>::const_iterator i = ides.begin();
		  i != ides.end(); ++i)
		mcinfo.trackIDs.push_back(i->trackID);
	      sort(mcinfo.trackIDs.begin(), mcinfo.trackIDs.end());
	    
	      // Get position of ionization for this hit.
	    
	      mcinfo.xyz = cheat::BackTracker::SimIDEsToXYZ(ides);
	    } // end loop over ihit
	  }// end loop oer planes
	}// end loop over TPCs
      }// end loop over cryostats

      // Loop over hits again and fill nearest neighbor information for real.
      for(unsigned int cstat = 0; cstat < ncstat; ++cstat){
	for(unsigned int tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc) {
	  int nplane = geom->Cryostat(cstat).TPC(tpc).Nplanes();
	  for(int plane = 0; plane < nplane; ++plane) {
	    for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator ihit = hitmap[cstat][tpc][plane].begin();
		ihit != hitmap[cstat][tpc][plane].end(); ++ihit) {
	      const art::Ptr<recob::Hit>& phit = ihit->second;
	      const recob::Hit& hit = *phit;
	      HitMCInfo& mcinfo = fHitMCMap[&hit];
	    
	      // Fill nearest neighbor information for this hit.
	    
	      for(int plane2 = 0; plane2 < nplane; ++plane2) {
		for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator jhit = hitmap[cstat][tpc][plane2].begin();
		    jhit != hitmap[cstat][tpc][plane2].end(); ++jhit) {
		  const art::Ptr<recob::Hit>& phit2 = jhit->second;
		  const recob::Hit& hit2 = *phit2;
		  const HitMCInfo& mcinfo2 = fHitMCMap[&hit2];
		
		  assert(mcinfo.xyz.size() == 3);
		  assert(mcinfo2.xyz.size() == 3);
		  double dx = mcinfo.xyz[0] - mcinfo2.xyz[0];
		  double dy = mcinfo.xyz[1] - mcinfo2.xyz[1];
		  double dz = mcinfo.xyz[2] - mcinfo2.xyz[2];
		  double dist2 = dx*dx + dy*dy + dz*dz;
		  if(dist2 < mcinfo.dist2[plane2]) {
		    mcinfo.dist2[plane2] = dist2;
		    mcinfo.pchit[plane2] = &hit2;
		  }
		}// end loop over jhit
	      }// end loop over plane2
	    }// end loop over ihit
	  }// end loop over plane
	}// end loop over tpc
      }// end loop over cryostats
    }// end if MC

    mf::LogDebug debug("SpacePointService");
    debug << "Total hits = " << hits.size() << "\n\n";

    for(unsigned int cstat = 0; cstat < ncstat; ++cstat){
      for(unsigned int tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc) {
	int nplane = hitmap[cstat][tpc].size();
	for(int plane = 0; plane < nplane; ++plane) {
	  debug << "TPC, Plane: " << tpc << ", " << plane 
		<< ", hits = " << hitmap[cstat][tpc][plane].size() << "\n";
	}
      }
    }// end loop over cryostats

    // Loop over TPCs.
    for(unsigned int cstat = 0; cstat < ncstat; ++cstat){
      for(unsigned int tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc) {

	// Make empty multimap from hit pointer on most-populated plane to space points 
	// that include that hit (used for filtering and merging).

	typedef const recob::Hit* sptkey_type;
	std::multimap<sptkey_type, recob::SpacePoint> sptmap;
	std::set<sptkey_type> sptkeys;              // Keys of multimap.

	// Sort maps in increasing order of number of hits.
	// This is so that we can do the outer loops over hits 
	// over the views with fewer hits.
      
	int nplane = hitmap[cstat][tpc].size();
	std::vector<int> index(nplane);

	for(int i=0; i<nplane; ++i)
	  index[i] = i;

	for(int i=0; i<nplane-1; ++i) {
	  for(int j=i+1; j<nplane; ++j) {
	    if(hitmap[cstat][tpc][index[i]].size() > hitmap[cstat][tpc][index[j]].size()) {
	      int temp = index[i];
	      index[i] = index[j];
	      index[j] = temp;
	    }
	  }
	}// end loop over i

	// If two-view space points are allowed, make a double loop
	// over hits and produce space points for compatible hit-pairs.
      
	if(fMinViews <= 2) {

	  // Loop over pairs of views.

	  for(int i=0; i<nplane-1; ++i) {
	    unsigned int plane1 = index[i];
	  
	    for(int j=i+1; j<nplane; ++j) {
	      unsigned int plane2 = index[j];
	    
	      // Get angle, pitch, and offset of plane2 wires.
	    
	      const geo::WireGeo& wgeo2 = geom->Cryostat(cstat).TPC(tpc).Plane(plane2).Wire(0);
	      double hl2 = wgeo2.HalfL();
	      double xyz21[3];
	      double xyz22[3];
	      wgeo2.GetCenter(xyz21, -hl2);
	      wgeo2.GetCenter(xyz22, hl2);
	      double s2 = (xyz22[1] - xyz21[1]) / (2.*hl2);
	      double c2 = (xyz22[2] - xyz21[2]) / (2.*hl2);
	      double dist2 = -xyz21[1] * c2 + xyz21[2] * s2;
	      double pitch2 = geom->WirePitch(0, 1, plane2, tpc, cstat);
	    
	      assert(hitmap[cstat][tpc][plane1].size() <= hitmap[cstat][tpc][plane2].size());

	      // Loop over pairs of hits.
	    
	      art::PtrVector<recob::Hit> hitvec;
	      hitvec.reserve(2);
	    
	      for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator
		    ihit1 = hitmap[cstat][tpc][plane1].begin();
		  ihit1 != hitmap[cstat][tpc][plane1].end(); ++ihit1) {
	      
		const art::Ptr<recob::Hit>& phit1 = ihit1->second;
		unsigned short channel1 = phit1->Channel();
	      
		// Get endpoint coordinates of this wire.
	      
		unsigned int tpc1a, plane1a, wire1, cstat1a;
		const geo::WireGeo& wgeo = geom->ChannelToWire(channel1, cstat1a, tpc1a, plane1a, wire1);
		assert(cstat1a == cstat);
		assert(tpc1a == tpc);
		assert(plane1a == plane1);
		double hl1 = wgeo.HalfL();
		double xyz1[3];
		double xyz2[3];
		wgeo.GetCenter(xyz1, -hl1);
		wgeo.GetCenter(xyz2, hl1);
	      
		// Find the plane2 wire numbers corresponding to the endpoints.
	      
		double wire21 = (-xyz1[1] * c2 + xyz1[2] * s2 - dist2) / pitch2;
		double wire22 = (-xyz2[1] * c2 + xyz2[2] * s2 - dist2) / pitch2;
	      
		int wmin = std::min(wire21, wire22);
		int wmax = std::max(wire21, wire22) + 1.;
	      
		for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator 
		      ihit2 = hitmap[cstat][tpc][plane2].lower_bound(wmin);
		    ihit2 != hitmap[cstat][tpc][plane2].upper_bound(wmax); ++ihit2) {
		
		  const art::Ptr<recob::Hit>& phit2 = ihit2->second;
		
		  // Check current pair of hits for compatibility.
		  // By construction, hits should always have compatible views 
		  // and times, but may not have compatible mc information.
		
		  hitvec.clear();
		  hitvec.push_back(phit1);
		  hitvec.push_back(phit2);
		  bool ok = compatible(hitvec, timeOffset, useMC, maxDT, maxS);
		  if(ok) {
		  
		    // Add a space point.
		  
		    ++n2;
		    if(filter || merge) {
		      sptkey_type key = &*phit2;
		      std::multimap<sptkey_type, recob::SpacePoint>::iterator it = 
			sptmap.insert(std::pair<sptkey_type, recob::SpacePoint>(key, recob::SpacePoint()));
		      sptkeys.insert(key);
		      recob::SpacePoint& spt = it->second;
		      fillSpacePoint(hitvec, timeOffset, spt);
		    }
		    else {
		      spts.push_back(recob::SpacePoint());
		      fillSpacePoint(hitvec, timeOffset, spts.back());
		    }
		  }
		}
	      }
	    }
	  }
	}// end if fMinViews <= 2

	// If three-view space points are allowed, make a triple loop
	// over hits and produce space points for compatible triplets.
      
	if(nplane >= 3 && fMinViews <= 3) {
	
	  // Loop over triplets of hits.
	
	  art::PtrVector<recob::Hit> hitvec;
	  hitvec.reserve(3);
	
	  unsigned int plane1 = index[0];
	  unsigned int plane2 = index[1];
	  unsigned int plane3 = index[2];
	
	  // Get angle, pitch, and offset of plane1 wires.
	
	  const geo::WireGeo& wgeo1 = geom->Cryostat(cstat).TPC(tpc).Plane(plane1).Wire(0);
	  double hl1 = wgeo1.HalfL();
	  double xyz11[3];
	  double xyz12[3];
	  wgeo1.GetCenter(xyz11, -hl1);
	  wgeo1.GetCenter(xyz12, hl1);
	  double s1 = (xyz12[1] - xyz11[1]) / (2.*hl1);
	  double c1 = (xyz12[2] - xyz11[2]) / (2.*hl1);
	  double dist1 = -xyz11[1] * c1 + xyz11[2] * s1;
	  double pitch1 = geom->WirePitch(0, 1, plane1, tpc, cstat);

	  // Get angle, pitch, and offset of plane2 wires.
	
	  const geo::WireGeo& wgeo2 = geom->Cryostat(cstat).TPC(tpc).Plane(plane2).Wire(0);
	  double hl2 = wgeo2.HalfL();
	  double xyz21[3];
	  double xyz22[3];
	  wgeo2.GetCenter(xyz21, -hl2);
	  wgeo2.GetCenter(xyz22, hl2);
	  double s2 = (xyz22[1] - xyz21[1]) / (2.*hl2);
	  double c2 = (xyz22[2] - xyz21[2]) / (2.*hl2);
	  double dist2 = -xyz21[1] * c2 + xyz21[2] * s2;
	  double pitch2 = geom->WirePitch(0, 1, plane2, tpc, cstat);

	  // Get angle, pitch, and offset of plane3 wires.
	
	  const geo::WireGeo& wgeo3 = geom->Cryostat(cstat).TPC(tpc).Plane(plane3).Wire(0);
	  double hl3 = wgeo3.HalfL();
	  double xyz31[3];
	  double xyz32[3];
	  wgeo3.GetCenter(xyz31, -hl3);
	  wgeo3.GetCenter(xyz32, hl3);
	  double s3 = (xyz32[1] - xyz31[1]) / (2.*hl3);
	  double c3 = (xyz32[2] - xyz31[2]) / (2.*hl3);
	  double dist3 = -xyz31[1] * c3 + xyz31[2] * s3;
	  double pitch3 = geom->WirePitch(0, 1, plane3, tpc, cstat);

	  // Get sine of angle differences.
	
	  double s12 = s1 * c2 - s2 * c1;   // sin(theta1 - theta2).
	  double s23 = s2 * c3 - s3 * c2;   // sin(theta2 - theta3).
	  double s31 = s3 * c1 - s1 * c3;   // sin(theta3 - theta1).
	
	  // Loop over hits in plane1.
	
	  for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator
		ihit1 = hitmap[cstat][tpc][plane1].begin();
	      ihit1 != hitmap[cstat][tpc][plane1].end(); ++ihit1) {
	  
	    unsigned int wire1 = ihit1->first;
	    const art::Ptr<recob::Hit>& phit1 = ihit1->second;
	    unsigned short channel1 = phit1->Channel();
	  
	    // Get endpoint coordinates of this wire from plane1.
	  
	    unsigned int tpc1a, plane1a, wire1a, cstat1a;
	    const geo::WireGeo& wgeo = geom->ChannelToWire(channel1, cstat1a, tpc1a, plane1a, wire1a);
	    assert(cstat1a == cstat);
	    assert(tpc1a == tpc);
	    assert(plane1a == plane1);
	    assert(wire1a == wire1);
	    double hl1 = wgeo.HalfL();
	    double xyz1[3];
	    double xyz2[3];
	    wgeo.GetCenter(xyz1, -hl1);
	    wgeo.GetCenter(xyz2, hl1);
	  
	    // Get corrected time and oblique coordinate of first hit.
	  
	    double t1 = phit1->PeakTime() - timeOffset[cstat][tpc][plane1];
	    double u1 = wire1 * pitch1 + dist1;
	  
	    // Find the plane2 wire numbers corresponding to the endpoints.
	  
	    double wire21 = (-xyz1[1] * c2 + xyz1[2] * s2 - dist2) / pitch2;
	    double wire22 = (-xyz2[1] * c2 + xyz2[2] * s2 - dist2) / pitch2;

	    int wmin = std::min(wire21, wire22);
	    int wmax = std::max(wire21, wire22) + 1.;
	  
	    for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator 
		  ihit2 = hitmap[cstat][tpc][plane2].lower_bound(wmin);
		ihit2 != hitmap[cstat][tpc][plane2].upper_bound(wmax); ++ihit2) {

	      int wire2 = ihit2->first;
	      const art::Ptr<recob::Hit>& phit2 = ihit2->second;
	    
	      // Get corrected time of second hit.
	    
	      double t2 = phit2->PeakTime() - timeOffset[cstat][tpc][plane2];
	    
	      // Check maximum time difference with first hit.
	    
	      bool dt12ok = std::abs(t1-t2) <= maxDT;
	      if(dt12ok) {
	      
		// Test first two hits for compatibility before looping 
		// over third hit.
	      
		hitvec.clear();
		hitvec.push_back(phit1);
		hitvec.push_back(phit2);
		bool h12ok = compatible(hitvec, timeOffset, useMC, maxDT, maxS);
		if(h12ok) {
		
		  // Get oblique coordinate of second hit.
		
		  double u2 = wire2 * pitch2 + dist2;
		
		  // Predict plane3 oblique coordinate and wire number.
		
		  double u3pred = (-u1*s23 - u2*s31) / s12;
		  double w3pred = (u3pred - dist3) / pitch3;
		  double w3delta = std::abs(maxS / (s12 * pitch3));
		  int w3min = std::ceil(w3pred - w3delta);
		  int w3max = std::floor(w3pred + w3delta);
		
		  for(std::map<unsigned int, art::Ptr<recob::Hit> >::const_iterator 
			ihit3 = hitmap[cstat][tpc][plane3].lower_bound(w3min);
		      ihit3 != hitmap[cstat][tpc][plane3].upper_bound(w3max); ++ihit3) {
		  
		    int wire3 = ihit3->first;
		    const art::Ptr<recob::Hit>& phit3 = ihit3->second;
		  
		    // Get corrected time of third hit.
		  
		    double t3 = phit3->PeakTime() - timeOffset[cstat][tpc][plane3];
		  
		    // Check time difference of third hit compared to first two hits.
		  
		    bool dt123ok = std::abs(t1-t3) <= maxDT && std::abs(t2-t3) <= maxDT;
		    if(dt123ok) {
		    
		      // Get oblique coordinate of third hit and check spatial separation.
		    
		      double u3 = wire3 * pitch3 + dist3;
		      double S = s23 * u1 + s31 * u2 + s12 * u3;
		      bool sok = std::abs(S) <= maxS;
		      if(sok) {
		      
			// Test triplet for compatibility.
		      
			hitvec.clear();
			hitvec.push_back(phit1);
			hitvec.push_back(phit2);
			hitvec.push_back(phit3);
			bool h123ok = compatible(hitvec, timeOffset, useMC, maxDT, maxS);
			if(h123ok) {
			
			  // Add a space point.
			
			  ++n3;
			  if(filter || merge) {
			    sptkey_type key = &*phit3;
			    std::multimap<sptkey_type, recob::SpacePoint>::iterator it = 
			      sptmap.insert(std::pair<sptkey_type, recob::SpacePoint>(key, recob::SpacePoint()));
			    sptkeys.insert(key);
			    recob::SpacePoint& spt = it->second;
			    fillSpacePoint(hitvec, timeOffset, spt);
			  }
			  else {
			    spts.push_back(recob::SpacePoint());
			    fillSpacePoint(hitvec, timeOffset, spts.back());
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}// end if fMinViews <= 3
          
	// Do Filtering.
      
	if(filter) {
	
	  // Transfer (some) space points from sptmap to spts.
	
	  spts.reserve(spts.size() + sptkeys.size());
	
	  // Loop over keys of space point map.
	  // Space points that have the same key are candidates for filtering.
	
	  for(std::set<sptkey_type>::const_iterator i = sptkeys.begin();
	      i != sptkeys.end(); ++i) {
	    sptkey_type key = *i;
	  
	    // Loop over space points corresponding to the current key.
	    // Choose the single best space point from among this group.
	  
	    double best_chisq = 0.;
	    const recob::SpacePoint* best_spt = 0;
	  
	    for(std::multimap<sptkey_type, recob::SpacePoint>::const_iterator j = sptmap.lower_bound(key);
		j != sptmap.upper_bound(key); ++j) {
	      const recob::SpacePoint& spt = j->second;
	      if(best_spt == 0 || spt.Chisq() < best_chisq) {
		best_spt = &spt;
		best_chisq = spt.Chisq();
	      }
	    }
	  
	    // Transfer best filtered space point to result vector.
	  
	    assert(best_spt != 0);
	    if(best_spt != 0) {
	      spts.push_back(*best_spt);
	      if(fMinViews <= 2)
		++n2filt;
	      else
		++n3filt;
	    }
	  }
	}// end if filtering
      
	// Do merging.
      
	else if(merge) {
	
	  // Transfer merged space points from sptmap to spts.
	
	  spts.reserve(spts.size() + sptkeys.size());
	
	  // Loop over keys of space point map.
	  // Space points that have the same key are candidates for merging.
	
	  for(std::set<sptkey_type>::const_iterator i = sptkeys.begin();
	      i != sptkeys.end(); ++i) {
	    sptkey_type key = *i;
	  
	    // Loop over space points corresponding to the current key.
	    // Make a collection of hits that is the union of the hits
	    // from each candidate space point.
	  
	    art::PtrVector<recob::Hit> merged_hits;
	    for(std::multimap<sptkey_type, recob::SpacePoint>::const_iterator j = sptmap.lower_bound(key);
		j != sptmap.upper_bound(key); ++j) {
	      const recob::SpacePoint& spt = j->second;
	    
	      // Loop over hits from this space points.
	      // Add each hit to the collection of all hits.
	    
	      const art::PtrVector<recob::Hit>& spt_hits = spt.AllHits();
	      for(art::PtrVector<recob::Hit>::const_iterator k = spt_hits.begin();
		  k != spt_hits.end(); ++k) {
		const art::Ptr<recob::Hit>& hit = *k;
		merged_hits.push_back(hit);
	      }
	    }
	  
	    // Remove duplicates.
	  
	    std::sort(merged_hits.begin(), merged_hits.end());
	    art::PtrVector<recob::Hit>::iterator it = 
	      std::unique(merged_hits.begin(), merged_hits.end());
	    merged_hits.erase(it, merged_hits.end());
	  
	    // Construct a complex space points using merged hits.
	  
	    spts.push_back(recob::SpacePoint());
	    fillComplexSpacePoint(merged_hits, timeOffset, spts.back());
	  
	    if(fMinViews <= 2)
	      ++n2filt;
	    else
	      ++n3filt;
	  }
	}// end if merging
	else {
	  n2filt = n2;
	  n3filt = n3;
	}
      }// end loop over tpcs
    }// end loop over cryostats
  
    debug << "\n2-hit space points = " << n2 << "\n"
	  << "3-hit space points = " << n3 << "\n"
	  << "2-hit filtered/merged space points = " << n2filt << "\n"
	  << "3-hit filtered/merged space points = " << n3filt;
  }

} // end namespace trkf
