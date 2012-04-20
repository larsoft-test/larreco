////////////////////////////////////////////////////////////////////////
// Class:       TrackKalmanCheater
// Module Type: producer
// File:        TrackKalmanCheater.cxx
//
// Generated at Wed Mar 28 13:43:51 2012 by Herbert Greenlee using artmod
// from art v1_00_11.
////////////////////////////////////////////////////////////////////////

#include <vector>
#include "TrackFinder/TrackKalmanCheater.h"
#include "TrackFinder/KHitContainerWireX.h"
#include "TrackFinder/SurfYZPlane.h"
#include "TrackFinder/PropYZPlane.h"
#include "TrackFinder/KalmanFilterService.h"
#include "Geometry/geo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "Simulation/SimListUtils.h"
#include "MCCheater/BackTracker.h"
#include "Utilities/AssociationUtil.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


/// Constructor.
///
/// Arguments:
///
/// p - Fcl parameters.
///
trkf::TrackKalmanCheater::TrackKalmanCheater(fhicl::ParameterSet const & pset) :
  fUseClusterHits(false),
  prop(new PropYZPlane(10.)),
  fNumEvent(0),
  fNumTrack(0)
{
  reconfigure(pset);
  produces<std::vector<recob::Track> >();
  produces<art::Assns<recob::Track, recob::Hit> >();

  // Report.

  mf::LogInfo("TrackKalmanCheater") 
    << "TrackKalmanCheater configured with the following parameters:\n"
    << "  UseClusterHits = " << fUseClusterHits << "\n"
    << "  GenModuleLabel = " << fGenModuleLabel << "\n"
    << "  HitModuleLabel = " << fHitModuleLabel << "\n"
    << "  ClusterModuleLabel = " << fClusterModuleLabel << "\n"
    << "  G4ModuleLabel = " << fG4ModuleLabel;
}

/// Destructor.
trkf::TrackKalmanCheater::~TrackKalmanCheater()
{}

/// Reconfigure method.
///
/// Arguments:
///
/// p - Fcl parameter set.
///
void trkf::TrackKalmanCheater::reconfigure(fhicl::ParameterSet const & pset)
{
  fUseClusterHits = pset.get<bool>("UseClusterHits");
  fGenModuleLabel = pset.get<std::string>("GenModuleLabel");
  fHitModuleLabel = pset.get<std::string>("HitModuleLabel");
  fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
  fG4ModuleLabel = pset.get<std::string>("G4ModuleLabel");
}

/// Begin job method.
void trkf::TrackKalmanCheater::beginJob()
{}

/// Produce method.
///
/// Arguments:
///
/// e - Art event.
///
/// This method extracts Hit from the event and produces and adds
/// Track objects.
///
void trkf::TrackKalmanCheater::produce(art::Event & evt)
{
  ++fNumEvent;

  // Get Services.

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<trkf::KalmanFilterService> kf;

  // Get SimChannels.
  // Make a vector where each channel in the detector is an entry
  std::vector<const sim::SimChannel*> simchanh;
  std::vector<const sim::SimChannel*> simchanv(geom->Nchannels(),0);
  evt.getView(fG4ModuleLabel, simchanh);
  for(size_t i = 0; i < simchanh.size(); ++i)
    simchanv[simchanh[i]->Channel()] = simchanh[i];

  // Get Hits.

  art::PtrVector<recob::Hit> hits;

  if(fUseClusterHits) {

    // Get clusters.

    art::Handle< std::vector<recob::Cluster> > clusterh;
    evt.getByLabel(fClusterModuleLabel, clusterh);

    // Get hits from all clusters.

    if(clusterh.isValid()) {
      int nclus = clusterh->size();

      for(int i = 0; i < nclus; ++i) {
	art::Ptr<recob::Cluster> pclus(clusterh, i);
	art::PtrVector<recob::Hit> clushits = pclus->Hits();
	int nhits = clushits.size();
	hits.reserve(hits.size() + nhits);

	for(art::PtrVector<recob::Hit>::const_iterator ihit = clushits.begin();
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

  // Sort hits into separate PtrVectors based on track id.

  std::map<int, art::PtrVector<recob::Hit> > hitmap;

  // Loop over hits.

  for(art::PtrVector<recob::Hit>::const_iterator ihit = hits.begin();
	ihit != hits.end(); ++ihit) {
    const recob::Hit& hit = **ihit;

    // Get track ids for this hit.

    unsigned int channel = hit.Channel();
    assert(channel < simchanv.size());
    const sim::SimChannel* simchan = simchanv[channel];
    std::vector<cheat::TrackIDE> tids = cheat::BackTracker::HitToTrackID(*simchan, *ihit);

    // Loop over track ids.

    for(std::vector<cheat::TrackIDE>::const_iterator itid = tids.begin();
	itid != tids.end(); ++itid) {
      int trackID = itid->trackID;

      // Add hit to PtrVector corresponding to this track id.

      hitmap[trackID].push_back(*ihit);
    }
  }

  // Extract geant mc particles.

  sim::ParticleList plist = sim::SimListUtils::GetParticleList(evt, fG4ModuleLabel);

  // Loop over geant particles.

  {
    mf::LogDebug log("TrackKalmanCheater");
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

	int trackid = part->TrackId();
	int stat = part->StatusCode();
	int nhit = 0;
	if(hitmap.count(trackid) != 0)
	  nhit = hitmap[trackid].size();
	log << "Trackid=" << trackid 
	    << ", pdgid=" << pdg 
	    << ", status=" << stat
	    << ", NHit=" << nhit << "\n";
	const TLorentzVector& pos = part->Position();
	double x = pos.X();
	double y = pos.Y();
	double z = pos.Z();
	const TLorentzVector& mom = part->Momentum();
	double px = mom.Px();
	double py = mom.Py();
	double pz = mom.Pz();
	double p = std::sqrt(px*px + py*py + pz*pz);
	log << "  x = " << pos.X()
	    << ", y = " << pos.Y()
	    << ", z = " << pos.Z() << "\n"
	    << "  px= " << mom.Px()
	    << ", py= " << mom.Py()
	    << ", pz= " << mom.Pz() << "\n";

	if(nhit > 0 && pz != 0.) {
	  const art::PtrVector<recob::Hit>& trackhits = hitmap[trackid];
	  assert(trackhits.size() > 0);

	  // Make a seed track (KHitsTrack).

	  boost::shared_ptr<const Surface> psurf(new SurfYZPlane(y, z, 0.));
	  TrackVector vec(5);
	  vec(0) = x;
	  vec(1) = 0.;
	  vec(2) = px / pz;
	  vec(3) = py / pz;
	  vec(4) = 1. / p;
	  TrackError err(5);
	  err.clear();
	  //err(0,0) = 1000.;
	  //err(1,1) = 1000.;
	  //err(2,2) = 10.;
	  //err(3,3) = 10.;
	  err(0,0) = 1.;
	  err(1,1) = 1.;
	  err(2,2) = 0.005;
	  err(3,3) = 0.005;
	  err(4,4) = 1.;
	  Surface::TrackDirection dir = (pz > 0. ? Surface::FORWARD : Surface::BACKWARD);
	  KHitsTrack trh(KETrack(psurf, vec, err, dir, pdg));
	  trh.setStat(KFitTrack::UNKNOWN);

	  // Fill KHitContainer with hits.

	  KHitContainerWireX cont;
	  cont.fill(trackhits, 0);

	  // Build track.

	  kf->buildTrack(trh, prop, Propagator::FORWARD, cont, 100.);
 	}
      }
    }
  }
}

/// End job method.
void trkf::TrackKalmanCheater::endJob()
{
  mf::LogInfo("TrackKalmanCheater") 
    << "TrackKalmanCheater statistics:\n"
    << "  Number of events = " << fNumEvent << "\n"
    << "  Number of tracks created = " << fNumTrack;
}

