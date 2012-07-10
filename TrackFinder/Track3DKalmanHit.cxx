/////////////////////////////////////////////////////////////////////
// Class:       Track3DKalmanHit
// Module Type: producer
// File:        Track3DKalmanHit.cxx
////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <vector>
#include <deque>
#include "TrackFinder/Track3DKalmanHit.h"
#include "TrackFinder/KHitContainerWireX.h"
#include "TrackFinder/SurfYZPlane.h"
#include "TrackFinder/PropYZPlane.h"
#include "TrackFinder/KHit.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"
#include "Utilities/AssociationUtil.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"

// Local functions.

namespace {

  //----------------------------------------------------------------------------
  // Fill a collection of hits used by space points.
  //
  // Arguments:
  //
  // spts - Space point collection.
  // hits - Hit collection to fill.
  //
  void SpacePointsToHits(const std::vector<recob::SpacePoint>& spts,
			 art::PtrVector<recob::Hit>& hits, 
			 trkf::SpacePointAlg const& spalg)
  {
    hits.clear();
    std::set<const recob::Hit*> used_hits;

    // Loop over space points.

    for(std::vector<recob::SpacePoint>::const_iterator ispt = spts.begin();
	ispt != spts.end(); ++ispt) {
      const recob::SpacePoint& spt = *ispt;

      // Loop over hits in space point.

      const art::PtrVector<recob::Hit> spthits = spalg.getAssociatedHits(spt);
      for(art::PtrVector<recob::Hit>::const_iterator ihit = spthits.begin();
	  ihit != spthits.end(); ++ihit) {
	const art::Ptr<recob::Hit> phit = *ihit;
	const recob::Hit& hit = *phit;

	// Check if hit should be added to collection.

	if(used_hits.count(&hit) == 0) {
	  used_hits.insert(&hit);
	  hits.push_back(phit);
	}
      }
    }
  }

  //----------------------------------------------------------------------------
  // Filter a collection of hits.
  //
  // Arguments:
  //
  // hits      - Hit collection from which hits should be removed.
  // used_hits - Hits to remove.
  //
  void FilterHits(art::PtrVector<recob::Hit>& hits,
		  art::PtrVector<recob::Hit>& used_hits)
  {
    if(used_hits.size() > 0) {

      // Make sure both hit collections are sorted.

      std::stable_sort(hits.begin(), hits.end());
      std::stable_sort(used_hits.begin(), used_hits.end());

      // Do set difference operation.

      art::PtrVector<recob::Hit>::iterator it =
	std::set_difference(hits.begin(), hits.end(),
			    used_hits.begin(), used_hits.end(),
			    hits.begin());

      // Truncate hit collection.

      hits.erase(it, hits.end());
    }
  }

  //----------------------------------------------------------------------------
  // Get preferred plane (find plane that has the most hits).
  //
  // Arguments:
  //
  // hits - Hit collection.
  //
  unsigned int GetPreferredPlane(const art::PtrVector<recob::Hit>& hits)
  {
    // Get geometry service.

    art::ServiceHandle<geo::Geometry> geom;

    // Count hits in each plane.

    std::vector<unsigned int> planehits(3, 0);

    // Loop over hits.

    for(art::PtrVector<recob::Hit>::const_iterator ih = hits.begin();
	ih != hits.end(); ++ih) {
      const recob::Hit& hit = **ih;
      unsigned int channel = hit.Channel();
      unsigned int cstat, tpc, plane, wire;
      geom->ChannelToWire(channel, cstat, tpc, plane, wire);
      assert(plane >= 0 && plane < planehits.size());
      ++planehits[plane];
    }

    // Figure out which plane has the most hits.

    unsigned int prefplane = 0;
    for(unsigned int i=0; i<planehits.size(); ++i) {
      if(planehits[i] > planehits[prefplane])
	prefplane = i;
    }
    return prefplane;
  }
}


//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// p - Fcl parameters.
///
trkf::Track3DKalmanHit::Track3DKalmanHit(fhicl::ParameterSet const & pset) :
  fHist(false),
  fUseClusterHits(false),
  fMaxTcut(0.),
  fMinSeedHits(0.),
  fMaxSeedChiDF(0.),
  fMinSeedSlope(0.),
  fKFAlg(pset.get<fhicl::ParameterSet>("KalmanFilterAlg")),
  fSeedFinderAlg(pset.get<fhicl::ParameterSet>("SeedFinderAlg")),
  fSpacePointAlg(pset.get<fhicl::ParameterSet>("SpacePointAlg")),
  fProp(0),
  fHIncChisq(0),
  fHPull(0),
  fNumEvent(0),
  fNumTrack(0)
{
  reconfigure(pset);
  produces<std::vector<recob::Track> >();
  produces<std::vector<recob::SpacePoint>            >();
  produces<art::Assns<recob::Track, recob::Hit> >();
  produces<art::Assns<recob::Track, recob::SpacePoint> >();
  produces<art::Assns<recob::SpacePoint, recob::Hit> >();

  // Report.

  mf::LogInfo("Track3DKalmanHit") 
    << "Track3DKalmanHit configured with the following parameters:\n"
    << "  UseClusterHits = " << fUseClusterHits << "\n"
    << "  HitModuleLabel = " << fHitModuleLabel << "\n"
    << "  ClusterModuleLabel = " << fClusterModuleLabel;
}

//----------------------------------------------------------------------------
/// Destructor.
trkf::Track3DKalmanHit::~Track3DKalmanHit()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// p - Fcl parameter set.
///
void trkf::Track3DKalmanHit::reconfigure(fhicl::ParameterSet const & pset)
{
  fHist = pset.get<bool>("Hist");
  fKFAlg.reconfigure(pset.get<fhicl::ParameterSet>("KalmanFilterAlg"));
  fSeedFinderAlg.reconfigure(pset.get<fhicl::ParameterSet>("SeedFinderAlg"));
  fSpacePointAlg.reconfigure(pset.get<fhicl::ParameterSet>("SpacePointAlg"));
  fUseClusterHits = pset.get<bool>("UseClusterHits");
  fHitModuleLabel = pset.get<std::string>("HitModuleLabel");
  fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
  fMaxTcut = pset.get<double>("MaxTcut");
  fMinSeedHits = pset.get<double>("MinSeedHits");
  fMaxSeedChiDF = pset.get<double>("MaxSeedChiDF");
  fMinSeedSlope = pset.get<double>("MinSeedSlope");
  if(fProp != 0)
    delete fProp;
  fProp = new PropYZPlane(fMaxTcut);
}

//----------------------------------------------------------------------------
/// Begin job method.
void trkf::Track3DKalmanHit::beginJob()
{
  if(fHist) {

    // Book histograms.

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir("hitkalman", "Track3DKalmanHit histograms");

    fHIncChisq = dir.make<TH1F>("IncChisq", "Incremental Chisquare", 100, 0., 20.);
    fHPull = dir.make<TH1F>("Pull", "Hit Pull", 100, -10., 10.);
  }
}

//----------------------------------------------------------------------------
/// Produce method.
///
/// Arguments:
///
/// e - Art event.
///
/// This method extracts Hit from the event and produces and adds
/// Track objects.
///
void trkf::Track3DKalmanHit::produce(art::Event & evt)
{
  ++fNumEvent;

  // Make a collection of tracks, plus associations, that will
  // eventually be inserted into the event.

  std::auto_ptr<std::vector<recob::Track> > tracks(new std::vector<recob::Track>);
  std::auto_ptr< art::Assns<recob::Track, recob::Hit> > th_assn(new art::Assns<recob::Track, recob::Hit>);
  std::auto_ptr< art::Assns<recob::Track, recob::SpacePoint> > tsp_assn(new art::Assns<recob::Track, recob::SpacePoint>);

  // Make a collection of space points, plus associations, that will
  // be inserted into the event.

  std::auto_ptr<std::vector<recob::SpacePoint> > spts(new std::vector<recob::SpacePoint>);
  std::auto_ptr< art::Assns<recob::SpacePoint, recob::Hit> > sph_assn(new art::Assns<recob::SpacePoint, recob::Hit>);

  // Make a collection of KGTracks where we will save our results.

  std::deque<KGTrack> kalman_tracks;

  // Get Services.

  art::ServiceHandle<geo::Geometry> geom;

  // Reset space point algorithm.

  fSpacePointAlg.clearHitMap();

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

  // The hit collection "hits" (just filled), initially containing all
  // hits, represents hits available for making tracks.  Now we will
  // fill a second hit collection called "seederhits", also initially
  // containing all hits, which will represent hits available for
  // making track seeds.  These collections are not necessarily the
  // same, since hits that are not suitable for seeds may still be
  // suitable for tracks.

  art::PtrVector<recob::Hit> seederhits = hits;

  // Start of loop.

  bool done = false;
  while(!done) {

    // Use remaining hits to make space points using the seed finder.

    std::vector<recob::SpacePoint> seedspts =
      fSeedFinderAlg.GetSpacePointsFromHitVector(seederhits);
    std::vector<std::vector<recob::SpacePoint> > seedsptvecs;
    std::vector<recob::Seed*> seeds;
    if(seedspts.size() > 0)
      seeds = fSeedFinderAlg.FindSeeds(seedspts, seedsptvecs);
    if(seeds.size() == 0 || !seeds.front()->IsValid()) {

      // Quit loop if we didn't find any new seeds.

      done = true;
      break;
    }
    else {

      // Found a seed.

      mf::LogDebug log("Track3DKalmanHit");

      // This method should find just one seed.  In any case, we only
      // look at the first seed found.

      // We don't actually look at the seed object, but only the space
      // points.

      //const recob::Seed& seed = *(seeds.front());
      const std::vector<recob::SpacePoint>& seedsptvec = seedsptvecs.front();

      // Extract hits used by this seed.

      art::PtrVector<recob::Hit> seedhits;
      SpacePointsToHits(seedsptvec, seedhits, fSeedFinderAlg.GetSpacePointAlg());

      // Filter hits used by seed from hits available to make future seeds.
      // No matter what, we will never use these hits for another seed.

      FilterHits(seederhits, seedhits);

      // Convert seed into initial KTrack on a z-plane surface.

      //double xyz[3];
      //double dir[3];
      //double err[3];   // Dummy.
      //seed.GetPoint(xyz, err);
      //seed.GetDirection(dir, err);


      double x0 = seedsptvec.front().XYZ()[0];
      double y0 = seedsptvec.front().XYZ()[1];
      double z0 = seedsptvec.front().XYZ()[2];

      double x1 = seedsptvec.back().XYZ()[0];
      double y1 = seedsptvec.back().XYZ()[1];
      double z1 = seedsptvec.back().XYZ()[2];

      double dx = x1 - x0;
      double dy = y1 - y0;
      double dz = z1 - z0;

      double x = (z0 < z1 ? x0 : x1);
      double y = (z0 < z1 ? y0 : y1);
      double z = (z0 < z1 ? z0 : z1);

      boost::shared_ptr<const Surface> psurf(new SurfYZPlane(y, z, 0.));
      TrackVector vec(5);
      vec(0) = x;
      vec(1) = 0.;
      vec(2) = dx / dz;
      vec(3) = dy / dz;
      vec(4) = 2.0;

      // Since we have no information, assume pz>0.

      Surface::TrackDirection trkdir = Surface::FORWARD;

      // Also assume muon.

      int pdg = 13;

      // Make KTrack object.

      KTrack trk(psurf, vec, trkdir, pdg);

      log << "Seed found with " << seedsptvec.size() 
	  << " space points and " << seedhits.size() << " hits.\n"
	  << "(x,y,z) = " << x << ", " << y << ", " << z << "\n"
	  << "(dx,dy,dz) = " << dx << ", " << dy << ", " << dz << "\n"
	  << "(x1, y1, z1)) = "
	  << seedsptvec.front().XYZ()[0] << ", "
	  << seedsptvec.front().XYZ()[1] << ", "
	  << seedsptvec.front().XYZ()[2] << "\n"
	  << "(xn, yn, zn)) = "
	  << seedsptvec.back().XYZ()[0] << ", "
	  << seedsptvec.back().XYZ()[1] << ", "
	  << seedsptvec.back().XYZ()[2] << "\n";

      // Cut on the seed slope.

      if(std::abs(dx) >= fMinSeedSlope * std::abs(dz)) {

	// Fill hit container with current seed hits.

	KHitContainerWireX seedcont;
	seedcont.fill(seedhits, -1);

	// Set the preferred plane to be the one with the most hits.

	unsigned int prefplane = GetPreferredPlane(seedhits);
	fKFAlg.setPlane(prefplane);
	log << "Preferred plane = " << prefplane << "\n";

	// Build and smooth seed track.

	KGTrack trg0;
	bool ok = fKFAlg.buildTrack(trk, trg0, fProp, Propagator::FORWARD, seedcont);
	if(ok) {
	  KGTrack trg1;
	  ok = fKFAlg.smoothTrack(trg0, &trg1, fProp);
	  if(ok) {

	    // Now we have the seed track in the form of a KGTrack.
	    // Do additional quality cuts.

	    size_t n = trg1.numHits();
	    ok = (n >= fMinSeedHits &&
		  trg0.startTrack().getChisq() <= n * fMaxSeedChiDF &&
		  trg0.endTrack().getChisq() <= n * fMaxSeedChiDF &&
		  trg1.startTrack().getChisq() <= n * fMaxSeedChiDF &&
		  trg1.endTrack().getChisq() <= n * fMaxSeedChiDF);
	    double mom0[3];
	    double mom1[3];
	    trg0.startTrack().getMomentum(mom0);
	    trg0.endTrack().getMomentum(mom1);
	    double dxdz0 = mom0[0] / mom0[2];
	    double dxdz1 = mom1[0] / mom1[2];
	    ok = ok && (std::abs(dxdz0) > fMinSeedSlope &&
			std::abs(dxdz1) > fMinSeedSlope);
	    if(ok) {

	      // Make a copy of the original hit collection of all
	      // available track hits.

	      art::PtrVector<recob::Hit> trackhits = hits;

	      // Do an extend + smooth loop here.

	      int niter = 6;

	      for(int ix = 0; ok && ix < niter; ++ix) {

		// Fill a collection of hits from the last good track
		// (initially the seed track).

		art::PtrVector<recob::Hit> goodhits;
		trg1.fillHits(goodhits);

		// Filter hits already on the track out of the available hits.

		FilterHits(trackhits, goodhits);

		// Fill hit container using filtered hits.

		KHitContainerWireX trackcont;
		trackcont.fill(trackhits, -1);

		// Extend the track.  It is not an error for the
		// extend operation to fail, meaning that no new hits
		// were added.

		fKFAlg.extendTrack(trg1, fProp, trackcont);

		// Smooth the extended track, and make a new
		// unidirectionally fit track in the opposite
		// direction.

		KGTrack trg2;
		ok = fKFAlg.smoothTrack(trg1, &trg2, fProp);
		if(ok) {
		  KETrack tremom;
		  bool pok = fKFAlg.fitMomentum(trg1, fProp, tremom);
		  if(pok)
		    fKFAlg.updateMomentum(tremom, fProp, trg2);
		  trg1 = trg2;
		}
	      }

	      // Do a final smooth.

	      if(ok) {
		ok = fKFAlg.smoothTrack(trg1, 0, fProp);
		if(ok) {
		  KETrack tremom;
		  bool pok = fKFAlg.fitMomentum(trg1, fProp, tremom);
		  if(pok)
		    fKFAlg.updateMomentum(tremom, fProp, trg1);

		  // Save this track.

		  ++fNumTrack;
		  kalman_tracks.push_back(trg1);

		  // Filter hits on the saved track out of the hits
		  // available for making either tracks or track
		  // seeds.

		  art::PtrVector<recob::Hit> goodhits;
		  trg1.fillHits(goodhits);
		  FilterHits(hits, goodhits);
		  FilterHits(seederhits, goodhits);
		}
	      }
	    }
	  }
	}
	if(ok) {
	  log << "Find track succeeded.\n";
	}
	else
	  log << "Find track failed.\n";
      }
    }
  }

  // Fill histograms.

  // First loop over tracks.

  for(std::deque<KGTrack>::const_iterator k = kalman_tracks.begin();
      k != kalman_tracks.end(); ++k) {
    const KGTrack& trg = *k;

    // Loop over measurements in this track.

    const std::multimap<double, KHitTrack>& trackmap = trg.getTrackMap();
    for(std::multimap<double, KHitTrack>::const_iterator ih = trackmap.begin();
	ih != trackmap.end(); ++ih) {
      const KHitTrack& trh = (*ih).second;
      const boost::shared_ptr<const KHitBase>& hit = trh.getHit();
      double chisq = hit->getChisq();
      fHIncChisq->Fill(chisq);
      const KHit<1>* ph1 = dynamic_cast<const KHit<1>*>(&*hit);
      if(ph1 != 0) {
	double pull = ph1->getResVector()(0) / std::sqrt(ph1->getResError()(0, 0));
	fHPull->Fill(pull);
      }
    }
  }

  // Process Kalman filter tracks into persistent objects.

  tracks->reserve(kalman_tracks.size());
  for(std::deque<KGTrack>::const_iterator k = kalman_tracks.begin();
      k != kalman_tracks.end(); ++k) {
    const KGTrack& kalman_track = *k;

    // Add Track object to collection.

    tracks->push_back(recob::Track());
    kalman_track.fillTrack(tracks->back(), tracks->size() - 1);

    // Make Track to Hit associations.  

    art::PtrVector<recob::Hit> trhits;
    kalman_track.fillHits(hits);
    util::CreateAssn(*this, evt, *tracks, trhits, *th_assn, tracks->size()-1);

    // Make space points from this track.

    int nspt = spts->size();
    kalman_track.fillSpacePoints(*spts, fSpacePointAlg);

    // Associate newly created space points with hits.
    // Also associate track with newly created space points.

    art::PtrVector<recob::SpacePoint> sptvec;

    // Loop over newly created space points.
      
    for(unsigned int ispt = nspt; ispt < spts->size(); ++ispt) {
      const recob::SpacePoint& spt = (*spts)[ispt];
      art::ProductID sptid = getProductID<std::vector<recob::SpacePoint> >(evt);
      art::Ptr<recob::SpacePoint> sptptr(sptid, ispt, evt.productGetter(sptid));
      sptvec.push_back(sptptr);

      // Make space point to hit associations.

      const art::PtrVector<recob::Hit>& sphits = 
	fSpacePointAlg.getAssociatedHits(spt);
      util::CreateAssn(*this, evt, *spts, sphits, *sph_assn, ispt);
    }

    // Make track to space point associations.

    util::CreateAssn(*this, evt, *tracks, sptvec, *tsp_assn, tracks->size()-1);
  }

  // Add tracks and associations to event.

  evt.put(tracks);
  evt.put(spts);
  evt.put(th_assn);
  evt.put(tsp_assn);
  evt.put(sph_assn);
}

//----------------------------------------------------------------------------
/// End job method.
void trkf::Track3DKalmanHit::endJob()
{
  mf::LogInfo("Track3DKalmanHit") 
    << "Track3DKalmanHit statistics:\n"
    << "  Number of events = " << fNumEvent << "\n"
    << "  Number of tracks created = " << fNumTrack;
}
