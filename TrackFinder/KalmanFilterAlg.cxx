///////////////////////////////////////////////////////////////////////
///
/// \file   KalmanFilterAlg.cxx
///
/// \brief  Kalman Filter.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "TrackFinder/KalmanFilterAlg.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

/// Constructor.
  
trkf::KalmanFilterAlg::KalmanFilterAlg(const fhicl::ParameterSet& pset) :
  fTrace(false),
  fMaxPErr(0.),
  fGoodPErr(0.),
  fMaxIncChisq(0.),
  fMinLHits(0),
  fMaxLDist(0.),
  fMaxPredDist(0.),
  fMaxPropDist(0.),
  fMinSortDist(0.),
  fMaxSortDist(0.),
  fMaxSamePlane(0),
  fGapDist(0.),
  fMaxNoiseHits(0),
  fPlane(-1)
{
  mf::LogInfo("KalmanFilterAlg") << "KalmanFilterAlg instantiated.";

  // Load fcl parameters.

  reconfigure(pset);
}

/// Destructor.
trkf::KalmanFilterAlg::~KalmanFilterAlg()
{}

/// Reconfigure method.
void trkf::KalmanFilterAlg::reconfigure(const fhicl::ParameterSet& pset)
{
  fTrace = pset.get<bool>("Trace");
  fMaxPErr = pset.get<double>("MaxPErr");
  fGoodPErr = pset.get<double>("GoodPErr");
  fMaxIncChisq = pset.get<double>("MaxIncChisq");
  fMinLHits = pset.get<int>("MinLHits");
  fMaxLDist = pset.get<double>("MaxLDist");
  fMaxPredDist = pset.get<double>("MaxPredDist");
  fMaxPropDist = pset.get<double>("MaxPropDist");
  fMinSortDist = pset.get<double>("MinSortDist");
  fMaxSortDist = pset.get<double>("MaxSortDist");
  fMaxSamePlane = pset.get<int>("MaxSamePlane");
  fGapDist = pset.get<double>("GapDist");
  fMaxNoiseHits = pset.get<double>("MaxNoiseHits");
}

/// Add hits to track.
///
/// Arguments:
///
/// trk      - Starting track.
/// trg      - Result global track.
/// prop     - Propagator.
/// dir      - Direction.
/// hits     - Candidate hits.
///
/// Returns: True if success.
///
/// This method makes a unidirectional Kalman fit in the specified
/// direction, visiting each surface of the passed KHitContainer and
/// updating the track.  In case of multiple measurements on the same
/// surface, keep (at most) the one with the smallest incremental
/// chisquare.  Any measurements that fail the incremental chisquare
/// cut are rejected.  Resolved hits (accepted or rejected) are moved
/// to the unused list in KHitContainer.  The container is sorted at
/// the start of the method, and may be resorted during the progress
/// of the fit.
///
bool trkf::KalmanFilterAlg::buildTrack(const KTrack& trk,
				       KGTrack& trg,
				       const Propagator* prop,
				       const Propagator::PropDirection dir,
				       KHitContainer& hits) const
{
  assert(prop != 0);

  // Direction must be forward or backward (unknown is not allowed).

  if(dir != Propagator::FORWARD && dir != Propagator::BACKWARD)
    throw cet::exception("KalmanFilterAlg") 
	<< "No direction for Kalman fit.\n";

  // Sort container using this seed track.

  hits.sort(trk, true, prop, dir);

  // Loop over measurements (KHitGroup) from sorted list.

  double tchisq = 0.;        // Cumulative chisquare.
  double path = 0.;          // Cumulative path distance.
  int step = 0;              // Step count.
  int nsame = 0;             // Number of consecutive measurements in same plane.
  int last_plane = -1;       // Last plane.

  // Make a copy of the starting track, in the form of a KFitTrack,
  // which we will update as we go.

  TrackError err;
  trk.getSurface()->getStartingError(err);
  KETrack tre(trk, err);
  KFitTrack trf(tre, path, tchisq);

  // Also use the starting track as the reference track for linearized
  // propagation, until the track is established with reasonably small
  // errors.

  KTrack ref(trk);
  KTrack* pref = &ref;

  mf::LogInfo log("KalmanFilterAlg");

  // Loop over measurement groups (KHitGroups).

  while(hits.getSorted().size() > 0) {
    ++step;
    if(fTrace) {
      log << "Build Step " << step << "\n";
      log << "KGTrack has " << trg.numHits() << " hits.\n";
      log << trf;
    }

    // Get an iterator for the next KHitGroup.

    std::list<KHitGroup>::iterator it;
    if(dir == Propagator::FORWARD)
      it = hits.getSorted().begin();
    else {
      assert(dir == Propagator::BACKWARD);
      it = hits.getSorted().end();
      --it;
    }
    const KHitGroup& gr = *it;

    if(fTrace) {
      double path_est = gr.getPath();
      log << "Next surface: " << *(gr.getSurface()) << "\n";
      log << "  Estimated distance = " << path_est << "\n";
    }

    // Get the next prediction surface.  If this KHitGroup is on the
    // preferred plane, use that as the prediction surface.
    // Otherwise, use the current track surface as the prediction
    // surface.

    boost::shared_ptr<const Surface> psurf = trf.getSurface();
    assert(gr.getPlane() >= 0);
    if(fPlane < 0 || gr.getPlane() < 0 || fPlane == gr.getPlane())
      psurf = gr.getSurface();

    // Propagate track to the prediction surface.

    boost::optional<double> dist = prop->noise_prop(trf, psurf, Propagator::UNKNOWN, true, pref);
    if(!!dist && std::abs(*dist) > fMaxPropDist)
      dist = boost::optional<double>(false, 0.);
    double ds = 0.;

    if(!!dist) {

      // Propagation succeeded.
      // Update cumulative path distance and track status.

      ds = *dist;
      path += ds;
      trf.setPath(path);
      if(dir == Propagator::FORWARD)
	trf.setStat(KFitTrack::FORWARD_PREDICTED);
      else {
	assert(dir == Propagator::BACKWARD);
	trf.setStat(KFitTrack::BACKWARD_PREDICTED);
      }
      if(fTrace) {
	log << "After propagation\n";
	log << "  Incremental distance = " << ds << "\n";
	log << "  Actual distance = " << path << "\n";
	log << "KGTrack has " << trg.numHits() << " hits.\n";
	log << trf;
      }

      // Loop over measurements in this group.

      const std::vector<boost::shared_ptr<const KHitBase> >& hits = gr.getHits();
      double best_chisq = 0.;
      boost::shared_ptr<const KHitBase> best_hit;
      for(std::vector<boost::shared_ptr<const KHitBase> >::const_iterator ihit = hits.begin();
	  ihit != hits.end(); ++ihit) {
	const KHitBase& hit = **ihit;

	// Update predction using current track hypothesis and get
	// incremental chisquare.

	bool ok = hit.predict(trf, prop);
	if(ok) {
	  double chisq = hit.getChisq();
	  double preddist = hit.getPredDistance();
	  if((pref != 0 || best_chisq < fMaxIncChisq) && abs(preddist) < fMaxPredDist &&
	     (best_hit.get() == 0 || chisq < best_chisq) ) {
	    best_hit = *ihit;
	    best_chisq = chisq;
	  }
	}
      }
      if(fTrace) {
	if(best_hit.get() != 0) {
	  log << "Hit after prediction\n";
	  log << *best_hit;
	}
      }

      // If we found a best measurement, and if the incremental
      // chisquare passes the cut, add it to the track and update 
      // fit information.

      if(best_hit.get() != 0) {
	ds += best_hit->getPredDistance();
	best_hit->update(trf);
	tchisq += best_chisq;
	trf.setChisq(tchisq);
	if(dir == Propagator::FORWARD)
	  trf.setStat(KFitTrack::FORWARD);
	else {
	  assert(dir == Propagator::BACKWARD);
	  trf.setStat(KFitTrack::BACKWARD);
	}

	// If the pointing error got too large, and there is no
	// reference track, quit.

	if(pref == 0 && trf.PointingError() > fMaxPErr) {
	  if(fTrace)
	    log << "Quitting because pointing error got too large.\n";
	  break;
	}
	  
	// Test number of consecutive measurements in the same plane.

	if(gr.getPlane() >= 0) {
	  if(gr.getPlane() == last_plane)
	    ++nsame;
	  else {
	    nsame = 1;
	    last_plane = gr.getPlane();
	  }
	}
	else {
	  nsame = 0;
	  last_plane = -1;
	}
	if(nsame <= fMaxSamePlane) {

	  // Make a KHitTrack and add it to the KGTrack.

	  KHitTrack trh(trf, best_hit);
	  trg.addTrack(trh);

	  // Decide if we want to kill the reference track.

	  if(pref != 0 && int(trg.numHits()) >= fMinLHits &&
	     (trf.PointingError() < fGoodPErr || path > fMaxLDist)) {
	    pref = 0;
	    if(fTrace)
	      log << "Killing reference track.\n";
	  }

	  if(fTrace) {
	    log << "After update\n";
	    log << "KGTrack has " << trg.numHits() << " hits.\n";
	    log << trf;
	  }
	}
      }
    }

    // The current KHitGroup is now resolved.
    // Move it to unused list.

    hits.getUnused().splice(hits.getUnused().end(), hits.getSorted(), it);

    // If the propagation distance was the wrong direction, resort the measurements.

    if(pref == 0 && !!dist && 
       ((dir == Propagator::FORWARD && (ds < fMinSortDist || ds > fMaxSortDist)) ||
	(dir == Propagator::BACKWARD && (-ds < fMinSortDist || -ds > fMaxSortDist)))) {
      if(fTrace)
	log << "Resorting measurements.\n";
      hits.sort(trf, true, prop, dir);
    }
  }

  // Clean track.

  cleanTrack(trg);

  // Set the fit status of the last added KHitTrack to optimal and get
  // the final chisquare and path length.

  double fchisq = 0.;
  path = 0.;
  if(trg.isValid()) {
    path = trg.endTrack().getPath() - trg.startTrack().getPath();
    if(dir == Propagator::FORWARD) {
      trg.endTrack().setStat(KFitTrack::OPTIMAL);
      fchisq = trg.endTrack().getChisq();
    }
    else {
      assert(dir == Propagator::BACKWARD);
      trg.startTrack().setStat(KFitTrack::OPTIMAL);
      fchisq = trg.startTrack().getChisq();
    }
  }

  // Summary.

  log << "KalmanFilterAlg build track summary.\n"
      << "Build direction = " << (dir == Propagator::FORWARD ? "FORWARD" : "BACKWARD") << "\n"
      << "Track has " << trg.numHits() << " hits.\n"
      << "Track length = " << path << "\n"
      << "Track chisquare = " << fchisq << "\n";

  // Done.  Return success if we added at least one measurement.

  return trg.numHits() > 0;
}

/// Smooth track.
///
/// Arguments:
///
/// trg  - Track to be smoothed.
/// trg1 - Track to receive result of unidirectional fit.
/// prop - Propagator.
///
/// Returns: True if success.
///
/// The starting track should be a global track that has been fit in
/// one direction.  Fit status should be optimal at (at least) one
/// end.  It is an error if the fit status is not optimal at either
/// end.  If the fit status is optimal at both ends, do nothing, but
/// return success.
///
/// If the second argument is non-null, save the result of the
/// unidirectional track fit produced as a byproduct of the smoothing
/// operation.  This track can be smoothed in order to iterate the
/// Kalman fit, etc.
///
/// The Kalman smoothing algorithm starts at the optimal end and fits
/// the track in the reverse direction, calculating optimal track
/// parameters at each measurement surface.
///
/// All measurements are included in the reverse fit.  No incremental
/// chisquare cut is applied.
///
/// If any measurement surface can not be reached because of a
/// measurement error, the entire smoothing operation is considered a
/// failure.  In that case, false is returned and the track is left in
/// an undefined state.
///
bool trkf::KalmanFilterAlg::smoothTrack(KGTrack& trg,
					KGTrack* trg1,
					const Propagator* prop) const
{
  assert(prop != 0);

  // Default result failure.

  bool result = false;

  // It is an error if the KGTrack is not valid.

  if(trg.isValid()) {

    // Examine the track endpoints and figure out which end of the track
    // to start from.  The fit always starts at the optimal end.  It is
    // an error if neither end point is optimal.  Do nothing and return
    // success if both end points are optimal.

    const KHitTrack& trh0 = trg.startTrack();
    const KHitTrack& trh1 = trg.endTrack();
    KFitTrack::FitStatus stat0 = trh0.getStat();
    KFitTrack::FitStatus stat1 = trh1.getStat();
    bool dofit = false;

    // Remember starting direction, track, and distance.

    Propagator::PropDirection dir = Propagator::UNKNOWN;
    const KTrack* trk = 0;
    double path = 0.;

    if(stat0 == KFitTrack::OPTIMAL) {
      if(stat1 == KFitTrack::OPTIMAL) {

	// Both ends optimal (do nothing, return success).

	dofit = false;
	result = true;

      }
      else {

	// Start optimal.

	dofit = true;
	dir = Propagator::FORWARD;
	trk = &trh0;
	path = 0.;
      }
    }
    else {
      if(stat1 == KFitTrack::OPTIMAL) {

	// End optimal.

	dofit = true;
	dir = Propagator::BACKWARD;
	trk = &trh1;
	path = trh1.getPath();
      }
      else {

	// Neither end optimal (do nothing, return failure).

	dofit = false;
	result = false;
      }
    }
    if(dofit) {
      assert(dir == Propagator::FORWARD || dir == Propagator::BACKWARD);
      assert(trk != 0);

      // Cumulative chisquare.

      double tchisq = 0.;

      // Construct starting KFitTrack with track information and distance
      // taken from the optimal end, but with "infinite" errors.

      TrackError err;
      trk->getSurface()->getStartingError(err);
      KETrack tre(*trk, err);
      KFitTrack trf(tre, path, tchisq);

      // Make initial reference track to be same as initial fit track.

      KTrack ref(trf);

      // Loop over KHitTracks contained in KGTrack.

      std::multimap<double, KHitTrack>::iterator it;
      std::multimap<double, KHitTrack>::iterator itend;
      if(dir == Propagator::FORWARD) {
	it = trg.getTrackMap().begin();
	itend = trg.getTrackMap().end();
      }
      else {
	assert(dir == Propagator::BACKWARD);
	it = trg.getTrackMap().end();
	itend = trg.getTrackMap().begin();
      }

      mf::LogInfo log("KalmanFilterAlg");

      // Loop starts here.

      result = true;             // Result success unless we find an error.
      int step = 0;              // Step count.
      while(dofit && it != itend) {
	++step;
	if(fTrace) {
	  log << "Smooth Step " << step << "\n";
	  log << "Reverse fit track:\n";
	  log << trf;
	}

	// For backward fit, decrement iterator at start of loop.

	if(dir == Propagator::BACKWARD)
	  --it;

	KHitTrack& trh = (*it).second;
	if(fTrace) {
	  log << "Forward track:\n";
	  log << trh;
	}

	// Extract measurement.

	const KHitBase& hit = *(trh.getHit());

	// Propagate KFitTrack to the next track surface.
	  
	boost::shared_ptr<const Surface> psurf = trh.getSurface();
	boost::optional<double> dist = prop->noise_prop(trf, psurf, Propagator::UNKNOWN,
							true, &ref);

	// Check if propagation succeeded.  If propagation fails, this
	// measurement will be dropped from the unidirectional fit
	// track.  This measurement will still be in the original
	// track, but with a status other than optimal.

	if(!!dist) {

	  // Propagation succeeded.
	  // Update cumulative path distance and track status.

	  double ds = *dist;
	  path += ds;
	  trf.setPath(path);
	  if(dir == Propagator::FORWARD)
	    trf.setStat(KFitTrack::FORWARD_PREDICTED);
	  else {
	    assert(dir == Propagator::BACKWARD);
	    trf.setStat(KFitTrack::BACKWARD_PREDICTED);
	  }
	  if(fTrace) {
	    log << "Reverse fit track after propagation:\n";
	    log << "  Propagation distance = " << ds << "\n";
	    log << trf;
	  }

	  // See if we have the proper information to calculate an optimal track
	  // at this surface (should normally be possible).

	  KFitTrack::FitStatus stat = trh.getStat();
	  KFitTrack::FitStatus newstat = trf.getStat();

	  if((newstat == KFitTrack::FORWARD_PREDICTED && stat == KFitTrack::BACKWARD) ||
	     (newstat == KFitTrack::BACKWARD_PREDICTED && stat == KFitTrack::FORWARD)) {

	    // Update stored KHitTrack to be optimal.

	    bool ok = trh.combineFit(trf);

	    // Update the stored path distance to be from the currently fitting track.

	    trh.setPath(trf.getPath());

	    // Update reference track.

	    ref = trh;

	    // If combination failed, abandon the fit and return failure.

	    if(!ok) {
	      dofit = false;
	      result = false;
	      break;
	    }
	    if(fTrace) {
	      log << "Combined track:\n";
	      log << trh;
	    }
	  }

	  // Update measurement predction using current track hypothesis.

	  bool ok = hit.predict(trf, prop, &ref);
	  if(!ok) {

	    // If prediction failed, abandon the fit and return failure.

	    dofit = false;
	    result = false;
	    break;
	  }
	  else {

	    // Prediction succeeded.
	    // Get incremental chisquare.
	    // Don't make cut, but do update cumulative chisquare.

	    double chisq = hit.getChisq();
	    tchisq += chisq;
	    trf.setChisq(tchisq);

	    // Update the reverse fitting track using the current measurement
	    // (both track parameters and status).

	    hit.update(trf);
	    if(dir == Propagator::FORWARD)
	      trf.setStat(KFitTrack::FORWARD);
	    else {
	      assert(dir == Propagator::BACKWARD);
	      trf.setStat(KFitTrack::BACKWARD);
	    }
	    if(fTrace) {
	      log << "Reverse fit track after update:\n";
	      log << trf;
	    }

	    // If unidirectional track pointer is not null, make a
	    // KHitTrack and save it in the unidirectional track.

	    if(trg1 != 0) {
	      KHitTrack trh1(trf, trh.getHit());
	      trg1->addTrack(trh1);
	    }
	  }
	}

	// For forward fit, increment iterator at end of loop.

	if(dir == Propagator::FORWARD)
	  ++it;

      }    // Loop over KHitTracks.

      // If fit was successful and the unidirectional track pointer
      // is not null, set the fit status of the last added KHitTrack
      // to optimal.

      if(result && trg1 != 0) {
	if(dir == Propagator::FORWARD)
	  trg1->endTrack().setStat(KFitTrack::OPTIMAL);
	else {
	  assert(dir == Propagator::BACKWARD);
	  trg1->startTrack().setStat(KFitTrack::OPTIMAL);
	}
      }

      // Recalibrate track map.

      trg.recalibrate();

    }      // Do fit.

    // Get the final chisquare.

    double fchisq = 0.5 * (trg.startTrack().getChisq() + trg.endTrack().getChisq());

    // Summary.

    mf::LogInfo log("KalmanFilterAlg");
    log << "KalmanFilterAlg smooth track summary.\n"
	<< "Smooth direction = " << (dir == Propagator::FORWARD ? "FORWARD" : "BACKWARD") << "\n"
	<< "Track has " << trg.numHits() << " hits.\n"
	<< "Track length = " << trg.endTrack().getPath() - trg.startTrack().getPath() << "\n"
	<< "Track chisquare = " << fchisq << "\n";
  }

  // Done.

  return result;
}

/// Add hits to existing track.
///
/// Arguments:
///
/// trg - Track to extend.
/// prop - Propagator.
/// hits - Hit collection to choose hits from.
///
/// This method extends a KGTrack by adding hits.  The KGTrack must
/// have previously been produced by a unidirectional Kalman fit (it
/// should be optimal at one end).  This method finds the optimal end
/// and extends the track in that direction.  If any hits are added,
/// the originally optimal end has its status reset to forward or
/// backward, and the new endpoint is optimal.  In any case, the final
/// result is unidirectionally fit KGTrack.
///
bool trkf::KalmanFilterAlg::extendTrack(KGTrack& trg,
					const Propagator* prop,
					KHitContainer& hits) const
{
  assert(prop != 0);

  // Default result failure.

  bool result = false;

  // Remember the original number of measurement.

  unsigned int nhits0 = trg.numHits();

  // It is an error if the KGTrack is not valid.

  if(trg.isValid()) {
    mf::LogInfo log("KalmanFilterAlg");

    // Examine the track endpoints and figure out which end of the
    // track to extend.  The track is always extended from the optimal
    // end.  It is an error if neither end point is optimal, or both
    // endoints are optimal.  Reset the status of the optimal, and
    // make a copy of the starting track fit.  Also get starting path
    // and chisquare.

    KHitTrack& trh0 = trg.startTrack();
    KHitTrack& trh1 = trg.endTrack();
    KFitTrack::FitStatus stat0 = trh0.getStat();
    KFitTrack::FitStatus stat1 = trh1.getStat();
    bool dofit = false;
    Propagator::PropDirection dir = Propagator::UNKNOWN;
    KFitTrack trf;
    double path = 0.;
    double tchisq = 0.;

    if(stat0 == KFitTrack::OPTIMAL) {
      if(stat1 == KFitTrack::OPTIMAL) {

	// Both ends optimal (do nothing, return failure).

	dofit = false;
	result = false;

      }
      else {

	// Start optimal.  Extend backward.

	dofit = true;
	dir = Propagator::BACKWARD;
	trh0.setStat(KFitTrack::BACKWARD);
	trf = trh0;
	path = trh0.getPath();
	tchisq = trh0.getChisq();
      }
    }
    else {
      if(stat1 == KFitTrack::OPTIMAL) {

	// End optimal.  Extend forward.

	dofit = true;
	dir = Propagator::FORWARD;
	trh1.setStat(KFitTrack::FORWARD);
	trf = trh1;
	path = trh1.getPath();
	tchisq = trh1.getChisq();
      }
      else {

	// Neither end optimal (do nothing, return failure).

	dofit = false;
	result = false;
      }
    }
    if(dofit) {

      // Sort hit container using starting track.

      hits.sort(trf, true, prop, dir);

      // Extend loop starts here.

      int step = 0;
      int nsame = 0;
      int last_plane = -1;
      while(hits.getSorted().size() > 0) {
	++step;
	if(fTrace) {
	  log << "Extend Step " << step << "\n";
	  log << "KGTrack has " << trg.numHits() << " hits.\n";
	  log << trf;
	}

	// Get an iterator for the next KHitGroup.

	std::list<KHitGroup>::iterator it;
	if(dir == Propagator::FORWARD)
	  it = hits.getSorted().begin();
	else {
	  assert(dir == Propagator::BACKWARD);
	  it = hits.getSorted().end();
	  --it;
	}
	const KHitGroup& gr = *it;

	if(fTrace) {
	  double path_est = gr.getPath();
	  log << "Next surface: " << *(gr.getSurface()) << "\n";
	  log << "  Estimated distance = " << path_est << "\n";
	}

	// Get the next prediction surface.  If this KHitGroup is on the
	// preferred plane, use that as the prediction surface.
	// Otherwise, use the current track surface as the prediction
	// surface.

	boost::shared_ptr<const Surface> psurf = trf.getSurface();
	assert(gr.getPlane() >= 0);
	if(fPlane < 0 || gr.getPlane() < 0 || fPlane == gr.getPlane())
	  psurf = gr.getSurface();

	// Propagate track to the prediction surface.

	boost::optional<double> dist = prop->noise_prop(trf, psurf, Propagator::UNKNOWN, true);
	if(!!dist && std::abs(*dist) > fMaxPropDist)
	  dist = boost::optional<double>(false, 0.);
	double ds = 0.;

	if(!!dist) {

	  // Propagation succeeded.
	  // Update cumulative path distance and track status.

	  ds = *dist;
	  path += ds;
	  trf.setPath(path);
	  if(dir == Propagator::FORWARD)
	    trf.setStat(KFitTrack::FORWARD_PREDICTED);
	  else {
	    assert(dir == Propagator::BACKWARD);
	    trf.setStat(KFitTrack::BACKWARD_PREDICTED);
	  }
	  if(fTrace) {
	    log << "After propagation\n";
	    log << "  Incremental distance = " << ds << "\n";
	    log << "  Actual distance = " << path << "\n";
	    log << "KGTrack has " << trg.numHits() << " hits.\n";
	    log << trf;
	  }

	  // Loop over measurements in this group.

	  const std::vector<boost::shared_ptr<const KHitBase> >& hits = gr.getHits();
	  double best_chisq = 0.;
	  boost::shared_ptr<const KHitBase> best_hit;
	  for(std::vector<boost::shared_ptr<const KHitBase> >::const_iterator ihit = hits.begin();
	      ihit != hits.end(); ++ihit) {
	    const KHitBase& hit = **ihit;

	    // Update predction using current track hypothesis and get
	    // incremental chisquare.

	    bool ok = hit.predict(trf, prop);
	    if(ok) {
	      double chisq = hit.getChisq();
	      double preddist = hit.getPredDistance();
	      if(chisq < fMaxIncChisq && abs(preddist) < fMaxPredDist &&
		 (best_hit.get() == 0 || chisq < best_chisq)) {
		best_hit = *ihit;
		best_chisq = chisq;
	      }
	    }
	  }
	  if(fTrace) {
	    if(best_hit.get() != 0) {
	      log << "Hit after prediction\n";
	      log << *best_hit;
	    }
	  }

	  // If we found a best measurement, and if the incremental
	  // chisquare passes the cut, add it to the track and update 
	  // fit information.

	  if(best_hit.get() != 0) {
	    ds += best_hit->getPredDistance();
	    best_hit->update(trf);
	    tchisq += best_chisq;
	    trf.setChisq(tchisq);
	    if(dir == Propagator::FORWARD)
	      trf.setStat(KFitTrack::FORWARD);
	    else {
	      assert(dir == Propagator::BACKWARD);
	      trf.setStat(KFitTrack::BACKWARD);
	    }

	    // If the pointing error got too large, quit.

	    if(trf.PointingError() > fMaxPErr) {
	      if(fTrace)
		log << "Quitting because pointing error got too large.\n";
	      break;
	    }

	    // Test number of consecutive measurements in the same plane.

	    if(gr.getPlane() >= 0) {
	      if(gr.getPlane() == last_plane)
		++nsame;
	      else {
		nsame = 1;
		last_plane = gr.getPlane();
	      }
	    }
	    else {
	      nsame = 0;
	      last_plane = -1;
	    }
	    if(nsame <= fMaxSamePlane) {

	      // Make a KHitTrack and add it to the KGTrack.

	      KHitTrack trh(trf, best_hit);
	      trg.addTrack(trh);

	      if(fTrace) {
		log << "After update\n";
		log << "KGTrack has " << trg.numHits() << " hits.\n";
		log << trf;
	      }
	    }
	  }
	}

	// The current KHitGroup is now resolved.
	// Move it to unused list.

	hits.getUnused().splice(hits.getUnused().end(), hits.getSorted(), it);

	// If the propagation distance was the wrong direction, resort the measurements.

	if(!!dist &&
	   ((dir == Propagator::FORWARD && (ds < fMinSortDist || ds > fMaxSortDist)) ||
	    (dir == Propagator::BACKWARD && (-ds < fMinSortDist || -ds > fMaxSortDist)))) {
	  if(fTrace)
	    log << "Resorting measurements.\n";
	  hits.sort(trf, true, prop, dir);
	}
      }
    }

    // Clean track.

    cleanTrack(trg);

    // Set the fit status of the last added KHitTrack to optimal and
    // get the final chisquare and path length.

    double fchisq = 0.;
    path = 0.;
    if(trg.isValid()) {
      path = trg.endTrack().getPath() - trg.startTrack().getPath();
      if(dir == Propagator::FORWARD) {
	trg.endTrack().setStat(KFitTrack::OPTIMAL);
	fchisq = trg.endTrack().getChisq();
      }
      else {
	assert(dir == Propagator::BACKWARD);
	trg.startTrack().setStat(KFitTrack::OPTIMAL);
	fchisq = trg.startTrack().getChisq();
      }
    }

    // Summary.

    log << "KalmanFilterAlg extend track summary.\n"
	<< "Extend direction = " << (dir == Propagator::FORWARD ? "FORWARD" : "BACKWARD") << "\n"
	<< "Track has " << trg.numHits() << " hits.\n"
	<< "Track length = " << path << "\n"
	<< "Track chisquare = " << fchisq << "\n";
  }

  // Done.

  result = (trg.numHits() > nhits0);
  return result;
}

/// Iteratively smooth a track.
///
/// Arguments:
///
/// nsmooth - Number of iterations.
/// trg     - Track to be smoothed.
/// prop    - Propagator.
///
/// Returns: True if success.
///
/// The initial track should have been unidirectionally fit.
///
bool trkf::KalmanFilterAlg::smoothTrackIter(int nsmooth,
					    KGTrack& trg,
					    const Propagator* prop) const
{
  bool ok = true;

  // Call smoothTrack method in a loop.

  for(int ismooth = 0; ok && ismooth < nsmooth-1; ++ismooth) {
    KGTrack trg1;
    ok = smoothTrack(trg, &trg1, prop);
    if(ok)
      trg = trg1;
  }

  // Do one final smooth without generating a new unidirectional
  // track.

  if(ok)
    ok = smoothTrack(trg, 0, prop);

  // Done.

  return ok;
}

/// Clean track by removing noise hits near endpoints.
///
/// Arguments:
///
/// trg - Track.
///
void trkf::KalmanFilterAlg::cleanTrack(KGTrack& trg) const
{
  // Get hold of a modifiable track map from trg.

  std::multimap<double, KHitTrack>& trackmap = trg.getTrackMap();

  // Do an indefinite loop until we no longer find any dirt.

  bool done = false;
  while(!done) {

    // If track map has fewer than fMaxNoiseHits tracks, then this is a
    // noise track.  Clear the map, making the whole track invalid.

    if(int(trackmap.size()) <= fMaxNoiseHits) {
      trackmap.clear();
      done = true;
      break;
    }

    // Make sure the first two and last two tracks belong to different
    // views.  If not, remove the first or last track.

    if(trackmap.size() >= 2) {

      // Check start.

      std::multimap<double, KHitTrack>::iterator it = trackmap.begin();
      int plane1 = (*it).second.getHit()->getMeasPlane();
      if(plane1 >= 0) {
	++it;
	int plane2 = (*it).second.getHit()->getMeasPlane();
	if(plane2 >= 0 && plane1 == plane2) {
	  trackmap.erase(trackmap.begin(), it);
	  done = false;
	  continue;
	}
      }

      // Check end.

      it = trackmap.end();
      --it;
      plane1 = (*it).second.getHit()->getMeasPlane();
      if(plane1 >= 0) {
	--it;
	int plane2 = (*it).second.getHit()->getMeasPlane();
	if(plane2 >= 0 && plane1 == plane2) {
	  ++it;
	  trackmap.erase(it, trackmap.end());
	  done = false;
	  continue;
	}
      }
    }

    // Loop over successive pairs of elements of track map.  Look for
    // adjacent elements with distance separation greater than
    // fGapDist.

    std::multimap<double, KHitTrack>::iterator it = trackmap.begin();
    std::multimap<double, KHitTrack>::iterator jt = trackmap.end();
    int nb = 0;                // Number of elements from begin to jt.
    int ne = trackmap.size();  // Number of elements it to end.
    bool found_noise = false;
    for(; it != trackmap.end(); ++it, ++nb, --ne) {
      if(jt == trackmap.end())
	jt = trackmap.begin();
      else {
	assert(nb >= 1);
	assert(ne >= 1);
	double disti = (*it).first;
	double distj = (*jt).first;
	double sep = disti - distj;
	assert(sep >= 0.);
	if(sep > fGapDist) {

	  // Found a gap.  See if we want to trim track.

	  if(nb <= fMaxNoiseHits) {

	    // Trim front.

	    found_noise = true;
	    trackmap.erase(trackmap.begin(), it);
	    break;
	  }
	  if(ne <= fMaxNoiseHits) {

	    // Trim back.

	    found_noise = true;
	    trackmap.erase(it, trackmap.end());
	    break;
	  }
	}
	++jt;
      }
    }

    // Decide if we are done.

    done = !found_noise;
  }
}
