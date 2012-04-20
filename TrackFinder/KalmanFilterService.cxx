///////////////////////////////////////////////////////////////////////
///
/// \file   KalmanFilterService.cxx
///
/// \brief  Kalman Filter.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include "TrackFinder/KalmanFilterService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

/// Constructor.
  
trkf::KalmanFilterService::KalmanFilterService(const fhicl::ParameterSet& pset,
					   art::ActivityRegistry& reg)
{
  mf::LogInfo("KalmanFilterService") << "KalmanFilterService instantiated.";
}

/// Destructor.
trkf::KalmanFilterService::~KalmanFilterService()
{}

/// Add hits to track.
///
/// Arguments:
///
/// trh      - Starting track.
/// prop     - Propagator.
/// dir      - Direction.
/// hits     - Candidate hits.
/// maxChisq - Incremental chisquare cut.
///
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
void trkf::KalmanFilterService::buildTrack(KHitsTrack& trh,
					   const Propagator* prop,
					   const Propagator::PropDirection dir,
					   KHitContainer& hits,
					   double maxChisq)
{
  // Remember starting path distance.

  assert(prop != 0);

  // Direction must be forward or backward (unknown is not allowed).

  if(dir != Propagator::FORWARD && dir != Propagator::BACKWARD)
    throw cet::exception("KalmanFilterService") 
	<< "No direction for Kalman fit.\n";

  // Sort container using this seed track.

  hits.sort(trh, true, prop, dir);

  // Loop over measurements (KHitGroup) from sorted list.

  double tchisq = trh.getChisq();   // Cumulative chisquare.
  double path = trh.getPath();     // Cumulative path distance.
  int step = 0;

  mf::LogInfo log("KalmanFilterService");

  while(hits.getSorted().size() > 0) {
    ++step;
    double xyz[3];
    double mom[3];
    trh.getPosition(xyz);
    trh.getMomentum(mom);
    log << "Step " << step << "\n";
    log << "  dist= " << path << "\n";
    log << "  chisq= " << tchisq << "\n";
    log << "   x= " << xyz[0] << ",  y= " << xyz[1] << ",  z= " << xyz[2] << "\n";
    log << "   px= " << mom[0] << ", py= " << mom[1] << ", pz= " << mom[2] << "\n";
    log << trh;

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

    double path_est = gr.getPath();
    log << "Next surface: " << *(gr.getSurface()) << "\n";
    log << "  Estimated distance = " << path_est << "\n";

    // Propagate track to the surface of the KHitGroup.

    boost::shared_ptr<const Surface> psurf = gr.getSurface();
    boost::optional<double> dist = prop->noise_prop(trh, psurf, Propagator::UNKNOWN, true);
    if(!!dist) {

      // Propagation succeeded.
      // Update cumulative path distance and track status.

      double ds = *dist;
      path += ds;
      trh.setPath(path);
      if(dir == Propagator::FORWARD)
	trh.setStat(KFitTrack::FORWARD_PREDICTED);
      else {
	assert(dir == Propagator::BACKWARD);
	trh.setStat(KFitTrack::BACKWARD_PREDICTED);
      }
      trh.getPosition(xyz);
      trh.getMomentum(mom);
      log << "After propagation\n";
      log << "   x= " << xyz[0] << ",  y= " << xyz[1] << ",  z= " << xyz[2] << "\n";
      log << "   px= " << mom[0] << ", py= " << mom[1] << ", pz= " << mom[2] << "\n";
      log << "  Actual distance = " << path << "\n";
      log << trh;

      // Loop over measurements in this group.

      const std::vector<boost::shared_ptr<const KHitBase> >& hits = gr.getHits();
      double best_chisq = 0.;
      boost::shared_ptr<const KHitBase> best_hit;
      for(std::vector<boost::shared_ptr<const KHitBase> >::const_iterator ihit = hits.begin();
	  ihit != hits.end(); ++ihit) {
	const KHitBase& hit = **ihit;

	// Update predction using current track hypothesis and get
	// incremental chisquare.

	bool ok = hit.predict(trh);
	if(ok) {
	  double chisq = hit.getChisq();
	  if(best_hit.get() == 0 || chisq < best_chisq) {
	    best_hit = *ihit;
	    best_chisq = chisq;
	  }
	}
      }
      if(best_hit.get() != 0) {
	log << "Hit after prediction\n";
	log << *best_hit;
      }

      // If we found a best measurement, and if the incremental
      // chisquare passes the cut, add it to the track and update 
      // fit information.

      if(best_hit.get() != 0 && best_chisq < maxChisq) {
	best_hit->update(trh);
	trh.addHit(best_hit);
	tchisq += best_chisq;
	trh.setChisq(tchisq);
	if(dir == Propagator::FORWARD)
	  trh.setStat(KFitTrack::FORWARD);
	else {
	  assert(dir == Propagator::BACKWARD);
	  trh.setStat(KFitTrack::BACKWARD);
	}
	trh.getPosition(xyz);
	trh.getMomentum(mom);
	log << "After update\n";
	log << "   x= " << xyz[0] << ",  y= " << xyz[1] << ",  z= " << xyz[2] << "\n";
	log << "   px= " << mom[0] << ", py= " << mom[1] << ", pz= " << mom[2] << "\n";
	log << trh;
      }
    }

    // The current KHitGroup is now resolved.
    // Move it to unused list.

    hits.getUnused().splice(hits.getUnused().end(), hits.getSorted(), it);
  }
}
