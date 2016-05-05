////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       ProjectionMatchingAlg
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "larreco/RecoAlg/ProjectionMatchingAlg.h"

#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

pma::ProjectionMatchingAlg::ProjectionMatchingAlg(const fhicl::ParameterSet& pset)
  : fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>())
{
	this->reconfigure(pset); 
}
// ------------------------------------------------------

pma::ProjectionMatchingAlg::~ProjectionMatchingAlg(void)
{
}
// ------------------------------------------------------

void pma::ProjectionMatchingAlg::reconfigure(const fhicl::ParameterSet& p)
{
	fOptimizationEps = p.get< double >("OptimizationEps");
	fFineTuningEps = p.get< double >("FineTuningEps");

	fTrkValidationDist2D = p.get< double >("TrkValidationDist2D");
	fHitTestingDist2D = p.get< double >("HitTestingDist2D");

	fMinTwoViewFraction = p.get< double >("MinTwoViewFraction");
	
	fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
	
	pma::Node3D::SetMargin(p.get< double >("NodeMargin3D"));

	pma::Element3D::SetOptFactor(geo::kZ, p.get< double >("HitWeightZ"));
	pma::Element3D::SetOptFactor(geo::kV, p.get< double >("HitWeightV"));
	pma::Element3D::SetOptFactor(geo::kU, p.get< double >("HitWeightU"));
}
// ------------------------------------------------------

double pma::ProjectionMatchingAlg::validate(const pma::Track3D& trk,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	unsigned int testView) const
{
	double max_d = fTrkValidationDist2D;
	double d2, max_d2 = max_d * max_d;
	unsigned int nAll = 0, nPassed = 0;

	std::vector< unsigned int > trkTPCs = trk.TPCs();
	std::vector< unsigned int > trkCryos = trk.Cryos();
	std::map< std::pair< unsigned int, unsigned int >, std::pair< TVector2, TVector2 > > ranges;
	std::map< std::pair< unsigned int, unsigned int >, double > wirePitch;
	for (auto c : trkCryos)
		for (auto t : trkTPCs)
		{
			ranges[std::pair< unsigned int, unsigned int >(t, c)] = trk.WireDriftRange(testView, t, c);
			wirePitch[std::pair< unsigned int, unsigned int >(t, c)] = fGeom->TPC(t, c).Plane(testView).WirePitch();
		}

	unsigned int tpc, cryo;
	std::map< std::pair< unsigned int, unsigned int >, std::vector< TVector2 > > all_close_points;

	for (const auto h : hits)
		if (h->WireID().Plane == testView)
	{
		tpc = h->WireID().TPC;
		cryo = h->WireID().Cryostat;
		std::pair< unsigned int, unsigned int > tpc_cryo(tpc, cryo);
		std::pair< TVector2, TVector2 > rect = ranges[tpc_cryo];

		if ((h->WireID().Wire > rect.first.X() - 10) &&  // chceck only hits in the rectangle around
		    (h->WireID().Wire < rect.second.X() + 10) && // the track projection, it is faster than
		    (h->PeakTime() > rect.first.Y() - 100) &&    // calculation of trk.Dist2(p2d, testView)
		    (h->PeakTime() < rect.second.Y() + 100))
		{
			TVector2 p2d(wirePitch[tpc_cryo] * h->WireID().Wire, fDetProp->ConvertTicksToX(h->PeakTime(), testView, tpc, cryo));

			d2 = trk.Dist2(p2d, testView, tpc, cryo);

			if (d2 < max_d2) all_close_points[tpc_cryo].push_back(p2d);
		}
	}

	auto const & channelStatus = art::ServiceHandle< lariov::ChannelStatusService >()->GetProvider();

	double step = 0.3;
	// then check how points-close-to-the-track-projection are distributed along the
	// track, namely: are there track sections crossing empty spaces, except dead wires?
	TVector3 p(trk.front()->Point3D());
	for (size_t i = 0; i < trk.Nodes().size() - 1; ++i)
	{
		auto const & node = *(trk.Nodes()[i]);
		tpc = node.TPC(); cryo = node.Cryo();

		TVector3 vNext(trk.Nodes()[i + 1]->Point3D());
		TVector3 vThis(node.Point3D());

		std::vector< TVector2 > const & points = all_close_points[std::pair< unsigned int, unsigned int >(tpc, cryo)];
		if (trk.Nodes()[i + 1]->TPC() == (int)tpc) // skip segments between tpc's, look only at those contained in tpc
		{
			TVector3 dc(vNext); dc -= vThis;
			dc *= step / dc.Mag();

			double f = pma::GetSegmentProjVector(p, vThis, vNext);
			double wirepitch = fGeom->TPC(tpc, cryo).Plane(testView).WirePitch();
			while ((f < 1.0) && node.SameTPC(p))
			{
				TVector2 p2d(fGeom->WireCoordinate(p.Y(), p.Z(), testView, tpc, cryo), p.X());
				geo::WireID wireID(cryo, tpc, testView, (int)p2d.X());
				if (fGeom->HasWire(wireID))
				{
				  raw::ChannelID_t ch = fGeom->PlaneWireToChannel(wireID);
				  if (channelStatus.IsGood(ch))
				  {
				        if (points.size())
					{
						p2d.Set(wirepitch * p2d.X(), p2d.Y());
						for (const auto & h : points)
						{
							d2 = pma::Dist2(p2d, h);
							if (d2 < max_d2) { nPassed++; break; }
						}
					}
					nAll++;
				  }
				  //else mf::LogVerbatim("ProjectionMatchingAlg")
				  //	<< "crossing BAD CHANNEL (wire #" << (int)p2d.X() << ")" << std::endl;
				}

				p += dc; f = pma::GetSegmentProjVector(p, vThis, vNext);
			}
		}
		p = vNext;
	}

	if (nAll > 0)
	{
		double v = nPassed / (double)nAll;
		mf::LogVerbatim("ProjectionMatchingAlg") << "  trk fraction ok: " << v;
		return v;
	}
	else return 1.0;
}
// ------------------------------------------------------

double pma::ProjectionMatchingAlg::validate(
	const TVector3& p0, const TVector3& p1,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	unsigned int testView, unsigned int tpc, unsigned int cryo) const
{
	double step = 0.3;
	double max_d = fTrkValidationDist2D;
	double d2, max_d2 = max_d * max_d;
	unsigned int nAll = 0, nPassed = 0;

	TVector3 p(p0);
	TVector3 dc(p1); dc -= p;
	dc *= step / dc.Mag();

	auto const & channelStatus = art::ServiceHandle< lariov::ChannelStatusService >()->GetProvider();

	double f = pma::GetSegmentProjVector(p, p0, p1);
	double wirepitch = fGeom->TPC(tpc, cryo).Plane(testView).WirePitch();
	while (f < 1.0)
	{
		TVector2 p2d(fGeom->WireCoordinate(p.Y(), p.Z(), testView, tpc, cryo), p.X());
		geo::WireID wireID(cryo, tpc, testView, (int)p2d.X());
		if (fGeom->HasWire(wireID))
		{
			raw::ChannelID_t ch = fGeom->PlaneWireToChannel(wireID);
			if (channelStatus.IsGood(ch))
			{
				p2d.Set(wirepitch * p2d.X(), p2d.Y());
				for (const auto & h : hits)
					if (h->WireID().Plane == testView)
				{
					d2 = pma::Dist2(p2d, pma::WireDriftToCm(h->WireID().Wire, h->PeakTime(), testView, tpc, cryo));
					if (d2 < max_d2) { nPassed++; break; }
				}
				nAll++;
			}
			//else mf::LogVerbatim("ProjectionMatchingAlg")
			//	<< "crossing BAD CHANNEL (wire #" << (int)p2d.X() << ")" << std::endl;
		}

		p += dc; f = pma::GetSegmentProjVector(p, p0, p1);
	}

	if (nAll > 3) // validate actually only if 2D projection in testView has some minimum length
	{
		double v = nPassed / (double)nAll;
		mf::LogVerbatim("ProjectionMatchingAlg") << "  segment fraction ok: " << v;
		return v;
	}
	else return 1.0;
}
// ------------------------------------------------------

double pma::ProjectionMatchingAlg::twoViewFraction(pma::Track3D& trk) const
{
	trk.SelectHits();
	trk.DisableSingleViewEnds();

	size_t idx = 0;
	while ((idx < trk.size() - 1) && !trk[idx]->IsEnabled()) idx++;
	double l0 = trk.Length(0, idx + 1);

	idx = trk.size() - 1;
	while ((idx > 1) && !trk[idx]->IsEnabled()) idx--;
	double l1 = trk.Length(idx - 1, trk.size() - 1);

	trk.SelectHits();

	return 1.0 - (l0 + l1) / trk.Length();
}
// ------------------------------------------------------

size_t pma::ProjectionMatchingAlg::getSegCount(size_t trk_size)
{
	int nSegments = (int)( 0.8 * trk_size / sqrt(trk_size) );

	if (nSegments > 1) return (size_t)nSegments;
	else return 1;
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildTrack(
	const std::vector< art::Ptr<recob::Hit> >& hits_1,
	const std::vector< art::Ptr<recob::Hit> >& hits_2) const
{
	pma::Track3D* trk = new pma::Track3D(); // track candidate
	trk->AddHits(hits_1);
	trk->AddHits(hits_2);

	mf::LogVerbatim("ProjectionMatchingAlg") << "track size: " << trk->size();
	std::vector< unsigned int > tpcs = trk->TPCs();
	for (size_t t = 0; t < tpcs.size(); ++t)
	{
		mf::LogVerbatim("ProjectionMatchingAlg") << "  tpc:" << tpcs[t];
	}
	mf::LogVerbatim("ProjectionMatchingAlg")
		<< "  #coll:" << trk->NHits(geo::kZ)
		<< "  #ind2:" << trk->NHits(geo::kV)
		<< "  #ind1:" << trk->NHits(geo::kU);

	size_t nSegments = getSegCount(trk->size());
	size_t nNodes = (size_t)( nSegments - 1 ); // n nodes to add

	mf::LogVerbatim("ProjectionMatchingAlg") << "  initialize trk";
	if (!trk->Initialize())
	{
		mf::LogWarning("ProjectionMatchingAlg") << "  initialization failed, delete trk";
		delete trk;
		return 0;
	}

	double f = twoViewFraction(*trk);
	if (f > fMinTwoViewFraction)
	{
		double g = 0.0;
		mf::LogVerbatim("ProjectionMatchingAlg") << "  optimize trk (" << nSegments << " seg)";
		if (nNodes)
		{
			g = trk->Optimize(nNodes, fOptimizationEps, false, true, 25, 10);   // build nodes
			mf::LogVerbatim("ProjectionMatchingAlg") << "  nodes done, g = " << g;
			trk->SelectAllHits();
		}
		g = trk->Optimize(0, fFineTuningEps);              // final tuning
		mf::LogVerbatim("ProjectionMatchingAlg") << "  tune done, g = " << g;

		trk->SortHits();
		// trk->ShiftEndsToHits(); // not sure if useful already here
		return trk;
	}
	else
	{
		mf::LogVerbatim("ProjectionMatchingAlg") << "  clusters do not match, f = " << f;
		delete trk;
		return 0;
	}
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildMultiTPCTrack(
	const std::vector< art::Ptr<recob::Hit> >& hits) const
{
	std::map< unsigned int, std::vector< art::Ptr<recob::Hit> > > hits_by_tpc;
	for (auto const & h : hits)
	{
		hits_by_tpc[h->WireID().TPC].push_back(h);
	}

	std::vector< pma::Track3D* > tracks;
	for(auto const & hsel : hits_by_tpc)
	{
		pma::Track3D* trk = buildTrack(hsel.second);
		if (trk) tracks.push_back(trk);
	}

	bool need_reopt = false;
	while (tracks.size() > 1)
	{
		need_reopt = true;

		pma::Track3D* first = tracks.front();
		pma::Track3D* best = 0;
		double d, dmin = 1.0e12;
		size_t t_best = 0, cfg = 0;
		for (size_t t = 1; t < tracks.size(); ++t)
		{
			pma::Track3D* second = tracks[t];

			d = pma::Dist2(first->front()->Point3D(), second->front()->Point3D());
			if (d < dmin) { dmin = d; best = second; t_best = t; cfg = 0; }

			d = pma::Dist2(first->front()->Point3D(), second->back()->Point3D());
			if (d < dmin) { dmin = d; best = second; t_best = t; cfg = 1; }

			d = pma::Dist2(first->back()->Point3D(), second->front()->Point3D());
			if (d < dmin) { dmin = d; best = second; t_best = t; cfg = 2; }

			d = pma::Dist2(first->back()->Point3D(), second->back()->Point3D());
			if (d < dmin) { dmin = d; best = second; t_best = t; cfg = 3; }
		}
		if (best)
		{
			switch (cfg)
			{
				default:
				case 0:
				case 1:
					mergeTracks(*best, *first, false);
					tracks[0] = best; delete first; break;

				case 2:
				case 3:
					mergeTracks(*first, *best, false);
					delete best; break;
			}
			tracks.erase(tracks.begin() + t_best);
		}
		else break; // should not happen
	}

	if (!tracks.empty())
	{
		pma::Track3D* trk = tracks.front();
		if (need_reopt)
		{
			double g = trk->Optimize(0, fFineTuningEps);
			mf::LogVerbatim("ProjectionMatchingAlg") << "  reopt after merging tpc parts: done, g = " << g;
		}
		return trk;
	}
	else return 0;
}

// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildShowerSeg(
		const std::vector< art::Ptr<recob::Hit> >& hits, 
		const art::Ptr<recob::Vertex>& vtx) const
{
	double vtxarray[3]; 
	vtx->XYZ(vtxarray);	

	if (!fGeom->HasTPC(fGeom->FindTPCAtPosition(vtxarray))) return 0;

	TVector3 vtxv3(vtxarray[0], vtxarray[1], vtxarray[2]);
	
	const size_t tpc = fGeom->FindTPCAtPosition(vtxarray).TPC;
	const size_t cryo = fGeom->FindCryostatAtPosition(vtxarray);
	const geo::TPCGeo& tpcgeom = fGeom->Cryostat(cryo).TPC(tpc);

	// use only hits from tpc where the vtx is
	std::vector<art::Ptr<recob::Hit> > hitstpc;
	for (size_t h = 0; h < hits.size(); ++h)
		if (hits[h]->WireID().TPC == tpc)
			hitstpc.push_back(hits[h]);

	if (!hitstpc.size()) return 0;
	
	std::vector<art::Ptr<recob::Hit> > hitstrk;
	size_t view = 0; size_t countviews = 0;
	while (view < 3)
	{
		mf::LogWarning("ProjectionMatchinAlg") << "collecting hits from view: " << view;
		if (!tpcgeom.HasPlane(view)) {++view; continue;}

		// select hits only for a single view
		std::vector<art::Ptr<recob::Hit> > hitsview;
		for (size_t h = 0; h < hitstpc.size(); ++h)
			if (hitstpc[h]->WireID().Plane == view)
				hitsview.push_back(hitstpc[h]);
		if (!hitsview.size()) {++view; continue;}

		// filter our small groups of hits, far from main shower		
		std::vector<art::Ptr<recob::Hit> > hitsfilter;
		TVector2 proj_pr = pma::GetProjectionToPlane(vtxv3, view, tpc, cryo);
	
		mf::LogWarning("ProjectionMatchinAlg") << "start filter out: ";
		FilterOutSmallParts(2.0, hitsview, hitsfilter, proj_pr);
		mf::LogWarning("ProjectionMatchingAlg") << "after filter out";

		for (size_t h = 0; h < hitsfilter.size(); ++h)
			hitstrk.push_back(hitsfilter[h]);

		if (hitsfilter.size() >= 5) countviews++;

		++view;
	}

	if (!hitstrk.size() || (countviews < 2)) 
	{
		mf::LogWarning("ProjectionMatchinAlg") << "too few hits, segment not built";
		return 0;
	}

	// hits are prepared, finally segment is built

	pma::Track3D* trk = new pma::Track3D();
	trk = buildSegment(hitstrk, vtxv3);
	
	// make shorter segment to estimate direction more precise
	ShortenSeg(*trk, tpcgeom);

	if (trk->size() < 10)
	{
		mf::LogWarning("ProjectionMatchingAlg") << "too short segment, delete segment";
		delete trk;
		return 0;
	}

	return trk;
}

// ------------------------------------------------------

void pma::ProjectionMatchingAlg::FilterOutSmallParts(
		double r2d,
		const std::vector< art::Ptr<recob::Hit> >& hits_in,
		std::vector< art::Ptr<recob::Hit> >& hits_out,
		const TVector2& vtx2d) const
{
	size_t min_size = hits_in.size() / 5;
	if (min_size < 3) min_size = 3;

	std::vector<size_t> used;
	std::vector< std::vector< art::Ptr<recob::Hit> > > sub_groups;
	std::vector< art::Ptr<recob::Hit> > close_hits;

	float mindist2 = 99999.99; size_t id_sub_groups_save = 0; size_t id = 0;
	while (GetCloseHits(r2d, hits_in, used, close_hits)) 
	{
		
		sub_groups.emplace_back(close_hits);

		for (size_t h = 0; h < close_hits.size(); ++h)
		{
			TVector2 hi_cm = pma::WireDriftToCm(close_hits[h]->WireID().Wire, 
			close_hits[h]->PeakTime(), 
			close_hits[h]->WireID().Plane, 
			close_hits[h]->WireID().TPC, 
			close_hits[h]->WireID().Cryostat);

			float dist2 = pma::Dist2(hi_cm, vtx2d);
			if (dist2 < mindist2)
			{
				id_sub_groups_save = id;
				mindist2 = dist2;
			}
		}
		
		id++;
	}

	for (size_t i = 0; i < sub_groups.size(); ++i)
	{
		if (i == id_sub_groups_save)
		{
			for (auto h : sub_groups[i]) hits_out.push_back(h);
		}
	}

	for (size_t i = 0; i < sub_groups.size(); ++i)
	{
		if ((i != id_sub_groups_save) && (hits_out.size() < 10) && (sub_groups[i].size() > min_size))
			for (auto h : sub_groups[i]) hits_out.push_back(h);
	}
}

// ------------------------------------------------------

bool pma::ProjectionMatchingAlg::GetCloseHits(
		double r2d, 
		const std::vector< art::Ptr<recob::Hit> >& hits_in, 
		std::vector<size_t>& used,
		std::vector< art::Ptr<recob::Hit> >& hits_out) const
{
	hits_out.clear();

	const double gapMargin = 5.0; // can be changed to f(id_tpc1, id_tpc2)
	size_t idx = 0;

	while ((idx < hits_in.size()) && Has(used, idx)) idx++;

	if (idx < hits_in.size())
	{		
		hits_out.push_back(hits_in[idx]);
		used.push_back(idx);

		unsigned int tpc = hits_in[idx]->WireID().TPC;
		unsigned int cryo = hits_in[idx]->WireID().Cryostat;
		unsigned int view = hits_in[idx]->WireID().Plane;
		double wirePitch = fGeom->TPC(tpc, cryo).Plane(view).WirePitch();
		double driftPitch = fDetProp->GetXTicksCoefficient(tpc, cryo);

		double r2d2 = r2d*r2d;
		double gapMargin2 = sqrt(2 * gapMargin*gapMargin);
		gapMargin2 = (gapMargin2 + r2d)*(gapMargin2 + r2d);

		bool collect = true;
		while (collect)
		{
			collect = false;
			for (size_t i = 0; i < hits_in.size(); i++)
				if (!Has(used, i))
				{
					art::Ptr<recob::Hit> const & hi = hits_in[i];
					TVector2 hi_cm(wirePitch * hi->WireID().Wire, driftPitch * hi->PeakTime());

					bool accept = false;
		
					for (auto const & ho : hits_out)					
					{

						TVector2 ho_cm(wirePitch * ho->WireID().Wire, driftPitch * ho->PeakTime());
						double d2 = pma::Dist2(hi_cm, ho_cm);
						
						if (d2 < r2d2) { accept = true; break; }
					}
					if (accept)
					{
						collect = true;
						hits_out.push_back(hi);
						used.push_back(i);
					}
				}
		}
		return true;
	}
	else return false;
}

// ------------------------------------------------------

void pma::ProjectionMatchingAlg::ShortenSeg(pma::Track3D& trk, const geo::TPCGeo& tpcgeom) const
{
	double mse = trk.GetMse();
	mf::LogWarning("ProjectionMatchingAlg") << "initial value of mse: " << mse;
	while ((mse > 0.5) && TestTrk(trk, tpcgeom))
	{
		mse = trk.GetMse();
		for (size_t h = 0; h < trk.size(); ++h)
			if (std::sqrt(pma::Dist2(trk[h]->Point3D(), trk.front()->Point3D())) 
					> 0.8*trk.Length()) trk[h]->SetEnabled(false);

		RemoveNotEnabledHits(trk);

		trk.Optimize(0.0001, false);
		trk.SortHits();

		mf::LogWarning("ProjectionMatchingAlg") << " mse: " << mse;
		if (mse == trk.GetMse()) break;
		else mse = trk.GetMse();
	}

	RemoveNotEnabledHits(trk);
}

// ------------------------------------------------------

bool pma::ProjectionMatchingAlg::TestTrk(pma::Track3D& trk, const geo::TPCGeo& tpcgeom) const
{
	bool test = false; 

	if (tpcgeom.HasPlane(0) && tpcgeom.HasPlane(1))
	{
		if ((trk.NEnabledHits(0) > 5) && (trk.NEnabledHits(1) > 5)) 
			test = true; 
	}
	else if (tpcgeom.HasPlane(0) && tpcgeom.HasPlane(2))
	{
		if ((trk.NEnabledHits(0) > 5) && (trk.NEnabledHits(2) > 5)) 
			test = true; 
	}
	else if (tpcgeom.HasPlane(1) && tpcgeom.HasPlane(2))
	{
		if ((trk.NEnabledHits(1) > 5) && (trk.NEnabledHits(2) > 5)) 
			test = true;
	}

	double length = 0.0;
	for (size_t h = 0; h < trk.size(); ++h)
	{
		if (!trk[h]->IsEnabled()) break;
		length = std::sqrt(pma::Dist2(trk.front()->Point3D(), trk[h]->Point3D()));
	}

	mf::LogWarning("ProjectionMatchingAlg") << "length of segment: " << length;
	if (length < 3.0) test = false; // cm
	
	return test;
}

// ------------------------------------------------------

bool pma::ProjectionMatchingAlg::Has(const std::vector<size_t>& v, size_t idx) const
{
    for (auto c : v) if (c == idx) return true;
    return false;
}

// ------------------------------------------------------

void pma::ProjectionMatchingAlg::RemoveNotEnabledHits(pma::Track3D& trk) const
{
	size_t h = 0;
	while (h < trk.size())
	{
		if (trk[h]->IsEnabled()) ++h;
		else (trk.release_at(h));
	}
}

// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildSegment(
	const std::vector< art::Ptr<recob::Hit> >& hits_1,
	const std::vector< art::Ptr<recob::Hit> >& hits_2) const
{
	pma::Track3D* trk = new pma::Track3D();
	trk->SetEndSegWeight(0.001F);
	trk->AddHits(hits_1);
	trk->AddHits(hits_2);

	if (trk->HasTwoViews() &&
	   (trk->TPCs().size() == 1)) // now only in single tpc
	{
		if (!trk->Initialize(0.001F))
		{
			mf::LogWarning("ProjectionMatchingAlg") << "initialization failed, delete segment";
			delete trk;
			return 0;
		}
		trk->Optimize(0, fFineTuningEps);

		trk->SortHits();
		return trk;
	}
	else
	{
		mf::LogWarning("ProjectionMatchingAlg") << "need at least two views in single tpc";
		delete trk;
		return 0;
	}
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildSegment(
	const std::vector< art::Ptr<recob::Hit> >& hits_1,
	const std::vector< art::Ptr<recob::Hit> >& hits_2,
	const TVector3& point) const
{
	pma::Track3D* trk = buildSegment(hits_1, hits_2);

	if (trk)
	{
		double dfront = pma::Dist2(trk->front()->Point3D(), point);
		double dback = pma::Dist2(trk->back()->Point3D(), point);
		if (dfront > dback) trk->Flip();

		trk->Nodes().front()->SetPoint3D(point);
		trk->Nodes().front()->SetFrozen(true);
		trk->Optimize(0, fFineTuningEps);

		trk->SortHits();
	}
	return trk;
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::buildSegment(
	const std::vector< art::Ptr<recob::Hit> >& hits,
	const TVector3& point) const
{
	pma::Track3D* trk = buildSegment(hits);

	if (trk)
	{
		double dfront = pma::Dist2(trk->front()->Point3D(), point);
		double dback = pma::Dist2(trk->back()->Point3D(), point);
		if (dfront > dback) trk->Flip();

		trk->Nodes().front()->SetPoint3D(point);
		trk->Nodes().front()->SetFrozen(true);
		trk->Optimize(0, fFineTuningEps);

		trk->SortHits();
	}
	return trk;
}
// ------------------------------------------------------

pma::Track3D* pma::ProjectionMatchingAlg::extendTrack(
	const pma::Track3D& trk,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	bool add_nodes) const
{
	pma::Track3D* copy = new pma::Track3D(trk);
	copy->AddHits(hits);

	mf::LogVerbatim("ProjectionMatchingAlg") << "ext. track size: " << copy->size()
		<< "  #coll:" << copy->NHits(geo::kZ)
		<< "  #ind2:" << copy->NHits(geo::kV)
		<< "  #ind1:" << copy->NHits(geo::kU);

	if (add_nodes)
	{
		size_t nSegments = getSegCount(copy->size());
		int nNodes = nSegments - copy->Nodes().size() + 1; // n nodes to add
		if (nNodes < 0) nNodes = 0;

		if (nNodes)
		{
			mf::LogVerbatim("ProjectionMatchingAlg") << "  add " << nNodes << " nodes";
			copy->Optimize(nNodes, fOptimizationEps);
		}
	}
	double g = copy->Optimize(0, fFineTuningEps);
	mf::LogVerbatim("ProjectionMatchingAlg") << "  reopt done, g = " << g;

	return copy;
}
// ------------------------------------------------------

bool pma::ProjectionMatchingAlg::chkEndpointHits(
	int wire, int wdir, double drift_x, int view,
	unsigned int tpc, unsigned int cryo,
	const pma::Track3D& trk,
	const std::vector< art::Ptr<recob::Hit> >& hits) const
{
	size_t nCloseHits = 0;
	int forwWires = 3, backWires = -1;
	double xMargin = 0.4;
	for (auto h : hits)
		if ((view == (int)h->WireID().Plane) &&
		    (tpc == h->WireID().TPC) &&
			(cryo == h->WireID().Cryostat))
	{
		bool found = false;
		for (size_t ht = 0; ht < trk.size(); ht++)
			if (trk[ht]->Hit2DPtr().key() == h.key())
		{
			found = true; break;
		}
		if (found) continue;

		int dw = wdir * (h->WireID().Wire - wire);
		if ((dw <= forwWires) && (dw >= backWires))
		{
			double x = fDetProp->ConvertTicksToX(h->PeakTime(), view, tpc, cryo);
			if (fabs(x - drift_x) < xMargin) nCloseHits++;
		}
	}
	if (nCloseHits > 1) return false;
	else return true;
}

bool pma::ProjectionMatchingAlg::addEndpointRef(pma::Track3D& trk,
	const std::map< unsigned int, std::vector< art::Ptr<recob::Hit> > >& hits,
	std::pair<int, int> const * wires, double const * xPos,
	unsigned int tpc, unsigned int cryo) const
{
	double x = 0.0, y = 0.0, z = 0.0;
	std::vector< std::pair<int, unsigned int> > wire_view;
	for (unsigned int i = 0; i < 3; i++)
		if (wires[i].first >= 0)
	{
		const auto hiter = hits.find(i);
		if (hiter != hits.end())
		{
			if (chkEndpointHits(wires[i].first, wires[i].second, xPos[i], i, tpc, cryo, trk, hiter->second))
			{
				x += xPos[i];
				wire_view.push_back(std::pair<int, unsigned int>(wires[i].first, i));
			}
		}
	}
	if (wire_view.size() > 1)
	{
		x /= wire_view.size();
		fGeom->IntersectionPoint(
			wire_view[0].first, wire_view[1].first,
			wire_view[0].second, wire_view[1].second,
			cryo, tpc, y, z);
		trk.AddRefPoint(x, y, z);
		mf::LogVerbatim("ProjectionMatchingAlg") << "trk tpc:" << tpc << " size:" << trk.size()
			<< " add ref.point (" << x << "; " << y << "; " << z << ")";
		return true;
	}
	else
	{
		mf::LogVerbatim("ProjectionMatchingAlg") << "trk tpc:" << tpc << " size:" << trk.size()
			<< " wire-plane-parallel track, but need two clean views of endpoint";
		return false;
	}
}

void pma::ProjectionMatchingAlg::guideEndpoints(
	pma::Track3D& trk, const std::map< unsigned int, std::vector< art::Ptr<recob::Hit> > >& hits) const
{
	unsigned int tpc = trk.FrontTPC(), cryo = trk.FrontCryo();
	if ((tpc != trk.BackTPC()) || (cryo != trk.BackCryo()))
	{
		mf::LogWarning("ProjectionMatchingAlg") << "Please, apply before TPC stitching.";
		return;
	}

	const double maxCosXZ = 0.992546; // 7 deg

	pma::Segment3D* segFront = trk.Segments().front();
	if (trk.Segments().size() > 2)
	{
		pma::Segment3D* segFront1 = trk.Segments()[1];
		if ((segFront->Length() < 0.8) && (segFront1->Length() > 5.0))
			segFront = segFront1;
	}
	TVector3 dirFront = segFront->GetDirection3D();
	TVector3 dirFrontXZ(dirFront.X(), 0., dirFront.Z());
	dirFrontXZ *= 1.0 / dirFrontXZ.Mag();

	pma::Segment3D* segBack = trk.Segments().back();
	if (trk.Segments().size() > 2)
	{
		pma::Segment3D* segBack1 = trk.Segments()[trk.Segments().size() - 2];
		if ((segBack->Length() < 0.8) && (segBack1->Length() > 5.0))
			segBack = segBack1;
	}
	TVector3 dirBack = segBack->GetDirection3D();
	TVector3 dirBackXZ(dirBack.X(), 0., dirBack.Z());
	dirBackXZ *= 1.0 / dirBackXZ.Mag();

	if ((fabs(dirFrontXZ.Z()) < maxCosXZ) && (fabs(dirBackXZ.Z()) < maxCosXZ))
	{
		return; // front & back are not parallel to wire planes => exit
	}

	unsigned int nPlanesFront = 0, nPlanesBack = 0;
	std::pair<int, int> wiresFront[3], wiresBack[3]; // wire index; index direction
	double xFront[3], xBack[3];

	for (unsigned int i = 0; i < 3; i++)
	{
		bool frontPresent = false, backPresent = false;
		if (fGeom->TPC(tpc, cryo).HasPlane(i))
		{
			int idxFront0 = trk.NextHit(-1, i);
			int idxBack0 = trk.PrevHit(trk.size(), i);
			if ((idxFront0 >= 0) && (idxBack0 >= 0))
			{
				int idxFront1 = trk.NextHit(idxFront0, i);
				int idxBack1 = trk.PrevHit(idxBack0, i);
				if ((idxFront1 >= 0) && (idxBack1 >= 0))
				{
					int wFront0 = trk[idxFront0]->Wire();
					int wBack0 = trk[idxBack0]->Wire();

					int wFront1 = trk[idxFront1]->Wire();
					int wBack1 = trk[idxBack1]->Wire();

					wiresFront[i].first = wFront0;
					wiresFront[i].second = wFront0 - wFront1;
					xFront[i] = fDetProp->ConvertTicksToX(trk[idxFront0]->PeakTime(), i, tpc, cryo);

					wiresBack[i].first = wBack0;
					wiresBack[i].second = wBack0 - wBack1;
					xBack[i] = fDetProp->ConvertTicksToX(trk[idxBack0]->PeakTime(), i, tpc, cryo);

					if (wiresFront[i].second)
					{
						if (wiresFront[i].second > 0) wiresFront[i].second = 1;
						else wiresFront[i].second = -1;

						frontPresent = true;
						nPlanesFront++;
					}

					if (wiresBack[i].second)
					{
						if (wiresBack[i].second > 0) wiresBack[i].second = 1;
						else wiresBack[i].second = -1;

						backPresent = true;
						nPlanesBack++;
					}
				}
			}
		}
		if (!frontPresent) { wiresFront[i].first = -1; }
		if (!backPresent) { wiresBack[i].first = -1; }
	}

	bool refAdded = false;
	if ((nPlanesFront > 1) && (fabs(dirFrontXZ.Z()) >= maxCosXZ))
	{
		refAdded |= addEndpointRef(trk, hits, wiresFront, xFront, tpc, cryo);
	}

	if ((nPlanesBack > 1) && (fabs(dirBackXZ.Z()) >= maxCosXZ))
	{
		refAdded |= addEndpointRef(trk, hits, wiresBack, xBack, tpc, cryo);
	}
	if (refAdded)
	{
		mf::LogVerbatim("ProjectionMatchingAlg") << "guide wire-plane-parallel track endpoints";
		double g = trk.Optimize(0, 0.1 * fFineTuningEps);
		mf::LogVerbatim("ProjectionMatchingAlg") << "  done, g = " << g;
	}
}
// ------------------------------------------------------

std::vector< pma::Hit3D* > pma::ProjectionMatchingAlg::trimTrackToVolume(
	pma::Track3D& trk, TVector3 p0, TVector3 p1) const
{
	std::vector< pma::Hit3D* > trimmedHits;



	return trimmedHits;
}
// ------------------------------------------------------

bool pma::ProjectionMatchingAlg::alignTracks(pma::Track3D& first, pma::Track3D& second) const
{
	unsigned int k = 0;
	double distFF = pma::Dist2(first.front()->Point3D(), second.front()->Point3D());
	double dist = distFF;

	double distFB = pma::Dist2(first.front()->Point3D(), second.back()->Point3D());
	if (distFB < dist) { k = 1; dist = distFB; }

	double distBF = pma::Dist2(first.back()->Point3D(), second.front()->Point3D());
	if (distBF < dist) { k = 2; dist = distBF; }

	double distBB = pma::Dist2(first.back()->Point3D(), second.back()->Point3D());
	if (distBB < dist) { k = 3; dist = distBB; }

	switch (k) // flip to get dst end before src start, do not merge if track's order reversed
	{
		case 0:	first.Flip(); break;
		case 1: mf::LogError("PMAlgTrackMaker") << "Tracks in reversed order."; return false;
		case 2: break;
		case 3: second.Flip(); break;
		default: mf::LogError("PMAlgTrackMaker") << "Should never happen."; return false;
	}
	return true;
}

void pma::ProjectionMatchingAlg::mergeTracks(pma::Track3D& dst, pma::Track3D& src, bool reopt) const
{
	if (!alignTracks(dst, src)) return;

	unsigned int tpc = src.FrontTPC();
	unsigned int cryo = src.FrontCryo();
	double lmean = dst.Length() / (dst.Nodes().size() - 1);
	if ((pma::Dist2(
			dst.Nodes().back()->Point3D(),
			src.Nodes().front()->Point3D()) > 0.5 * lmean) ||
	    (tpc != dst.BackTPC()) || (cryo != dst.BackCryo()))
	{
		dst.AddNode(src.Nodes().front()->Point3D(), tpc, cryo);
		if (src.Nodes().front()->IsFrozen())
			dst.Nodes().back()->SetFrozen(true);
	}
	for (size_t n = 1; n < src.Nodes().size(); n++)
	{
		pma::Node3D* node = src.Nodes()[n];

		dst.AddNode(src.Nodes()[n]->Point3D(),
		            node->TPC(), node->Cryo());

		if (node->IsFrozen())
			dst.Nodes().back()->SetFrozen(true);
	}
	for (size_t h = 0; h < src.size(); h++)
	{
		dst.push_back(src[h]->Hit2DPtr());
	}
	if (reopt)
	{
		double g = dst.Optimize(0, fFineTuningEps);
		mf::LogVerbatim("ProjectionMatchingAlg") << "  reopt after merging done, g = " << g;
	}
	else
	{
		dst.MakeProjection();
	}

	dst.SortHits();
	dst.ShiftEndsToHits();

	dst.MakeProjection();
	dst.SortHits();
}
// ------------------------------------------------------

double pma::ProjectionMatchingAlg::selectInitialHits(pma::Track3D& trk, unsigned int view, unsigned int* nused) const
{
	for (size_t i = 0; i < trk.size(); i++)
	{
		pma::Hit3D* hit = trk[i];
		if (hit->View2D() == view)
		{
			if ((hit->GetDistToProj() > 0.5) || // more than 0.5cm away away from the segment
			    (hit->GetSegFraction() < -1.0)) // projects before segment start (to check!!!)
				hit->TagOutlier(true);
			else hit->TagOutlier(false);
		}
	}

	unsigned int nhits = 0;
	double last_x, dx = 0.0, last_q, dq = 0.0, dqdx = 0.0;
	int ih = trk.NextHit(-1, view);

	pma::Hit3D* hit = trk[ih];
	pma::Hit3D* lastHit = hit;

	if ((ih >= 0) && (ih < (int)trk.size()))
	{
		hit->TagOutlier(true);

		ih = trk.NextHit(ih, view);
		while ((dx < 2.5) && (ih >= 0) && (ih < (int)trk.size()))
		{
			hit = trk[ih];

			if (abs(hit->Wire() - lastHit->Wire()) > 2)
				break; // break on gap in wire direction

			last_x = trk.HitDxByView(ih, view);
			last_q = hit->SummedADC();
			if (dx + last_x < 3.0)
			{
				dq += last_q;
				dx += last_x;
				nhits++;
			}
			else break;

			lastHit = hit;
			ih = trk.NextHit(ih, view);
		}
		while ((ih >= 0) && (ih < (int)trk.size()))
		{
			hit = trk[ih];
			hit->TagOutlier(true);
			ih = trk.NextHit(ih, view);
		}
	}
	else { mf::LogError("ProjectionMatchingAlg") << "Initial part selection failed."; }

	if (!nhits) { mf::LogError("ProjectionMatchingAlg") << "Initial part too short to select useful hits."; }

	if (dx > 0.0) dqdx = dq / dx;

  if (nused) (*nused) = nhits;

	//std::cout << "nhits=" << nhits << ", dq=" << dq << ", dx=" << dx << std::endl;
	return dqdx;
}
// ------------------------------------------------------

void pma::ProjectionMatchingAlg::setTrackTag(pma::Track3D& trk) const
{
	double length = trk.Length();
	double meanAngle = trk.GetMeanAng();

	unsigned int cryo = trk.FrontCryo();
	unsigned int tpc = trk.FrontTPC();

	double sizeU = 0.0, sizeV = 0.0, sizeZ = 0.0, maxSize = 0.0;
	if (fGeom->TPC(tpc, cryo).HasPlane(geo::kU))
	{
		sizeU = fabs(trk.Nodes().front()->Projection2D(geo::kU).X() - trk.Nodes().back()->Projection2D(geo::kU).X());
		maxSize = sizeU;
	}
	if (fGeom->TPC(tpc, cryo).HasPlane(geo::kV))
	{
		sizeV = fabs(trk.Nodes().front()->Projection2D(geo::kV).X() - trk.Nodes().back()->Projection2D(geo::kV).X());
		if (sizeV > maxSize) maxSize = sizeV;
	}
	if (fGeom->TPC(tpc, cryo).HasPlane(geo::kZ))
	{
		sizeZ = fabs(trk.Nodes().front()->Projection2D(geo::kZ).X() - trk.Nodes().back()->Projection2D(geo::kZ).X());
		if (sizeZ > maxSize) maxSize = sizeZ;
	}

	// skip views with short wire span in projection (EM tagging)
	double hmse, mse_sum = 0.0;
	size_t nEnabled, n_sum = 0;
	if (sizeU > 0.5 * maxSize)
	{
		hmse = trk.GetMse(geo::kU); nEnabled = trk.NEnabledHits(geo::kU);
		mse_sum += nEnabled * hmse; n_sum += nEnabled;
	}
	if (sizeV > 0.5 * maxSize)
	{
		hmse = trk.GetMse(geo::kV); nEnabled = trk.NEnabledHits(geo::kV);
		mse_sum += nEnabled * hmse; n_sum += nEnabled;
	}
	if (sizeZ > 0.5 * maxSize)
	{
		hmse = trk.GetMse(geo::kZ); nEnabled = trk.NEnabledHits(geo::kZ);
		mse_sum += nEnabled * hmse; n_sum += nEnabled;
	}
	if (n_sum > 0) mse_sum /= n_sum;
	else mse_sum = 0.0;

	if ( // (length < 80.0) &&  // tag only short tracks as EM shower-like
	     // ((mse_sum > 0.0005 * length) || (meanAngle < 3.0))) // select only very chaotic EM parts
		((mse_sum > 0.0001 * length) || (meanAngle < 3.02)))
	{
		trk.SetTag(pma::Track3D::kEmLike);
	}
	else
	{
		trk.SetTag(pma::Track3D::kTrackLike);
	}
}
// ------------------------------------------------------


