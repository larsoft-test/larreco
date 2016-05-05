void PMAlgTrackMaker::produce(art::Event& evt)
{

	fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();

	reset(evt); // set default values, clear containers at the beginning of each event

	pma::TrkCandidateColl result;

	std::unique_ptr< std::vector< recob::Track > > tracks(new std::vector< recob::Track >);
	std::unique_ptr< std::vector< recob::SpacePoint > > allsp(new std::vector< recob::SpacePoint >);
	std::unique_ptr< std::vector< recob::Vertex > > vtxs(new std::vector< recob::Vertex >);  // interaction vertices
	std::unique_ptr< std::vector< recob::Vertex > > kinks(new std::vector< recob::Vertex >); // kinks on tracks (no new particles start in kinks)
	std::unique_ptr< std::vector< recob::Vertex > > nodes(new std::vector< recob::Vertex >); // pma nodes
	std::unique_ptr< std::vector< anab::T0 > > t0s(new std::vector< anab::T0 >);

	std::unique_ptr< art::Assns< recob::Track, recob::Hit > > trk2hit_oldway(new art::Assns< recob::Track, recob::Hit >); // ****** REMEMBER to remove when FindMany improved ******
	std::unique_ptr< art::Assns< recob::Track, recob::Hit, recob::TrackHitMeta > > trk2hit(new art::Assns< recob::Track, recob::Hit, recob::TrackHitMeta >);

	std::unique_ptr< art::Assns< recob::Track, recob::SpacePoint > > trk2sp(new art::Assns< recob::Track, recob::SpacePoint >);
	std::unique_ptr< art::Assns< recob::Track, anab::T0 > > trk2t0(new art::Assns< recob::Track, anab::T0 >);

	std::unique_ptr< art::Assns< recob::SpacePoint, recob::Hit > > sp2hit(new art::Assns< recob::SpacePoint, recob::Hit >);
	std::unique_ptr< art::Assns< recob::Vertex, recob::Track > > vtx2trk(new art::Assns< recob::Vertex, recob::Track >); // one or more tracks (particles) start in the vertex
	std::unique_ptr< art::Assns< recob::Track, recob::Vertex > > trk2kink(new art::Assns< recob::Track, recob::Vertex >); // one or more kinks on the track


	std::unique_ptr< std::vector< recob::PFParticle > > pfps(new std::vector< recob::PFParticle >);

    std::unique_ptr< art::Assns<recob::PFParticle, recob::Cluster> > pfp2clu(new art::Assns<recob::PFParticle, recob::Cluster>);
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Vertex> > pfp2vtx(new art::Assns<recob::PFParticle, recob::Vertex>);
	std::unique_ptr< art::Assns< recob::PFParticle, recob::Track > > pfp2trk(new art::Assns< recob::PFParticle, recob::Track >);

	bool sortHitsClustersOK = false;
	switch (fCluMatchingAlg)
	{
		default: // try to match from all clusters in the event
		case 1: sortHitsClustersOK = sortHits(evt); break;

		case 2: // take clusters-hit assns from PFP, keep all hits for validation
		case 3: sortHitsClustersOK = sortHitsPfp(evt); break;
	}

	if (sortHitsClustersOK)
	{
		int retCode = 0;
		switch (fCluMatchingAlg)
		{
			default:
			case 1: retCode = fromMaxCluster(evt, result); break; // try to match from all clusters in the event
			case 2: retCode = fromPfpClusterSubset(evt, result); break; // each trk matched only from clusters assigned to PFP
			case 3: retCode = fromPfpDirect(evt, result); break; // no pattern recognition, just take clusters assigned to PFP
		}
		switch (retCode)
		{
			case -2: mf::LogError("Summary") << "problem"; break;
			case -1: mf::LogWarning("Summary") << "no input"; break;
			case  0: mf::LogVerbatim("Summary") << "no tracks done"; break;
			default:
				if (retCode < 0) mf::LogVerbatim("Summary") << "unknown result";
				else if (retCode == 1) mf::LogVerbatim("Summary") << retCode << " track ready";
				else mf::LogVerbatim("Summary") << retCode << " tracks ready";
				break;
		}

		if (!result.empty()) // ok, there is something to save
		{
			size_t spStart = 0, spEnd = 0;
			double sp_pos[3], sp_err[6];
			for (size_t i = 0; i < 6; i++) sp_err[i] = 1.0;

			//double dQdxFlipThr = 0.0;
			//if (fFlipToBeam) dQdxFlipThr = 0.4;

            // use the following to create PFParticle <--> Track associations;
			// note: these are assns to existing PFParticles, that are used for CluMatchingAlg = 2 or 3.
            std::map< size_t, std::vector< art::Ptr<recob::Track> > > pfPartToTrackVecMap;

			if (fFlipToBeam) result.flipTreesToCoordinate(2);        // flip the tracks / trees to the beam direction (Z)
			else if (fFlipDownward) result.flipTreesToCoordinate(1); // flip the tracks / trees to point downward (-Y)

			if (fAutoFlip_dQdx) result.flipTreesByDQdx(); // flip the tracks / trees to get best dQ/dx sequences

			tracks->reserve(result.size());
			for (fTrkIndex = 0; fTrkIndex < (int)result.size(); ++fTrkIndex)
			{
				pma::Track3D* trk = result[fTrkIndex].Track();
				if (!(trk->HasTwoViews() && (trk->Nodes().size() > 1)))
				{
					mf::LogWarning("PMAlgTrackMaker") << "Skip degenerated track, code needs to be corrected.";
					continue;
				}

				trk->SelectHits();  // just in case, set all to enabled
				unsigned int itpc = trk->FrontTPC(), icryo = trk->FrontCryo();
				if (fGeom->TPC(itpc, icryo).HasPlane(geo::kU)) trk->CompleteMissingWires(geo::kU);
				if (fGeom->TPC(itpc, icryo).HasPlane(geo::kV)) trk->CompleteMissingWires(geo::kV);
				if (fGeom->TPC(itpc, icryo).HasPlane(geo::kZ)) trk->CompleteMissingWires(geo::kZ);

				tracks->push_back(convertFrom(*trk));

				double xShift = trk->GetXShift();
				if (xShift > 0.0)
				{
					double tisk2time = 1.0; // what is the coefficient, offset?
					double t0time = tisk2time * xShift / fDetProp->GetXTicksCoefficient(trk->FrontTPC(), trk->FrontCryo());

					// TriggBits=3 means from 3d reco (0,1,2 mean something else)
					t0s->push_back(anab::T0(t0time, 0, 3, tracks->back().ID()));
					util::CreateAssn(*this, evt, *tracks, *t0s, *trk2t0, t0s->size() - 1, t0s->size());
				}

				size_t trkIdx = tracks->size() - 1; // stuff for assns:
				art::ProductID trkId = getProductID< std::vector<recob::Track> >(evt);
				art::Ptr<recob::Track> trkPtr(trkId, trkIdx, evt.productGetter(trkId));

				// which idx from start, except disabled, really....
				unsigned int hIdxs[trk->size()];
				for (size_t h = 0, cnt = 0; h < trk->size(); h++)
				{
					if ((*trk)[h]->IsEnabled()) hIdxs[h] = cnt++;
					else hIdxs[h] = 0;
				}

				art::PtrVector< recob::Hit > sp_hits;
				spStart = allsp->size();
				for (int h = trk->size() - 1; h >= 0; h--)
				{
					pma::Hit3D* h3d = (*trk)[h];
					if (!h3d->IsEnabled()) continue;

					recob::TrackHitMeta metadata(hIdxs[h], h3d->Dx());
					trk2hit->addSingle(trkPtr, h3d->Hit2DPtr(), metadata);
					trk2hit_oldway->addSingle(trkPtr, h3d->Hit2DPtr()); // ****** REMEMBER to remove when FindMany improved ******

					double hx = h3d->Point3D().X() + xShift;
					double hy = h3d->Point3D().Y();
					double hz = h3d->Point3D().Z();

					if ((h == 0) || (sp_pos[0] != hx) || (sp_pos[1] != hy) || (sp_pos[2] != hz))
					{
						if (sp_hits.size()) // hits assigned to the previous sp
						{
							util::CreateAssn(*this, evt, *allsp, sp_hits, *sp2hit);
							sp_hits.clear();
						}
						sp_pos[0] = hx; sp_pos[1] = hy; sp_pos[2] = hz;
						allsp->push_back(recob::SpacePoint(sp_pos, sp_err, 1.0));
					}
					sp_hits.push_back(h3d->Hit2DPtr());
				}

				if (sp_hits.size()) // hits assigned to the last sp
				{
					util::CreateAssn(*this, evt, *allsp, sp_hits, *sp2hit);
				}
				spEnd = allsp->size();

				if (spEnd > spStart) util::CreateAssn(*this, evt, *tracks, *allsp, *trk2sp, spStart, spEnd);

                // if there is a PFParticle collection then recover PFParticle and add info to map
                if (!fMakePFPs && (result[fTrkIndex].Key() > -1))
                {
                    size_t trackIdx = tracks->size() - 1;
                    art::ProductID trackId = getProductID< std::vector<recob::Track> >(evt);
                    art::Ptr<recob::Track> trackPtr(trackId, trackIdx, evt.productGetter(trackId));
                    pfPartToTrackVecMap[result[fTrkIndex].Key()].push_back(trackPtr);
                }
			}

			auto pfpid = getProductID< std::vector<recob::PFParticle> >(evt);
			auto vid = getProductID< std::vector<recob::Vertex> >(evt);
			auto kid = getProductID< std::vector<recob::Vertex> >(evt, kKinksName);
			auto const* kinkGetter = evt.productGetter(kid);

			auto tid = getProductID< std::vector<recob::Track> >(evt);
			auto const* trkGetter = evt.productGetter(tid);

			auto vsel = fPMAlgVertexing.getVertices(result, fSaveOnlyBranchingVtx); // vtx pos's with vector of connected track idxs
			auto ksel = fPMAlgVertexing.getKinks(result); // pairs of kink position - associated track idx 
			std::map< size_t, art::Ptr<recob::Vertex> > frontVtxs; // front vertex ptr for each track index

			if (fRunVertexing) // save vertices and vtx-trk assns
			{
				double xyz[3];
				for (auto const & v : vsel)
				{
					xyz[0] = v.first.X(); xyz[1] = v.first.Y(); xyz[2] = v.first.Z();
					mf::LogVerbatim("Summary")
						<< "  vtx:" << xyz[0] << ":" << xyz[1] << ":" << xyz[2]
						<< "  (" << v.second.size() << " tracks)";

					size_t vidx = vtxs->size();
					vtxs->push_back(recob::Vertex(xyz, vidx));

					art::Ptr<recob::Vertex> vptr(vid, vidx, evt.productGetter(vid));
					if (vptr.isNull()) mf::LogWarning("PMAlgTrackMaker") << "Vertex ptr is null.";
					if (!v.second.empty())
					{
						for (const auto & vEntry : v.second)
						{
							size_t tidx = vEntry.first;
							bool isFront = vEntry.second;

							if (isFront) frontVtxs[tidx] = vptr; // keep ptr of the front vtx

							art::Ptr<recob::Track> tptr(tid, tidx, trkGetter);
							vtx2trk->addSingle(vptr, tptr);
						}
					}
					else mf::LogWarning("PMAlgTrackMaker") << "No tracks found at this vertex.";
				}
				mf::LogVerbatim("Summary") << vtxs->size() << " vertices ready";

				for (auto const & k : ksel)
				{
					xyz[0] = k.first.X(); xyz[1] = k.first.Y(); xyz[2] = k.first.Z();
					mf::LogVerbatim("Summary") << "  kink:" << xyz[0] << ":" << xyz[1] << ":" << xyz[2];

					size_t kidx = kinks->size();
					size_t tidx = k.second; // track idx on which this kink was found

					kinks->push_back(recob::Vertex(xyz, tidx)); // save index of track (will have color of trk in evd)

					art::Ptr<recob::Track> tptr(tid, tidx, trkGetter);
					art::Ptr<recob::Vertex> kptr(kid, kidx, kinkGetter);
					trk2kink->addSingle(tptr, kptr);
				}
				mf::LogVerbatim("Summary") << ksel.size() << " kinks ready";
			}

			if (fSavePmaNodes)
			{
				double xyz[3];
				for (size_t t = 0; t < result.size(); ++t)
				{
					auto const & trk = *(result[t].Track());
					for (auto const * node : trk.Nodes())
					{
						xyz[0] = node->Point3D().X(); xyz[1] = node->Point3D().Y(); xyz[2] = node->Point3D().Z();
						nodes->push_back(recob::Vertex(xyz, t));
					}
				}
			}

			if (fMakePFPs) // create new collection of PFParticles to save hierarchy as found by this module
			{
				// first particle, to be replaced with nu reco when possible
				pfps->emplace_back(0, 0, 0, std::vector< size_t >());

				result.setParentDaughterConnections();
				for (size_t t = 0; t < result.size(); ++t)
				{
					size_t parentIdx = 0;
					if (result[t].Parent() >= 0) parentIdx = (size_t)result[t].Parent() + 1;

					std::vector< size_t > daughterIdxs;
					for (size_t idx : result[t].Daughters()) daughterIdxs.push_back(idx + 1);

					size_t pfpidx = pfps->size();
					pfps->emplace_back(0, pfpidx, parentIdx, daughterIdxs);

					art::Ptr<recob::PFParticle> pfpptr(pfpid, pfpidx, evt.productGetter(pfpid));
					art::Ptr<recob::Track> tptr(tid, t, trkGetter);
					pfp2trk->addSingle(pfpptr, tptr);

					if (fRunVertexing) // vertexing was used, so add assns to front vertex of each particle
					{
						art::Ptr<recob::Vertex> vptr = frontVtxs[t];
						if (!vptr.isNull()) pfp2vtx->addSingle(pfpptr, vptr);
						else mf::LogWarning("PMAlgTrackMaker") << "Front vertex for PFParticle is missing.";
					}
				}
				mf::LogVerbatim("Summary") << pfps->size() << " PFParticles created";
			}
			else
			{
	            // if we have used existing PFParticles then do the associations here
    	        if (!pfPartToTrackVecMap.empty())
    	        {
					art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
					evt.getByLabel(fCluModuleLabel, pfParticleHandle);
    	            for (const auto & pfParticleItr : pfPartToTrackVecMap)
    	            {
    	                art::Ptr<recob::PFParticle> pfParticle(pfParticleHandle, pfParticleItr.first);
    	                mf::LogVerbatim("PMAlgTrackMaker") << "PFParticle key: " << pfParticle.key()
							<< ", self: " << pfParticle->Self() << ", #tracks: " << pfParticleItr.second.size();

    	                if (!pfParticle.isNull()) util::CreateAssn(*this, evt, pfParticle, pfParticleItr.second, *pfp2trk);
						else mf::LogError("PMAlgTrackMaker") << "Error in PFParticle lookup, pfparticle index: "
							<< pfParticleItr.first << ", key: " << pfParticle.key();
    	            }
    	        }
			}

			// data prods done, delete all pma::Track3D's
			for (auto t : result.tracks()) t.DeleteTrack();
		}
	}
	else mf::LogError("PMAlgTrackMaker") << "Hits not found in the event.";

	evt.put(std::move(tracks));
	evt.put(std::move(allsp));
	evt.put(std::move(vtxs));
	evt.put(std::move(kinks), kKinksName);
	evt.put(std::move(nodes), kNodesName);
	evt.put(std::move(t0s));

	evt.put(std::move(trk2hit_oldway)); // ****** REMEMBER to remove when FindMany improved ******
	evt.put(std::move(trk2hit));
	evt.put(std::move(trk2sp));
	evt.put(std::move(trk2t0));

	evt.put(std::move(sp2hit));
	evt.put(std::move(vtx2trk));
	evt.put(std::move(trk2kink), kKinksName);

	evt.put(std::move(pfps));
	evt.put(std::move(pfp2clu));
	evt.put(std::move(pfp2vtx));
	evt.put(std::move(pfp2trk));
}
// ------------------------------------------------------
// ------------------------------------------------------
// ------------------------------------------------------

int PMAlgTrackMaker::fromMaxCluster(const art::Event& evt, pma::TrkCandidateColl & result)
{
	art::Handle< std::vector<recob::Cluster> > cluListHandle;
	if (evt.getByLabel(fCluModuleLabel, cluListHandle))
	{
		initial_clusters.clear();
		tried_clusters.clear();
		used_clusters.clear();

        std::vector< art::Ptr<recob::Cluster> > clusterVec;
        art::fill_ptr_vector(clusterVec, cluListHandle); // use all clusters

		tpc_track_map tracks; // track parts in tpc's

		// find reasonably large parts
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			fromMaxCluster_tpc(tracks[tpc_iter->TPC], clusterVec, fMinSeedSize1stPass, tpc_iter->TPC, tpc_iter->Cryostat);
		}

		// loop again to find small things
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			fromMaxCluster_tpc(tracks[tpc_iter->TPC], clusterVec, fMinSeedSize2ndPass, tpc_iter->TPC, tpc_iter->Cryostat);
		}

		// try correcting track ends:
		//   - 3D ref.points for clean endpoints of wire-plae parallel tracks
		//   - single-view sections spuriously merged on 2D clusters level
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			guideEndpoints(tracks[tpc_iter->TPC]);
			reassignSingleViewEnds(tracks[tpc_iter->TPC], clusterVec);
		}

		if (fMergeWithinTPC)
		{
			for (auto tpc_iter = fGeom->begin_TPC_id();
			          tpc_iter != fGeom->end_TPC_id();
			          tpc_iter++)
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "Merge co-linear tracks within TPC " << tpc_iter->TPC << ".";
				while (mergeCoLinear(tracks[tpc_iter->TPC]))
				{
					mf::LogVerbatim("PMAlgTrackMaker") << "  found co-linear tracks";
				}
			}
		}

		if (fStitchBetweenTPCs)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Stitch co-linear tracks between TPCs.";
			mergeCoLinear(tracks);
		}

		for (auto const & tpc_entry : tracks) // put tracks in the single collection
			for (auto & trk : tpc_entry.second.tracks())
				if (trk.Track()->HasTwoViews() && (trk.Track()->Nodes().size() > 1))
		{
			fProjectionMatchingAlg.setTrackTag(*(trk.Track())); // tag EM-like tracks
			result.push_back(trk);
		}

		if (fRunVertexing)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Vertex finding / track-vertex reoptimization.";
			fPMAlgVertexing.run(result);

			//reassignSingleViewEnds(result); // final check for correct hit-track assignments
		}

		if (fMatchT0inAPACrossing)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Find co-linear APA-crossing tracks with any T0.";
			matchCoLinearAnyT0(result);
		}

		listUsedClusters(clusterVec);
	}
	else
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "no clusters";
		return -1;
	}

	return result.size();
}


void PMAlgTrackMaker::fromMaxCluster_tpc(
	pma::TrkCandidateColl & result,
	const std::vector< art::Ptr<recob::Cluster> >& clusters,
	size_t minBuildSize, unsigned int tpc, unsigned int cryo,
	int pfParticleIdx)
{
	initial_clusters.clear();

	size_t minSizeCompl = minBuildSize / 8;  // smaller minimum required in complementary views
	if (minSizeCompl < 2) minSizeCompl = 2;  // but at least two hits!

	int max_first_idx = 0;
	while (max_first_idx >= 0) // loop over clusters, any view, starting from the largest
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "Find max cluster...";
		max_first_idx = maxCluster(clusters, minBuildSize, geo::kUnknown, tpc, cryo); // any view
		if (max_first_idx >= 0)
		{
			geo::View_t first_view = clusters[max_first_idx]->View();

			pma::TrkCandidate candidate = matchCluster(fCluHits.at(clusters[max_first_idx].key()),
				clusters, minSizeCompl, tpc, cryo, first_view, max_first_idx, pfParticleIdx);

			if (candidate.IsGood()) result.push_back(candidate);
		}
		else mf::LogVerbatim("PMAlgTrackMaker") << "small clusters only";
	}

	initial_clusters.clear();
}

int PMAlgTrackMaker::matchCluster(const pma::TrkCandidate& trk,
	const std::vector< art::Ptr<recob::Cluster> >& clusters,
	size_t minSize, double fraction,
	unsigned int preferedView, unsigned int testView,
	unsigned int tpc, unsigned int cryo)
{
	double f, fmax = 0.0;
	unsigned int n, max = 0;
	int idx = -1;
	for (size_t i = 0; i < clusters.size(); ++i)
	{
		unsigned int view = clusters[i]->View();
		unsigned int nhits = fCluHits.at(clusters[i].key()).size();

		if (has(used_clusters, i) ||                             // don't try already used clusters
			has(trk.Clusters(), i) ||                            // don't try clusters from this candidate
		    (view == testView) ||                                // don't use clusters from validation view
		    ((preferedView != geo::kUnknown)&&(view != preferedView)) || // only prefered view if specified
		    (nhits < minSize))                                   // skip small clusters
		    continue;

		n = fProjectionMatchingAlg.testHits(*(trk.Track()), fCluHits.at(clusters[i].key()));
		f = n / (double)nhits;
		if ((f > fraction) && (n > max))
		{
			max = n; fmax = f; idx = i;
		}
	}

	if (idx >= 0) mf::LogVerbatim("PMAlgTrackMaker") << "max matching hits: " << max << " (" << fmax << ")";
	else mf::LogVerbatim("PMAlgTrackMaker") << "no clusters to extend the track";

	return idx;
}

pma::TrkCandidate PMAlgTrackMaker::matchCluster(
	const std::vector< art::Ptr<recob::Hit> > & first_hits,
	const std::vector< art::Ptr<recob::Cluster> >& clusters,
	size_t minSizeCompl, unsigned int tpc, unsigned int cryo,
	geo::View_t first_view, int first_clu_idx, int pfParticleIdx)
{
	pma::TrkCandidate result;

	geo::View_t sec_view_a, sec_view_b;
	switch (first_view)
	{
		case geo::kU: sec_view_a = geo::kZ; sec_view_b = geo::kV; break;
		case geo::kV: sec_view_a = geo::kZ; sec_view_b = geo::kU; break;
		case geo::kZ: sec_view_a = geo::kV; sec_view_b = geo::kU; break;
		default: mf::LogError("PMAlgTrackMaker") << "Not a 2D view.";
			return result;
	}

	tried_clusters[geo::kU].clear();
	tried_clusters[geo::kV].clear();
	tried_clusters[geo::kZ].clear();

	if (first_clu_idx >= 0)
	{
		tried_clusters[first_view].push_back((size_t)first_clu_idx);
		initial_clusters.push_back((size_t)first_clu_idx);
	}

	unsigned int nFirstHits = first_hits.size();
	mf::LogVerbatim("PMAlgTrackMaker") << std::endl << "--- start new candidate ---";
	mf::LogVerbatim("PMAlgTrackMaker") << "use plane  *** " << first_view << " ***  size: " << nFirstHits;

	float x, xmax = fDetProp->ConvertTicksToX(first_hits.front()->PeakTime(), first_view, tpc, cryo), xmin = xmax;
	for (size_t j = 1; j < first_hits.size(); ++j)
	{
		x = fDetProp->ConvertTicksToX(first_hits[j]->PeakTime(), first_view, tpc, cryo);
		if (x > xmax) { xmax = x; }
		if (x < xmin) { xmin = x; }
	}

	fCandidates.clear(); // temporary set of possible solutions of the selected cluster and clusters in complementary views

	size_t imatch = 0;
	bool try_build = true;
	while (try_build) // loop over complementary views
	{
		pma::TrkCandidate candidate;
		if (first_clu_idx >= 0) candidate.Clusters().push_back((size_t)first_clu_idx);
		candidate.SetKey(pfParticleIdx);

		int idx, max_sec_a_idx, max_sec_b_idx;
		max_sec_a_idx = maxCluster(first_clu_idx, clusters, xmin, xmax, minSizeCompl, sec_view_a, tpc, cryo);
		max_sec_b_idx = maxCluster(first_clu_idx, clusters, xmin, xmax, minSizeCompl, sec_view_b, tpc, cryo);

		unsigned int nSecHitsA = 0, nSecHitsB = 0;
		if (max_sec_a_idx >= 0) nSecHitsA = fCluHits.at(clusters[max_sec_a_idx].key()).size();
		if (max_sec_b_idx >= 0) nSecHitsB = fCluHits.at(clusters[max_sec_b_idx].key()).size();

		unsigned int testView = geo::kUnknown;
		if ((nSecHitsA > nSecHitsB) && (nSecHitsA >= minSizeCompl))
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "--> " << imatch++ << " match with:";
			mf::LogVerbatim("PMAlgTrackMaker") << "use plane  *** " << sec_view_a << " ***  size: " << nSecHitsA;
			tried_clusters[sec_view_a].push_back(max_sec_a_idx);
			idx = max_sec_a_idx; testView = sec_view_b;
		}
		else if (nSecHitsB >= minSizeCompl)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "--> " << imatch++ << " match with:";
			mf::LogVerbatim("PMAlgTrackMaker") << "use plane  *** " << sec_view_b << " ***  size: " << nSecHitsB;
			tried_clusters[sec_view_b].push_back(max_sec_b_idx);
			idx = max_sec_b_idx; testView = sec_view_a;
		}
		else try_build = false;

		if (try_build)
		{
			if (!fGeom->TPC(tpc, cryo).HasPlane(testView)) testView = geo::kUnknown;

			double m0 = 0.0, v0 = 0.0;
			double mseThr = 0.15, validThr = 0.7; // cuts for a good track candidate

			candidate.Clusters().push_back(idx);
			candidate.SetTrack(fProjectionMatchingAlg.buildTrack(first_hits, fCluHits.at(clusters[idx].key())));

			if (candidate.IsValid() && // no track if hits from 2 views do not alternate
			    fProjectionMatchingAlg.isContained(*(candidate.Track()))) // sticks out of TPC's?
			{
				m0 = candidate.Track()->GetMse();
				if (m0 < mseThr) // check validation only if MSE is passing - thanks for Tracy for noticing this
					v0 = validate(*(candidate.Track()), testView);
			}
			if (candidate.Track() && (m0 < mseThr) && (v0 > validThr)) // good candidate, try to extend it
			{
				mf::LogVerbatim("PMAlgTrackMaker")
					<< "  good track candidate, MSE = " << m0 << ", v = " << v0;

				candidate.SetMse(m0);
				candidate.SetValidation(v0);
				candidate.SetGood(true);

				size_t minSize = 5;      // min size for clusters matching
				double fraction = 0.5;   // min fraction of close hits

				idx = 0;
				while (idx >= 0) // try to collect matching clusters, use **any** plane except validation
				{
					idx = matchCluster(candidate, clusters, minSize, fraction, geo::kUnknown, testView, tpc, cryo);
					if (idx >= 0)
					{
						// try building extended copy:
						//                src,        hits,    valid.plane, add nodes
						if (extendTrack(candidate, fCluHits.at(clusters[idx].key()),  testView,    true))
						{
							candidate.Clusters().push_back(idx);
						}
						else idx = -1;
					}
				}

				mf::LogVerbatim("PMAlgTrackMaker") << "merge clusters from the validation plane";
				fraction = 0.7; // only well matching the existing track

				idx = 0;
				bool extended = false;
				while ((idx >= 0) && (testView != geo::kUnknown))
				{	//                     match clusters from the plane used previously for the validation
					idx = matchCluster(candidate, clusters, minSize, fraction, testView, geo::kUnknown, tpc, cryo);
					if (idx >= 0)
					{
						// validation not checked here, no new nodes:
						if (extendTrack(candidate, fCluHits.at(clusters[idx].key()), geo::kUnknown, false))
						{
							candidate.Clusters().push_back(idx);
							extended = true;
						}
						else idx = -1;
					}
				}
				// need to calculate again only if trk was extended w/o checking validation:
				if (extended) candidate.SetValidation(validate(*(candidate.Track()), testView));
			}
			else
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "track REJECTED, MSE = " << m0 << "; v = " << v0;
				candidate.SetGood(false); // save also bad matches to avoid trying again the same pair of clusters
			}
			fCandidates.push_back(candidate);
		}
		else
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "no matching clusters";
		}
	} // end loop over complementary views

	if (fCandidates.size()) // return best candidate, release other tracks and clusters
	{
		int best_trk = -1;
		double f, max_f = 0., min_mse = 10., max_v = 0.;
		for (size_t t = 0; t < fCandidates.size(); t++)
			if (fCandidates[t].IsGood() &&
			    (fCandidates[t].Track()->Nodes().size() > 1) &&
			    fCandidates[t].Track()->HasTwoViews())
		{
			f = fProjectionMatchingAlg.twoViewFraction(*(fCandidates[t].Track()));

			if ((f > max_f) || ((f == max_f) &&
				((fCandidates[t].Validation() > max_v) || (fCandidates[t].Mse() < min_mse))))
			{
				max_f = f;
				min_mse = fCandidates[t].Mse();
				max_v = fCandidates[t].Validation();
				best_trk = t;
			}
		}

		if ((best_trk > -1) && fCandidates[best_trk].IsGood() && (max_f > fMinTwoViewFraction))
		{
			fCandidates[best_trk].Track()->ShiftEndsToHits();

			for (auto c : fCandidates[best_trk].Clusters())
				used_clusters.push_back(c);

			result = fCandidates[best_trk];
		}

		for (size_t t = 0; t < fCandidates.size(); t++)
		{
			if (int(t) != best_trk) fCandidates[t].DeleteTrack();
		}
		fCandidates.clear();
	}

	return result;
}
// ------------------------------------------------------

int PMAlgTrackMaker::maxCluster(
	const std::vector< art::Ptr<recob::Cluster> >& clusters,
	size_t min_clu_size,
	geo::View_t view, unsigned int tpc, unsigned int cryo)
{
	int idx = -1;
	size_t s_max = 0, s;

	for (size_t i = 0; i < clusters.size(); ++i)
	{
		const auto & v = fCluHits.at(clusters[i].key());

		if (!v.size() ||
		    has(used_clusters, i) ||
		    has(initial_clusters, i) ||
		    has(tried_clusters[view], i) ||
		   ((view != geo::kUnknown) && (clusters[i]->View() != view)))
		continue;

		if ((v.front()->WireID().TPC == tpc) &&
		    (v.front()->WireID().Cryostat == cryo))
		{
			s = v.size();
			if ((s >= min_clu_size) && (s > s_max))
			{
				s_max = s; idx = i;
			}
		}
	}
	return idx;
}
// ------------------------------------------------------

int PMAlgTrackMaker::maxCluster(int first_idx_tag,
	const std::vector< art::Ptr<recob::Cluster> >& clusters,
	float xmin, float xmax, size_t min_clu_size,
	geo::View_t view, unsigned int tpc, unsigned int cryo)
{
	int idx = -1;
	size_t s_max = 0, s;
	double fraction = 0.0;
	float x;

	size_t first_idx = 0;
	bool has_first = false;
	if (first_idx_tag >= 0)
	{
		first_idx = (size_t)first_idx_tag;
		has_first = true;
	}

	for (size_t i = 0; i < clusters.size(); ++i)
	{
		if (has(used_clusters, i) ||
		    has(initial_clusters, i) ||
		    has(tried_clusters[view], i) ||
		    (fCluHits.at(clusters[i].key()).size() <  min_clu_size) ||
		    (clusters[i]->View() != view)) continue;

		bool pair_checked = false;
		for (auto const & c : fCandidates.tracks())
			if (has_first && has(c.Clusters(), first_idx) && has(c.Clusters(), i))
			{
				pair_checked = true; break;
			}
		if (pair_checked) continue;
		    
		const auto & v = fCluHits.at(clusters[i].key());

		if ((v.front()->WireID().TPC == tpc) &&
		    (v.front()->WireID().Cryostat == cryo))
		{
			s = 0;
			for (size_t j = 0; j < v.size(); ++j)
			{
				x = fDetProp->ConvertTicksToX(v[j]->PeakTime(), view, tpc, cryo);
				if ((x >= xmin) && (x <= xmax)) s++;
			}

			if (s > s_max)
			{
				s_max = s; idx = i; fraction = s / (double)v.size();
			}
		}
	}
	if (fraction > 0.4) return idx;
	else return -1;
}
// ------------------------------------------------------
// ------------------------------------------------------

int PMAlgTrackMaker::fromPfpClusterSubset(const art::Event& evt, pma::TrkCandidateColl & result)
{
	bool skipPdg = true;
	if (!fTrackingSkipPdg.empty() && (fTrackingSkipPdg.front() == 0))
		skipPdg = false;

	bool selectPdg = true;
	if (!fTrackingOnlyPdg.empty() && (fTrackingOnlyPdg.front() == 0))
		selectPdg = false;

	// Code from Tracy merged with recent additions to PMA. Still to be changed in order to
	// skip not reasonalbe parts in this configuration.
    if (!fPfpClusters.empty() && !fCluHits.empty())
    {
		tpc_track_map tracks; // track parts in tpc's

		// Armed with all of this information we can begin looping through the PFParticles
		for (const auto & pfpCluEntry : fPfpClusters)
		{
			int pfPartIdx = pfpCluEntry.first;
			int pdg = fPfpPdgCodes[pfPartIdx];

			if (skipPdg && has(fTrackingSkipPdg, pdg)) continue;
			if (selectPdg && !has(fTrackingOnlyPdg, pdg)) continue;

			mf::LogVerbatim("PMAlgTrackMaker") << "Process clusters from PFP:" << pfPartIdx << ", pdg:" << pdg;

			const auto & clusterVec = pfpCluEntry.second;

			initial_clusters.clear();
			tried_clusters.clear();
			used_clusters.clear();

			size_t minBuildSize = 2;
			for (auto tpc_iter = fGeom->begin_TPC_id();
			          tpc_iter != fGeom->end_TPC_id();
			          tpc_iter++)
			{
				fromMaxCluster_tpc(tracks[tpc_iter->TPC], clusterVec, minBuildSize, tpc_iter->TPC, tpc_iter->Cryostat, pfPartIdx);
			}
   
			// used for development
			listUsedClusters(clusterVec);
		}

		// try correcting track ends:
		//   - 3D ref.points for clean endpoints of wire-plae parallel tracks
		//   - single-view sections spuriously merged on 2D clusters level
		for (auto tpc_iter = fGeom->begin_TPC_id();
		          tpc_iter != fGeom->end_TPC_id();
		          tpc_iter++)
		{
			guideEndpoints(tracks[tpc_iter->TPC]);
			reassignSingleViewEnds(tracks[tpc_iter->TPC], std::vector< art::Ptr<recob::Cluster> >());
		}

		// merge co-linear parts inside each tpc
		if (fMergeWithinTPC)
		{
			for (auto tpc_iter = fGeom->begin_TPC_id();
			          tpc_iter != fGeom->end_TPC_id();
			          tpc_iter++)
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "Merge co-linear tracks within TPC " << tpc_iter->TPC << ".";
				while (mergeCoLinear(tracks[tpc_iter->TPC]))
				{
					mf::LogVerbatim("PMAlgTrackMaker") << "  found co-linear tracks";
				}
			}
		}

		// merge co-linear parts between tpc's
		if (fStitchBetweenTPCs)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Stitch co-linear tracks between TPCs.";
			mergeCoLinear(tracks);
		}

		for (auto const & tpc_entry : tracks)
			for (auto & trk : tpc_entry.second.tracks())
				if (trk.Track()->HasTwoViews() && (trk.Track()->Nodes().size() > 1))
		{
			fProjectionMatchingAlg.setTrackTag(*(trk.Track()));
			result.push_back(trk);
		}

		if (fRunVertexing)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Vertex finding / track-vertex reoptimization.";
			fPMAlgVertexing.run(result);
		}

		if (fMatchT0inAPACrossing)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Find co-linear APA-crossing tracks with any T0.";
			matchCoLinearAnyT0(result);
		}
    }
    else
    {
        mf::LogWarning("PMAlgTrackMaker") << "no clusters, no pfparticles";
        return -1;
    }
    
    return result.size();
}

// ------------------------------------------------------

int PMAlgTrackMaker::fromPfpDirect(const art::Event& evt, pma::TrkCandidateColl & result)
{
    if (!fPfpClusters.empty() && !fCluHits.empty())
    {
			// build pm tracks
			buildTrks(result);

			guideEndpoints(result); // add 3D ref.points for clean endpoints of wire-plae parallel tracks

			if (fRunVertexing) fPMAlgVertexing.run(result);

			// build segment of shower
			buildShSeg(result);
    }
    else
    {
        mf::LogWarning("PMAlgTrackMaker") << "no clusters, no pfparticles";
        return -1;
    }
    
    return result.size();
}
// ------------------------------------------------------

void PMAlgTrackMaker::buildTrks(pma::TrkCandidateColl & result)
{
		bool skipPdg = true;
		if (!fTrackingSkipPdg.empty() && (fTrackingSkipPdg.front() == 0))
			skipPdg = false;

		bool selectPdg = true;
		if (!fTrackingOnlyPdg.empty() && (fTrackingOnlyPdg.front() == 0))
			selectPdg = false;

		for (const auto & pfpCluEntry : fPfpClusters)
		{
			int pfPartIdx = pfpCluEntry.first;
			int pdg = fPfpPdgCodes[pfPartIdx];

			if (pdg == 11) continue;
			if (skipPdg && has(fTrackingSkipPdg, pdg)) continue;
			if (selectPdg && !has(fTrackingOnlyPdg, pdg)) continue;

			mf::LogVerbatim("PMAlgTrackMaker") << "Process clusters from PFP:" << pfPartIdx << ", pdg:" << pdg;

			std::vector< art::Ptr<recob::Hit> > allHits;

			pma::TrkCandidate candidate;
			for (const auto & c : pfpCluEntry.second)
			{
				candidate.Clusters().push_back(c.key());

				allHits.reserve(allHits.size() + fCluHits.at(c.key()).size());
				for (const auto & h : fCluHits.at(c.key()))
					allHits.push_back(h);
			}
			candidate.SetKey(pfpCluEntry.first);

			candidate.SetTrack(fProjectionMatchingAlg.buildMultiTPCTrack(allHits));

			if (candidate.IsValid() &&
			    candidate.Track()->HasTwoViews() &&
			    (candidate.Track()->Nodes().size() > 1))
			{
	   			result.push_back(candidate);
			}
			else
			{
				candidate.DeleteTrack();
			}
		}
}

// ------------------------------------------------------

void PMAlgTrackMaker::buildShSeg(pma::TrkCandidateColl & result)
{
		bool skipPdg = true;
		if (!fTrackingSkipPdg.empty() && (fTrackingSkipPdg.front() == 0))
			skipPdg = false;

		bool selectPdg = true;
		if (!fTrackingOnlyPdg.empty() && (fTrackingOnlyPdg.front() == 0))
			selectPdg = false;

		for (const auto & pfpCluEntry : fPfpClusters)
		{
			int pfPartIdx = pfpCluEntry.first;
			int pdg = fPfpPdgCodes[pfPartIdx];

			if (pdg != 11) continue;
			if (skipPdg && has(fTrackingSkipPdg, pdg)) continue;
			if (selectPdg && !has(fTrackingOnlyPdg, pdg)) continue;

			mf::LogVerbatim("PMAlgTrackMaker") << "Process clusters from PFP:" << pfPartIdx << ", pdg:" << pdg;

			std::vector< art::Ptr<recob::Hit> > allHits;

			pma::TrkCandidate candidate;
			for (const auto & c : pfpCluEntry.second)
			{
				candidate.Clusters().push_back(c.key());

				allHits.reserve(allHits.size() + fCluHits.at(c.key()).size());
				for (const auto & h : fCluHits.at(c.key()))
					allHits.push_back(h);
			}

			candidate.SetKey(pfpCluEntry.first);

			mf::LogVerbatim("PMAlgTrackMaker") << "building..." << ", pdg:" << pdg;

			auto search = fPfpVtx.find(pfPartIdx);
			if (search != fPfpVtx.end()) 
			{
				candidate.SetTrack(fProjectionMatchingAlg.buildShowerSeg(allHits, fPfpVtx[pfPartIdx]));
				if (candidate.IsValid()
						&& candidate.Track()->HasTwoViews() 
						&& (candidate.Track()->Nodes().size() > 1)) 
				{
					result.push_back(candidate);
				}
				else
				{
					candidate.DeleteTrack();
				}
			}
		}
}

// ------------------------------------------------------

void PMAlgTrackMaker::listUsedClusters(const std::vector< art::Ptr<recob::Cluster> >& clusters) const
{
	mf::LogVerbatim("PMAlgTrackMaker") << std::endl << "----------- matched clusters: -----------";
	for (size_t i = 0; i < clusters.size(); ++i)
		if (has(used_clusters, i))
		{
			mf::LogVerbatim("PMAlgTrackMaker")
				<< "    tpc: " << clusters[i]->Plane().TPC
				<< ";\tview: " << clusters[i]->View()
				<< ";\tsize: " << fCluHits.at(clusters[i].key()).size();
		}
	mf::LogVerbatim("PMAlgTrackMaker") << "--------- not matched clusters: ---------";
	for (size_t i = 0; i < clusters.size(); ++i)
		if (!has(used_clusters, i))
		{
			mf::LogVerbatim("PMAlgTrackMaker")
				<< "    tpc: " << clusters[i]->Plane().TPC
				<< ";\tview: " << clusters[i]->View()
				<< ";\tsize: " << fCluHits.at(clusters[i].key()).size();
		}
	mf::LogVerbatim("PMAlgTrackMaker") << "-----------------------------------------";
}

