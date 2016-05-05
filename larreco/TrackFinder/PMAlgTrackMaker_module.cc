////////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       PMAlgTrackMaker
// Module Type: producer
// File:        PMAlgTrackMaker_module.cc
// Author:      D.Stefan (Dorota.Stefan@ncbj.gov.pl) and R.Sulej (Robert.Sulej@cern.ch), May 2015
//
// Creates 3D tracks and vertices using Projection Matching Algorithm,
// please see RecoAlg/ProjectionMatchingAlg.h for basics of the PMA algorithm and its settings.
//
// Progress:
//    May-June 2015:   track finding and validation, growing tracks by iterative merging of matching
//                     clusters, no attempts to build multi-track structures, however cosmic tracking
//                     works fine as they are sets of independent tracks
//    June-July 2015:  merging track parts within a single tpc and stitching tracks across tpc's
//    August 2015:     optimization of track-vertex structures (so 3D vertices are also produced)
//    November 2015:   use track-shower splitting at 2D level, then tag low-E EM cascades in 3D
//                     note: the splitter is not finished and not as good as we want it
//    January 2016:    output of track-vertex finding as a tree of PFParticles, refined vertexing
//                     code, put vertex at front of each track, flip whole track structures to
//                     selected, direction (beam/down/dQdx), use any pattern reco stored in
//                     PFParticles as a source of associated clusters
//    Mar-Apr 2016:    kinks finding, EM shower direction reconstructed for PFPaarticles tagged as
//                     electrons
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/PFParticle.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/TrackHitMeta.h"
#include "lardata/RecoBase/Vertex.h"
#include "lardata/RecoBase/SpacePoint.h"
#include "lardata/AnalysisBase/T0.h" 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larreco/RecoAlg/ProjectionMatchingAlg.h"
#include "larreco/RecoAlg/PMAlgVertexing.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "larreco/RecoAlg/PMAlg/PmaTrkCandidate.h"

#include "TTree.h"
#include "TMath.h"

#include <memory>

namespace trkf {

typedef std::map< size_t, std::vector<double> > dedx_map;
typedef std::map< unsigned int, std::vector< art::Ptr<recob::Hit> > > view_hitmap;
typedef std::map< unsigned int, view_hitmap > tpc_view_hitmap;
typedef std::map< unsigned int, tpc_view_hitmap > cryo_tpc_view_hitmap;

typedef std::map< size_t, pma::TrkCandidateColl > tpc_track_map;

class PMAlgTrackMaker : public art::EDProducer {
public:
  explicit PMAlgTrackMaker(fhicl::ParameterSet const & p);

  PMAlgTrackMaker(PMAlgTrackMaker const &) = delete;
  PMAlgTrackMaker(PMAlgTrackMaker &&) = delete;
  PMAlgTrackMaker & operator = (PMAlgTrackMaker const &) = delete;
  PMAlgTrackMaker & operator = (PMAlgTrackMaker &&) = delete;

  void beginJob() override;

  void reconfigure(fhicl::ParameterSet const& p) override;

  void produce(art::Event & e) override;


private:
  // *** methods to create tracks from clusters ***

  // loop over all clusters and build as much as possible to find
  // the logic implemented here for sure is not exhaustive, it was
  // checked on long cosmic tracks and low energy stopping tracks
  // and seems to be a good example how to use the algorithm
  int fromMaxCluster(const art::Event& evt, pma::TrkCandidateColl & result);

  void fromMaxCluster_tpc(
    pma::TrkCandidateColl & result,
    const std::vector< art::Ptr<recob::Cluster> >& clusters,
    size_t minBuildSize, unsigned int tpc, unsigned int cryo,
	int pfParticleIdx = -1);

  bool extendTrack(
	pma::TrkCandidate& candidate,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	unsigned int testView, bool add_nodes);

  int matchCluster(
    const pma::TrkCandidate& trk,
	const std::vector< art::Ptr<recob::Cluster> >& clusters,
	size_t minSize, double fraction,
	unsigned int preferedView, unsigned int testView,
	unsigned int tpc, unsigned int cryo);

  pma::TrkCandidate matchCluster(
	const std::vector< art::Ptr<recob::Hit> > & first_hits,
	const std::vector< art::Ptr<recob::Cluster> >& clusters,
	size_t minSizeCompl, unsigned int tpc, unsigned int cryo,
	geo::View_t first_view, int first_clu_idx, int pfParticleIdx);

	void buildTrks(pma::TrkCandidateColl & result);

	void buildShSeg(pma::TrkCandidateColl & result);

  // display what was used and what is left
  void listUsedClusters(const std::vector< art::Ptr<recob::Cluster> >& clusters) const;
  // ------------------------------------------------------

  // build tracks from clusters associated to PFParticles (use internal pattern recognition
  // on the subset of clusters selected with PFParticle)
  int fromPfpClusterSubset(const art::Event& evt, pma::TrkCandidateColl & result);
  // ------------------------------------------------------

  // build tracks from straight from clusters associated to PFParticle (no pattern recognition)
  int fromPfpDirect(const art::Event& evt, pma::TrkCandidateColl & result);
  // ------------------------------------------------------


  // ************* some common functionality **************
  void reset(const art::Event& evt);

  cryo_tpc_view_hitmap fHitMap;
  std::vector< std::vector< art::Ptr<recob::Hit> > > fCluHits;
  std::map< int, std::vector< art::Ptr<recob::Cluster> > > fPfpClusters;
  std::map< int, int > fPfpPdgCodes;
	std::map< int, art::Ptr<recob::Vertex> > fPfpVtx;
  bool sortHits(const art::Event& evt);
  bool sortHitsPfp(const art::Event& evt);

  bool has(const std::vector<size_t>& v, size_t idx) const
  {
  	for (auto c : v) if (c == idx) return true;
  	return false;
  }
  bool has(const std::vector<int>& v, int i) const
  {
  	for (auto c : v) if (c == i) return true;
  	return false;
  }

  std::vector< size_t > used_clusters, initial_clusters;
  std::map< unsigned int, std::vector<size_t> > tried_clusters;

  // temporary set of possible solutions of the selected cluster and clusters in complementary views
  pma::TrkCandidateColl fCandidates;

  int maxCluster(
    const std::vector< art::Ptr<recob::Cluster> >& clusters,
    size_t min_clu_size,
    geo::View_t view, unsigned int tpc, unsigned int cryo);
  int maxCluster(
	int first_idx_tag,
    const std::vector< art::Ptr<recob::Cluster> >& clusters,
    float tmin, float tmax, size_t min_clu_size,
    geo::View_t view, unsigned int tpc, unsigned int cryo);

  void freezeBranchingNodes(pma::TrkCandidateColl & tracks);
  void releaseAllNodes(pma::TrkCandidateColl & tracks);

  bool areCoLinear(
	pma::Track3D* trk1, pma::Track3D* trk2,
	double& dist, double& cos, bool& reverseOrder,
	double distThr, double distThrMin,
	double distProjThr,
	double cosThr);
  bool mergeCoLinear(pma::TrkCandidateColl & tracks);
  void mergeCoLinear(tpc_track_map& tracks);

  bool areCoLinear(
		double& cos3d,
		TVector3 f0, TVector3 b0, TVector3 f1, TVector3 b1,
		double distProjThr);
  void matchCoLinearAnyT0(pma::TrkCandidateColl & tracks);

  bool reassignHits(
	const std::vector< art::Ptr<recob::Hit> >& hits,
	const std::vector< art::Ptr<recob::Cluster> > & clusters,
	pma::TrkCandidateColl & tracks, size_t trk_idx, double dist2);
  bool reassignHits(
	const std::vector< art::Ptr<recob::Hit> >& hits,
	pma::TrkCandidateColl & tracks, size_t trk_idx, double dist2);

  double collectSingleViewFront(pma::Track3D & trk,
	std::vector< art::Ptr<recob::Hit> > & hits);
  double collectSingleViewEnd(pma::Track3D & trk,
	std::vector< art::Ptr<recob::Hit> > & hits);

  bool reassignSingleViewEnds(pma::TrkCandidateColl & tracks,
	const std::vector< art::Ptr<recob::Cluster> > & clusters);
  bool reassignSingleViewEnds(pma::TrkCandidateColl & tracks);
  void guideEndpoints(pma::TrkCandidateColl & tracks);

  double validate(pma::Track3D& trk, unsigned int testView);
  recob::Track convertFrom(const pma::Track3D& src);
  // ------------------------------------------------------

  art::ServiceHandle< geo::Geometry > fGeom;
  const detinfo::DetectorProperties* fDetProp;

  // ******************* tree output **********************
  int fEvNumber;        // event number
  int fTrkIndex;        // track index in the event
  int fPidTag;          // Tag: 0=trk-like, 1=cascade-like
  double fLength;       // track length
  double fHitsMse;      // MSE of hits: mean dist^2 of hit to 2D track projection
  double fSegAngMean;   // Mean segment-segment 3D angle.
  TTree* fTree_trk;     // overall info

  // ******************** fcl parameters **********************
  std::string fHitModuleLabel; // label for hits collection (used for trk validation)
  std::string fCluModuleLabel; // label for input cluster collection
  int fCluMatchingAlg;         // which algorithm for cluster association

  std::vector<int> fTrackingOnlyPdg; // make tracks only for this pdg's when using input from PFParticles
  std::vector<int> fTrackingSkipPdg; // skip tracks with this pdg's when using input from PFParticles

  bool fMakePFPs;              // output track-vertex net as a tree of PFParticles

  size_t fMinSeedSize1stPass;  // min. cluster size used to start building a track in the 1st pass
  size_t fMinSeedSize2ndPass;  // min. cluster size used to start building a track in the 2nd pass

  bool fFlipToBeam;            // set the track direction to increasing Z values
  bool fFlipDownward;          // set the track direction to decreasing Y values
  bool fAutoFlip_dQdx;         // set the track direction to increasing dQ/dx

  bool fMergeWithinTPC;          // merge witnin single TPC; finds tracks best matching by angle, with limits:
  double fMergeTransverseShift;  //   - max. transverse displacement [cm] between tracks
  double fMergeAngle;            //   - max. angle [degree] between tracks (nearest segments)

  bool fStitchBetweenTPCs;       // stitch between TPCs; finds tracks best matching by angle, with limits:
  double fStitchDistToWall;      //   - max. track endpoint distance [cm] to TPC boundary
  double fStitchTransverseShift; //   - max. transverse displacement [cm] between tracks
  double fStitchAngle;           //   - max. angle [degree] between tracks (nearest segments)

  bool fMatchT0inAPACrossing;    // match T0 of APA-crossing tracks, TPC stitching limits are used

  pma::ProjectionMatchingAlg fProjectionMatchingAlg;
  double fMinTwoViewFraction;    // ProjectionMatchingAlg parameter used also in the module

  pma::PMAlgVertexing fPMAlgVertexing;
  bool fRunVertexing;          // run vertex finding
  bool fSaveOnlyBranchingVtx;  // for debugging, save only vertices which connect many tracks
  bool fSavePmaNodes;          // for debugging, save only track nodes

  // ********** instance names (collections, assns) ************
  static const std::string kKinksName;        // kinks on tracks
  static const std::string kNodesName;        // pma nodes
};
// -------------------------------------------------------------
const std::string PMAlgTrackMaker::kKinksName = "kink";
const std::string PMAlgTrackMaker::kNodesName = "node";
// -------------------------------------------------------------

PMAlgTrackMaker::PMAlgTrackMaker(fhicl::ParameterSet const & p) :
	fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>()),
	fProjectionMatchingAlg(p.get< fhicl::ParameterSet >("ProjectionMatchingAlg")),
	fPMAlgVertexing(p.get< fhicl::ParameterSet >("PMAlgVertexing"))
{
	this->reconfigure(p);

	produces< std::vector<recob::Track> >();
	produces< std::vector<recob::SpacePoint> >();
	produces< std::vector<recob::Vertex> >(); // no instance name for interaction vertices
	produces< std::vector<recob::Vertex> >(kKinksName); // collection of kinks on tracks
	produces< std::vector<recob::Vertex> >(kNodesName); // collection of pma nodes
	produces< std::vector<anab::T0> >();

	produces< art::Assns<recob::Track, recob::Hit> >(); // ****** REMEMBER to remove when FindMany improved ******
	produces< art::Assns<recob::Track, recob::Hit, recob::TrackHitMeta> >();

	produces< art::Assns<recob::Track, recob::SpacePoint> >();
	produces< art::Assns<recob::SpacePoint, recob::Hit> >();
	produces< art::Assns<recob::Vertex, recob::Track> >(); // no instance name for assns of tracks to interaction vertices
	produces< art::Assns<recob::Track, recob::Vertex> >(kKinksName);  // assns of kinks to tracks
	produces< art::Assns<recob::Track, anab::T0> >();

	produces< std::vector<recob::PFParticle> >();
	produces< art::Assns<recob::PFParticle, recob::Cluster> >();
	produces< art::Assns<recob::PFParticle, recob::Vertex> >();
	produces< art::Assns<recob::PFParticle, recob::Track> >();
}
// ------------------------------------------------------

void PMAlgTrackMaker::beginJob()
{
	art::ServiceHandle<art::TFileService> tfs;
	fTree_trk = tfs->make<TTree>("PMAlgTrackMaker_trk", "tracks overall info");
	fTree_trk->Branch("fEvNumber", &fEvNumber, "fEvNumber/I");
	fTree_trk->Branch("fTrkIndex", &fTrkIndex, "fTrkIndex/I");
	fTree_trk->Branch("fLength", &fLength, "fLength/D");
	fTree_trk->Branch("fHitsMse", &fHitsMse, "fHitsMse/D");
	fTree_trk->Branch("fSegAngMean", &fSegAngMean, "fSegAngMean/D");
	fTree_trk->Branch("fPidTag", &fPidTag, "fPidTag/I");
}

void PMAlgTrackMaker::reconfigure(fhicl::ParameterSet const& pset)
{
	fHitModuleLabel = pset.get< std::string >("HitModuleLabel");
	fCluModuleLabel = pset.get< std::string >("ClusterModuleLabel");
	fCluMatchingAlg = pset.get< int >("CluMatchingAlg");

	fTrackingOnlyPdg = pset.get< std::vector<int> >("TrackingOnlyPdg");
	fTrackingSkipPdg = pset.get< std::vector<int> >("TrackingSkipPdg");

	fMakePFPs = pset.get< bool >("MakePFPs");

	fMinSeedSize1stPass = pset.get< size_t >("MinSeedSize1stPass");
	fMinSeedSize2ndPass = pset.get< size_t >("MinSeedSize2ndPass");

	fFlipToBeam = pset.get< bool >("FlipToBeam");
	fFlipDownward = pset.get< bool >("FlipDownward");
	fAutoFlip_dQdx = pset.get< bool >("AutoFlip_dQdx");

	fMergeWithinTPC = pset.get< bool >("MergeWithinTPC");
	fMergeTransverseShift = pset.get< double >("MergeTransverseShift");
	fMergeAngle = pset.get< double >("MergeAngle");

	fStitchBetweenTPCs = pset.get< bool >("StitchBetweenTPCs");
	fStitchDistToWall = pset.get< double >("StitchDistToWall");
	fStitchTransverseShift = pset.get< double >("StitchTransverseShift");
	fStitchAngle = pset.get< double >("StitchAngle");

	fMatchT0inAPACrossing = pset.get< bool >("MatchT0inAPACrossing");

	fProjectionMatchingAlg.reconfigure(pset.get< fhicl::ParameterSet >("ProjectionMatchingAlg"));
	fMinTwoViewFraction = pset.get< double >("ProjectionMatchingAlg.MinTwoViewFraction");

	fPMAlgVertexing.reconfigure(pset.get< fhicl::ParameterSet >("PMAlgVertexing"));
	fRunVertexing = pset.get< bool >("RunVertexing");
	fSaveOnlyBranchingVtx = pset.get< bool >("SaveOnlyBranchingVtx");
	fSavePmaNodes = pset.get< bool >("SavePmaNodes");
}

void PMAlgTrackMaker::reset(const art::Event& evt)
{
	fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();

	fHitMap.clear();
	fCluHits.clear();
	fPfpClusters.clear();
	fPfpPdgCodes.clear();
	fPfpVtx.clear();
	fEvNumber = evt.id().event();
	fTrkIndex = 0;
	fPidTag = 0;
	fLength = 0.0;
	fHitsMse = 0.0;
	fSegAngMean = 0.0;

	fPMAlgVertexing.reset();
}
// ------------------------------------------------------

recob::Track PMAlgTrackMaker::convertFrom(const pma::Track3D& src)
{
	std::vector< TVector3 > xyz, dircos;
	xyz.reserve(src.size()); dircos.reserve(src.size());

	std::vector< std::vector<double> > dst_dQdx; // [view][dQ/dx]
	dst_dQdx.push_back(std::vector<double>()); // kU
	dst_dQdx.push_back(std::vector<double>()); // kV
	dst_dQdx.push_back(std::vector<double>()); // kZ

	unsigned int cryo = src.FrontCryo();
	unsigned int tpc = src.FrontTPC();

	std::map< unsigned int, dedx_map > src_dQdx;
	if (fGeom->TPC(tpc, cryo).HasPlane(geo::kU))
	{
		src_dQdx[geo::kU] = dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kU], geo::kU);
	}
	if (fGeom->TPC(tpc, cryo).HasPlane(geo::kV))
	{
		src_dQdx[geo::kV] = dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kV], geo::kV);
	}
	if (fGeom->TPC(tpc, cryo).HasPlane(geo::kZ))
	{
		src_dQdx[geo::kZ] = dedx_map();
		src.GetRawdEdxSequence(src_dQdx[geo::kZ], geo::kZ);
	}

	TVector3 p3d;
	double xshift = src.GetXShift();
	bool has_shift = (xshift != 0.0);
	for (size_t i = 0; i < src.size(); i++)
		if (src[i]->IsEnabled())
	{
		p3d = src[i]->Point3D();
		if (has_shift) p3d.SetX(p3d.X() + xshift);
		xyz.push_back(p3d);

		if (i < src.size() - 1)
		{
			size_t j = i + 1;
			double mag = 0.0;
			TVector3 dc(0., 0., 0.);

			while ((mag == 0.0) && (j < src.size()))
			{
				dc = src[j]->Point3D();
				dc -= src[i]->Point3D();
				mag = dc.Mag();
				j++;
			}

			if (mag > 0.0) dc *= 1.0 / mag;
			else if (!dircos.empty()) dc = dircos.back();

			dircos.push_back(dc);
		}
		else dircos.push_back(dircos.back());

		double dQ = 0., dx = 0.;
		dst_dQdx[geo::kU].push_back(0.);
		dst_dQdx[geo::kV].push_back(0.);
		dst_dQdx[geo::kZ].push_back(0.);

		double dQdx;
		for (auto const& m : src_dQdx)
		{
			auto it = m.second.find(i);
			if (it != m.second.end())
			{
				dQ = it->second[5];
				dx = it->second[6];
				if (dx > 0.) dQdx = dQ/dx;
				else dQdx = 0.;

				size_t backIdx = dst_dQdx[m.first].size() - 1;
				dst_dQdx[m.first][backIdx] = dQdx;

				break;
			}
		}
	}

	fLength = src.Length();
	fHitsMse = src.GetMse();
	fSegAngMean = src.GetMeanAng();

	// 0 is track-like (long and/or very straight, well matching 2D hits);
	// 0x10000 is EM shower-like trajectory
	if (src.GetTag() == pma::Track3D::kEmLike) fPidTag = 0x10000;
	else fPidTag = 0;

	fTree_trk->Fill();

	if (xyz.size() != dircos.size())
	{
		mf::LogError("PMAlgTrackMaker") << "pma::Track3D to recob::Track conversion problem.";
	}

	return recob::Track(xyz, dircos, dst_dQdx, std::vector< double >(2, util::kBogusD), fTrkIndex + fPidTag);
}
// ------------------------------------------------------

double PMAlgTrackMaker::validate(pma::Track3D& trk, unsigned int testView)
{
	if ((trk.FirstElement()->GetDistToWall() < -3.0) ||
	    (trk.LastElement()->GetDistToWall() < -3.0))
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "first or last node too far out of its initial TPC";
		return 0.0;
	}

	if (testView != geo::kUnknown)
		mf::LogVerbatim("PMAlgTrackMaker") << "validation in plane: " << testView;
	else return 1.0;

	std::vector< art::Ptr<recob::Hit> >& hits = fHitMap[trk.FrontCryo()][trk.FrontTPC()][testView];

	// always validate (needed for disambiguation postponed to 3D step):
	return fProjectionMatchingAlg.validate(trk, hits, testView);

	// in case of usual disambig on hit level one may validate only if there are at least a few hits:
/*	if (hits.size() > 10) return fProjectionMatchingAlg.validate(trk, hits, testView);
	else
	{
		mf::LogWarning("PMAlgTrackMaker") << "   too few hits (" << hits.size() << ")";
		return 1.0;
	}
*/
}
// ------------------------------------------------------

bool PMAlgTrackMaker::extendTrack(pma::TrkCandidate& candidate,
	const std::vector< art::Ptr<recob::Hit> >& hits,
	unsigned int testView, bool add_nodes)
{
	double m_max = 2.0 * candidate.Mse(); // max acceptable MSE value
	if (m_max < 0.05) m_max = 0.05;     // this is still good, low MSE value

	double v_min1 = 0.98 * candidate.Validation();
	double v_min2 = 0.9 * candidate.Validation();

	pma::Track3D* copy = fProjectionMatchingAlg.extendTrack(*(candidate.Track()), hits, add_nodes);
	double m1 = copy->GetMse();
	double v1 = validate(*copy, testView);

	if (((m1 < candidate.Mse()) && (v1 >= v_min2)) ||
	    ((m1 < 0.5) && (m1 <= m_max) && (v1 >= v_min1)))
	{
		mf::LogVerbatim("PMAlgTrackMaker")
			<< "  track EXTENDED, MSE = " << m1 << ", v = " << v1;
		candidate.SetTrack(copy);  // replace with the new track (deletes old one)
		copy->SortHits();          // sort hits in the new track

		candidate.SetMse(m1);      // save info
		candidate.SetValidation(v1);

		return true;
	}
	else
	{
		mf::LogVerbatim("PMAlgTrackMaker")
			<< "  track NOT extended, MSE = " << m1 << ", v = " << v1;
		delete copy;
		return false;
	}
}
// ------------------------------------------------------

void PMAlgTrackMaker::freezeBranchingNodes(pma::TrkCandidateColl & tracks)
{
	for (auto const & trk : tracks.tracks())
		for (auto node : trk.Track()->Nodes())
			if (node->IsBranching()) node->SetFrozen(true);
}
// ------------------------------------------------------

void PMAlgTrackMaker::releaseAllNodes(pma::TrkCandidateColl & tracks)
{
	for (auto const & trk : tracks.tracks())
		for (auto node : trk.Track()->Nodes())
			node->SetFrozen(false);
}
// ------------------------------------------------------

bool PMAlgTrackMaker::areCoLinear(pma::Track3D* trk1, pma::Track3D* trk2,
	double& dist, double& cos3d, bool& reverseOrder,
	double distThr, double distThrMin,
	double distProjThr,
	double cosThr)
{
	double lmax;
	double l1 = trk1->Length();
	double l2 = trk2->Length();

	if (l1 > l2) lmax = l1;
	else lmax = l2;

	double d = lmax * distThr;
	if (d < distThrMin) d = distThrMin;

	unsigned int k = 0;
	double distFF = pma::Dist2(trk1->front()->Point3D(), trk2->front()->Point3D());
	dist = distFF;

	double distFB = pma::Dist2(trk1->front()->Point3D(), trk2->back()->Point3D());
	if (distFB < dist) { k = 1; dist = distFB; }

	double distBF = pma::Dist2(trk1->back()->Point3D(), trk2->front()->Point3D());
	if (distBF < dist) { k = 2; dist = distBF; }

	double distBB = pma::Dist2(trk1->back()->Point3D(), trk2->back()->Point3D());
	if (distBB < dist) { k = 3; dist = distBB; }

	dist = sqrt(dist);
	cos3d = 0.0;

	if (dist < d)
	{
		pma::Track3D* tmp = 0;
		switch (k) // swap or flip to get trk1 end before trk2 start
		{
			case 0:	trk1->Flip(); break;
			case 1: tmp = trk1;	trk1 = trk2; trk2 = tmp; break;
			case 2: break;
			case 3: trk2->Flip(); break;
			default: mf::LogError("PMAlgTrackMaker") << "Should never happen.";
		}
		if (k == 1) reverseOrder = true;
		else reverseOrder = false;

		size_t nodeEndIdx = trk1->Nodes().size() - 1;

		TVector3 endpoint1 = trk1->back()->Point3D();
		TVector3 trk2front0 = trk2->Nodes()[0]->Point3D();
		TVector3 trk2front1 = trk2->Nodes()[1]->Point3D();
		TVector3 proj1 = pma::GetProjectionToSegment(endpoint1, trk2front0, trk2front1);
		double distProj1 = sqrt( pma::Dist2(endpoint1, proj1) );

		TVector3 endpoint2 = trk2->front()->Point3D();
		TVector3 trk1back0 = trk1->Nodes()[nodeEndIdx]->Point3D();
		TVector3 trk1back1 = trk1->Nodes()[nodeEndIdx - 1]->Point3D();
		TVector3 proj2 = pma::GetProjectionToSegment(endpoint2, trk1back1, trk1back0);
		double distProj2 = sqrt( pma::Dist2(endpoint2, proj2) );

		TVector3 dir1 = trk1->Segments().back()->GetDirection3D();
		TVector3 dir2 = trk2->Segments().front()->GetDirection3D();

		cos3d = dir1 * dir2;

		if ((cos3d > cosThr) && (distProj1 < distProjThr) && (distProj2 < distProjThr))
			return true;
		else // check if parallel to wires & colinear in 2D
		{
			const double maxCosXZ = 0.996195; // 5 deg

			TVector3 dir1_xz(dir1.X(), 0., dir1.Z());
			dir1_xz *= 1.0 / dir1_xz.Mag();

			TVector3 dir2_xz(dir2.X(), 0., dir2.Z());
			dir2_xz *= 1.0 / dir2_xz.Mag();

			if ((fabs(dir1_xz.Z()) > maxCosXZ) && (fabs(dir2_xz.Z()) > maxCosXZ))
			{
				endpoint1.SetY(0.);
				trk2front0.SetY(0.);
				trk2front1.SetY(0.);
				proj1 = pma::GetProjectionToSegment(endpoint1, trk2front0, trk2front1);
				distProj1 = sqrt( pma::Dist2(endpoint1, proj1) );

				endpoint2.SetY(0.);
				trk1back0.SetY(0.);
				trk1back1.SetY(0.);
				proj2 = pma::GetProjectionToSegment(endpoint2, trk1back1, trk1back0);
				distProj2 = sqrt( pma::Dist2(endpoint2, proj2) );
			
				double cosThrXZ = cos(0.5 * acos(cosThr));
				double distProjThrXZ = 0.5 * distProjThr;
				double cosXZ = dir1_xz * dir2_xz;
				if ((cosXZ > cosThrXZ) && (distProj1 < distProjThrXZ) && (distProj2 < distProjThrXZ))
					return true;
			}
		}
	}
	return false;
}
// ------------------------------------------------------

bool PMAlgTrackMaker::mergeCoLinear(pma::TrkCandidateColl & tracks)
{
	double distThr = 0.25;    // max gap as a fraction of the longer track length
	double distThrMin = 0.5;  // lower limit of max gap threshold [cm]

	double distProjThr = fMergeTransverseShift;
	double cosThr = cos(TMath::Pi() * fMergeAngle / 180.0);

	bool foundMerge = false;

	std::sort(tracks.tracks().begin(), tracks.tracks().end(), pma::bTrack3DLonger());

	bool r;
	double d, dmin, c, cmax, l, lbest;
	size_t t = 0, u = 0;
	while (t < tracks.size())
	{
		pma::Track3D* trk1 = tracks[t].Track();

		pma::Track3D* trk2 = 0;
		pma::Track3D* best_trk2 = 0;
		dmin = 1.0e12; cmax = 0; lbest = 0;
		for (u = t + 1; u < tracks.size(); u++)
		{
			trk2 = tracks[u].Track();
			if (areCoLinear(trk1, trk2, d, c, r, distThr, distThrMin, distProjThr, cosThr))
			{
				l = std::sqrt(pma::Dist2(trk2->front()->Point3D(), trk2->back()->Point3D()));
				if (((c > cmax) && (d < dmin + 0.5 * lbest)) ||
				    ((d < dmin) && (l > 1.5 * lbest)))
				{
					cmax = c; dmin = d;
					best_trk2 = trk2;
					lbest = l;
				}
			}
			trk2 = 0;
		}
		trk2 = best_trk2;

		if (trk2)
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "Merge track ("
				<< trk1->size() << ") with track (" << trk2->size() << ")";
			if (r)
			{
				fProjectionMatchingAlg.mergeTracks(*trk2, *trk1, true);
				tracks[t].SetTrack(trk2); // deletes old trk1
			}
			else
			{
				fProjectionMatchingAlg.mergeTracks(*trk1, *trk2, true);
				tracks[u].DeleteTrack();
			}
			tracks.erase_at(u);
			foundMerge = true;
		}
		else t++;
	}

	return foundMerge;
}
// ------------------------------------------------------

void PMAlgTrackMaker::mergeCoLinear(tpc_track_map& tracks)
{
	double distThr = 0.25;    // max gap as a fraction of the longer track length
	double distThrMin = 2.5;  // lower limit of max gap threshold [cm]

	double distProjThr = fStitchTransverseShift;
	double cosThr = cos(TMath::Pi() * fStitchAngle / 180.0);

	double wallDistThr = fStitchDistToWall;
	double dfront1, dback1, dfront2, dback2;

	//for (auto & tpc_entry : tracks) freezeBranchingNodes(tpc_entry.second);

	for (auto & tpc_entry1 : tracks)
	{
		unsigned int tpc1 = tpc_entry1.first;
		pma::TrkCandidateColl & tracks1 = tpc_entry1.second;

		size_t t = 0;
		while (t < tracks1.size())
		{
			bool r, reverse = false;
			double l, lbest = 0, d, dmin = 1.0e12, c, cmax = 0.0;
			pma::Track3D* best_trk2 = 0;
			unsigned int best_tpc = 0;
			size_t best_idx = 0;

			pma::Track3D* trk1 = tracks1[t].Track();
			dfront1 = trk1->Nodes().front()->GetDistToWall();
			dback1 = trk1->Nodes().back()->GetDistToWall();
			if ((dfront1 < wallDistThr) || (dback1 < wallDistThr))
			{
				for (auto & tpc_entry2 : tracks)
				{
					unsigned int tpc2 = tpc_entry2.first;
					if (tpc2 == tpc1) continue;

					pma::TrkCandidateColl & tracks2 = tpc_entry2.second;

					for (size_t u = 0; u < tracks2.size(); u++)
					{
						pma::Track3D* trk2 = tracks2[u].Track();
						dfront2 = trk2->Nodes().front()->GetDistToWall();
						dback2 = trk2->Nodes().back()->GetDistToWall();
						if ((dfront2 < wallDistThr) || (dback2 < wallDistThr))
						{
							if (areCoLinear(trk1, trk2, d, c, r, distThr, distThrMin, distProjThr, cosThr))
							{
								l = std::sqrt(pma::Dist2(trk2->front()->Point3D(), trk2->back()->Point3D()));
								if (((c > cmax) && (d < dmin + 0.5 * lbest)) || (0.75 * l < dmin))
								{
									cmax = c; dmin = d; lbest = l;
									best_trk2 = trk2;
									best_tpc = tpc2;
									best_idx = u;
									reverse = r;
								}
							}
						}
					}
				}
			}

			if (best_trk2)
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "Merge track ("
					<< tpc1 << ":" << tracks1.size() << ":" << trk1->size() << ") with track ("
					<< best_tpc  << ":" << tracks[best_tpc].size() << ":" << best_trk2->size() << ")";
				if (reverse)
				{
					fProjectionMatchingAlg.mergeTracks(*best_trk2, *trk1, true);
					tracks1[t].SetTrack(best_trk2);
				}
				else
				{
					fProjectionMatchingAlg.mergeTracks(*trk1, *best_trk2, true);
					tracks[best_tpc][best_idx].DeleteTrack();
				}
				tracks[best_tpc].erase_at(best_idx);
			}
			else t++;
		}
	}

	//for (auto & tpc_entry : tracks) releaseAllNodes(tpc_entry.second);
}
// ------------------------------------------------------

bool PMAlgTrackMaker::areCoLinear(double& cos3d,
		TVector3 f0, TVector3 b0, TVector3 f1, TVector3 b1,
		double distProjThr)
{
	TVector3 s0 = b0 - f0, s1 = b1 - f1;
	cos3d = s0 * s1 / (s0.Mag() * s1.Mag());

	TVector3 proj0 = pma::GetProjectionToSegment(b0, f1, b1);
	double distProj0 = sqrt( pma::Dist2(b0, proj0) );

	TVector3 proj1 = pma::GetProjectionToSegment(b1, f0, b0);
	double distProj1 = sqrt( pma::Dist2(b1, proj1) );

	double d = sqrt( pma::Dist2(b0, f1) );
	double dThr = (1 + 0.02 * d) * distProjThr;

	mf::LogVerbatim("PMAlgTrackMaker")
		<< "   dThr:" << dThr << " d0:" << distProj0 << " d1:" << distProj1 << " c:" << cos3d;

	if ((distProj0 < dThr) && (distProj1 < dThr))
		return true;
	else return false;
}
// ------------------------------------------------------

void PMAlgTrackMaker::matchCoLinearAnyT0(pma::TrkCandidateColl& tracks)
{
	double distProjThr = fStitchTransverseShift;
	double cosThr = cos(TMath::Pi() * fStitchAngle / 180.0);
	double xApaDistDiffThr = 1.0;

	for (size_t u = 0; u < tracks.size(); u++)
	{
		pma::Track3D* trk1 = tracks[u].Track();

		unsigned int tpcFront1 = trk1->FrontTPC(), cryoFront1 = trk1->FrontCryo();
		unsigned int tpcBack1 = trk1->BackTPC(), cryoBack1 = trk1->BackCryo();

		unsigned int firstPlane = 0;
		while ((firstPlane < 3) && !fGeom->TPC(tpcFront1, cryoFront1).HasPlane(firstPlane)) firstPlane++;

		double dxFront1 = trk1->front()->Point3D().X();
		dxFront1 -= fGeom->TPC(tpcFront1, cryoFront1).PlaneLocation(firstPlane)[0];

		firstPlane = 0;
		while ((firstPlane < 3) && !fGeom->TPC(tpcBack1, cryoBack1).HasPlane(firstPlane)) firstPlane++;

		double dxBack1 = trk1->back()->Point3D().X();
		dxBack1 -= fGeom->TPC(tpcBack1, cryoBack1).PlaneLocation(firstPlane)[0];

		//size_t best_idx = 0;
		pma::Track3D* best_trk2 = 0;
		double dx1 = 0.0, dx2 = 0.0, c, cmax = cosThr;
		bool reverse = false, flip1 = false, flip2 = false;
		TVector3 f0, b0, f1, b1;
		for (size_t t = u + 1; t < tracks.size(); t++)
		{
			pma::Track3D* trk2 = tracks[t].Track();
			if (trk2->GetXShift() != 0.0) continue;

			unsigned int tpcFront2 = trk2->FrontTPC(), cryoFront2 = trk2->FrontCryo();
			unsigned int tpcBack2 = trk2->BackTPC(), cryoBack2 = trk2->BackCryo();

			firstPlane = 0;
			while ((firstPlane < 3) && !fGeom->TPC(tpcFront2, cryoFront2).HasPlane(firstPlane)) firstPlane++;

			double dxFront2 = trk2->front()->Point3D().X();
			dxFront2 -= fGeom->TPC(tpcFront2, cryoFront2).PlaneLocation(firstPlane)[0];

			firstPlane = 0;
			while ((firstPlane < 3) && !fGeom->TPC(tpcBack2, cryoBack2).HasPlane(firstPlane)) firstPlane++;

			double dxBack2 = trk2->back()->Point3D().X();
			if (dxBack2 > 0.0) dxBack2 -= fGeom->TPC(tpcBack2, cryoBack2).PlaneLocation(firstPlane)[0];

			mf::LogVerbatim("PMAlgTrackMaker")
				<< "   xf1:" << dxFront1 << " xb1:" << dxBack1
				<< "   xf2:" << dxFront2 << " xb2:" << dxBack2;

			if ((cryoFront1 == cryoFront2) && (dxFront1 * dxFront2 < 0.0) &&
			    (fabs(dxFront1 + dxFront2) < xApaDistDiffThr))
			{
				if (fabs(dxFront1) < fabs(dxFront2)) dxFront2 = -dxFront1;
				else dxFront1 = -dxFront2;

				f0 = trk1->Nodes()[1]->Point3D(); f0.SetX(f0.X() - dxFront1);
				b0 = trk1->Nodes()[0]->Point3D(); b0.SetX(b0.X() - dxFront1);
				f1 = trk2->Nodes()[0]->Point3D(); f1.SetX(f1.X() - dxFront2);
				b1 = trk2->Nodes()[1]->Point3D(); b1.SetX(b1.X() - dxFront2);
				if (areCoLinear(c, f0, b0, f1, b1, distProjThr) && (c > cmax))
				{
					cmax = c; reverse = false; flip1 = true; flip2 = false;
					//best_idx = t;
					best_trk2 = trk2;
					dx1 = dxFront1; dx2 = dxFront2;
				}
			}
			else if ((cryoFront1 == cryoBack2) && (dxFront1 * dxBack2 < 0.0) &&
			    (fabs(dxFront1 + dxBack2) < xApaDistDiffThr))
			{
				if (fabs(dxFront1) < fabs(dxBack2)) dxBack2 = -dxFront1;
				else dxFront1 = -dxBack2;

				f0 = trk1->Nodes()[1]->Point3D(); f0.SetX(f0.X() - dxFront1);
				b0 = trk1->Nodes()[0]->Point3D(); b0.SetX(b0.X() - dxFront1);
				f1 = trk2->Nodes()[trk2->Nodes().size() - 1]->Point3D(); f1.SetX(f1.X() - dxBack2);
				b1 = trk2->Nodes()[trk2->Nodes().size() - 2]->Point3D(); b1.SetX(b1.X() - dxBack2);
				if (areCoLinear(c, f0, b0, f1, b1, distProjThr) && (c > cmax))
				{
					cmax = c; reverse = true; flip1 = false; flip2 = false;
					//best_idx = t;
					best_trk2 = trk2;
					dx1 = dxFront1; dx2 = dxBack2;
				}
			}
			else if ((cryoBack1 == cryoFront2) && (dxBack1 * dxFront2 < 0.0) &&
			    (fabs(dxBack1 + dxFront2) < xApaDistDiffThr))
			{
				if (fabs(dxBack1) < fabs(dxFront2)) dxFront2 = -dxBack1;
				else dxBack1 = -dxFront2;

				f0 = trk1->Nodes()[trk1->Nodes().size() - 2]->Point3D(); f0.SetX(f0.X() - dxBack1);
				b0 = trk1->Nodes()[trk1->Nodes().size() - 1]->Point3D(); b0.SetX(b0.X() - dxBack1);
				f1 = trk2->Nodes()[0]->Point3D(); f1.SetX(f1.X() - dxFront2);
				b1 = trk2->Nodes()[1]->Point3D(); b1.SetX(b1.X() - dxFront2);
				if (areCoLinear(c, f0, b0, f1, b1, distProjThr) && (c > cmax))
				{
					cmax = c; reverse = false; flip1 = false; flip2 = false;
					//best_idx = t;
					best_trk2 = trk2;
					dx1 = dxBack1; dx2 = dxFront2;
				}
			}
			else if ((cryoBack1 == cryoBack2) && (dxBack1 * dxBack2 < 0.0) &&
			    (fabs(dxBack1 + dxBack2) < xApaDistDiffThr))
			{
				if (fabs(dxBack1) < fabs(dxBack2)) dxBack2 = -dxBack1;
				else dxBack1 = -dxBack2;

				f0 = trk1->Nodes()[trk1->Nodes().size() - 2]->Point3D(); f0.SetX(f0.X() - dxBack1);
				b0 = trk1->Nodes()[trk1->Nodes().size() - 1]->Point3D(); b0.SetX(b0.X() - dxBack1);
				f1 = trk2->Nodes()[trk2->Nodes().size() - 1]->Point3D(); f1.SetX(f1.X() - dxBack2);
				b1 = trk2->Nodes()[trk2->Nodes().size() - 2]->Point3D(); b1.SetX(b1.X() - dxBack2);
				if (areCoLinear(c, f0, b0, f1, b1, distProjThr) && (c > cmax))
				{
					cmax = c; reverse = false; flip1 = false; flip2 = true;
					//best_idx = t;
					best_trk2 = trk2;
					dx1 = dxBack1; dx2 = dxBack2;
				}
			}
		}
		if (best_trk2)
		{
			if (flip1) trk1->Flip();
			if (flip2) best_trk2->Flip();

			trk1->GetRoot()->ApplyXShiftInTree(-dx1);
			best_trk2->GetRoot()->ApplyXShiftInTree(-dx2);

			if (reverse)
			{
				best_trk2->SetSubsequentTrack(trk1);
				trk1->SetPrecedingTrack(best_trk2);
			}
			else
			{
				trk1->SetSubsequentTrack(best_trk2);
				best_trk2->SetPrecedingTrack(trk1);
			}
		}
	}
}
// ------------------------------------------------------

bool PMAlgTrackMaker::reassignHits(const std::vector< art::Ptr<recob::Hit> > & hits,
	const std::vector< art::Ptr<recob::Cluster> > & clusters,
	pma::TrkCandidateColl & tracks, size_t trk_idx, double dist2)
{
	pma::Track3D* trk1 = tracks[trk_idx].Track();

	bool result = false;
	if ((hits.size() > 1) || (dist2 > 1.0)) // min. 2 hits or single hit separated from the rest
	{
		pma::Track3D* best_trk = 0;

		size_t best_u = 0, n_max = 0;
		for (size_t u = 0; u < tracks.size(); u++)
			if (trk_idx != u)
		{
			pma::Track3D* trk2 = tracks[u].Track();
			size_t n = fProjectionMatchingAlg.testHits(*trk2, hits);
			if (n > n_max) { n_max = n; best_u = u; best_trk = trk2; }
		}

		if (best_trk && (n_max >= hits.size() / 3)) // /2
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "  Reassign(v1) " << n_max << " hits." << std::endl;

			trk1->RemoveHits(hits);
			trk1->CleanupTails();
			trk1->ShiftEndsToHits();

			pma::Track3D* ext = fProjectionMatchingAlg.extendTrack(*best_trk, hits,	false);
			ext->SortHits(); ext->ShiftEndsToHits();
			if (fProjectionMatchingAlg.isContained(*ext))
			{
				tracks[best_u].SetTrack(ext); // and this deletes best_trk stored at best_u
				result = true;
			}
			else delete ext;
		}
		else if (!clusters.empty() && (hits.size() >= fMinSeedSize2ndPass))
		{
			size_t minSizeCompl = hits.size() / 8;  // much smaller minimum required in complementary views
			if (minSizeCompl < 3) minSizeCompl = 3; // but at least three hits!

			geo::View_t first_view = (geo::View_t)hits.front()->WireID().Plane;
			unsigned int tpc = hits.front()->WireID().TPC;
			unsigned int cryo = hits.front()->WireID().Cryostat;

			pma::TrkCandidate candidate = matchCluster(
				hits, clusters, minSizeCompl, tpc, cryo, first_view, -1, -1);

			if (candidate.IsGood())
			{
				mf::LogVerbatim("PMAlgTrackMaker") << "  Add new track, cut hits from source track." << std::endl;
				tracks.push_back(candidate);

				trk1->RemoveHits(hits);
				trk1->CleanupTails();
				trk1->ShiftEndsToHits();
			}
		}
	}
	else if ((hits.size() == 1) || (dist2 > 2.25)) // dist > 1.5cm
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "  Cut single-view isolated hit." << std::endl;
		trk1->RemoveHits(hits);
		trk1->CleanupTails();
		trk1->ShiftEndsToHits();
	}
	return result;
}

bool PMAlgTrackMaker::reassignHits(const std::vector< art::Ptr<recob::Hit> > & hits,
	pma::TrkCandidateColl & tracks, size_t trk_idx, double dist2)
{
	pma::Track3D* trk1 = tracks[trk_idx].Track();

	bool result = false;
	if ((hits.size() > 1) || (dist2 > 1.0)) // min. 2 hits or single hit separated from the rest
	{
		pma::Track3D* best_trk = 0;

		size_t n_max = 0;
		for (size_t u = 0; u < tracks.size(); u++)
			if (trk_idx != u)
		{
			pma::Track3D* trk2 = tracks[u].Track();
			size_t n = fProjectionMatchingAlg.testHits(*trk2, hits, 0.5);
			if (n > n_max) { n_max = n; best_trk = trk2; }
		}

		if (best_trk && (n_max >= (size_t)(0.8 * hits.size()))) // most hits!
		{
			mf::LogVerbatim("PMAlgTrackMaker") << "  Reassign(v2) " << n_max << " hits." << std::endl;

			trk1->RemoveHits(hits);
			trk1->CleanupTails();
			trk1->ShiftEndsToHits();

			best_trk->AddHits(hits);

			result = true;
		}
	}
	else if ((hits.size() == 1) || (dist2 > 2.25)) // dist > 1.5cm
	{
		mf::LogVerbatim("PMAlgTrackMaker") << "  Cut single-view isolated hit." << std::endl;
		trk1->RemoveHits(hits);
		trk1->CleanupTails();
		trk1->ShiftEndsToHits();

		result = true;
	}

	if (result)
	{
		// reopt trees
	}

	return result;
}

double PMAlgTrackMaker::collectSingleViewEnd(pma::Track3D & trk,
	std::vector< art::Ptr<recob::Hit> > & hits)
{
	size_t idx = 0;
	while ((idx < trk.size() - 1) && !(trk[idx]->IsEnabled()))
	{
		hits.push_back(trk[idx++]->Hit2DPtr());
	}

	double d2 = 0.0;
	if (idx > 0)
	{
		if ((idx < trk.size() - 1) &&
		    (trk[idx]->View2D() == trk[idx - 1]->View2D()))
		{
			double dprev = pma::Dist2(trk[idx]->Point3D(), trk[idx - 1]->Point3D());
			double dnext = pma::Dist2(trk[idx]->Point3D(), trk[idx + 1]->Point3D());
			if (dprev < dnext)
			{
				hits.push_back(trk[idx++]->Hit2DPtr());
			}
		}
		d2 = pma::Dist2(trk[idx]->Point3D(), trk[idx - 1]->Point3D());
	}
	return d2;
}

double PMAlgTrackMaker::collectSingleViewFront(pma::Track3D & trk,
	std::vector< art::Ptr<recob::Hit> > & hits)
{
	size_t idx = trk.size() - 1;
	while ((idx > 0) && !(trk[idx]->IsEnabled()))
	{
		hits.push_back(trk[idx--]->Hit2DPtr());
	}

	double d2 = 0.0;
	if (idx < trk.size() - 1)
	{
		if ((idx > 0) &&
		    (trk[idx]->View2D() == trk[idx + 1]->View2D()))
		{
			double dprev = pma::Dist2(trk[idx]->Point3D(), trk[idx + 1]->Point3D());
			double dnext = pma::Dist2(trk[idx]->Point3D(), trk[idx - 1]->Point3D());
			if (dprev < dnext)
			{
				hits.push_back(trk[idx--]->Hit2DPtr());
			}
		}
		d2 = pma::Dist2(trk[idx]->Point3D(), trk[idx + 1]->Point3D());
	}
	return d2;
}

bool PMAlgTrackMaker::reassignSingleViewEnds(pma::TrkCandidateColl & tracks,
	const std::vector< art::Ptr<recob::Cluster> > & clusters)
{
	bool result = false;
	for (size_t t = 0; t < tracks.size(); t++)
	{
		pma::Track3D & trk = *(tracks[t].Track());
		if (trk.size() < 6) continue;

		trk.DisableSingleViewEnds();

		std::vector< art::Ptr<recob::Hit> > hits;

		double d2 = collectSingleViewEnd(trk, hits);
		result |= reassignHits(hits, clusters, tracks, t, d2);

		hits.clear();

		d2 = collectSingleViewFront(trk, hits);
		result |= reassignHits(hits, clusters, tracks, t, d2);

		trk.SelectHits();
	}
	return result;
}

bool PMAlgTrackMaker::reassignSingleViewEnds(pma::TrkCandidateColl & tracks)
{
	bool result = false;
	for (size_t t = 0; t < tracks.size(); t++)
	{
		pma::Track3D & trk = *(tracks[t].Track());
		if (trk.size() < 6) continue;

		trk.DisableSingleViewEnds();

		std::vector< art::Ptr<recob::Hit> > hits;

		double d2 = collectSingleViewEnd(trk, hits);
		result |= reassignHits(hits, tracks, t, d2);

		hits.clear();

		d2 = collectSingleViewFront(trk, hits);
		result |= reassignHits(hits, tracks, t, d2);

		trk.SelectHits();
	}
	return result;
}
// ------------------------------------------------------

void PMAlgTrackMaker::guideEndpoints(pma::TrkCandidateColl& tracks)
{
	for (auto const & t : tracks.tracks())
	{
		auto & trk = *(t.Track());
		fProjectionMatchingAlg.guideEndpoints(trk, fHitMap[trk.FrontCryo()][trk.FrontTPC()]);
	}
}
// ------------------------------------------------------

bool PMAlgTrackMaker::sortHits(const art::Event& evt)
{
	fHitMap.clear(); fCluHits.clear();

	art::Handle< std::vector<recob::Hit> > allHitListHandle;
	art::Handle< std::vector<recob::Cluster> > cluListHandle, splitCluHandle;
	std::vector< art::Ptr<recob::Hit> > allhitlist;
	if (evt.getByLabel(fHitModuleLabel, splitCluHandle) &&   // clusters that tag em-like hits
	    evt.getByLabel(fCluModuleLabel, allHitListHandle) && // all hits associated to both cluster sets
	    evt.getByLabel(fCluModuleLabel, cluListHandle))      // clusters used to build 3D tracks
	{
		art::fill_ptr_vector(allhitlist, allHitListHandle);

		mf::LogVerbatim("PMAlgTrackMaker") << "Sort all hits for validation...";
		unsigned int cryo, tpc, view;
		for (auto const& h : allhitlist) // all hits used for validation
		{
			cryo = h->WireID().Cryostat;
			tpc = h->WireID().TPC;
			view = h->WireID().Plane;

			fHitMap[cryo][tpc][view].push_back(h);
		}
		mf::LogVerbatim("PMAlgTrackMaker") << "...done.";

		mf::LogVerbatim("PMAlgTrackMaker") << "Filter track-like clusters...";
		fCluHits.reserve(cluListHandle->size());
		art::FindManyP< recob::Hit > fbp(cluListHandle, evt, fCluModuleLabel);
		art::FindManyP< recob::Hit > fem(splitCluHandle, evt, fHitModuleLabel);
		for (size_t i = 0; i < cluListHandle->size(); ++i)
		{
			auto v = fbp.at(i);

			fCluHits.push_back(std::vector< art::Ptr<recob::Hit> >());

			for (auto const & h : v)
			{
				bool trkLike = true;
				if (fCluModuleLabel != fHitModuleLabel)
				{
					for (size_t j = 0; j < splitCluHandle->size(); ++j)
					{
						auto u = fem.at(j);
						for (auto const & g : u) // is hit clustered in one of em-like?
						{
							if (g.key() == h.key())
							{
								trkLike = false; break;
							}
						}
					}
				}
				if (trkLike) fCluHits.back().push_back(h);
			}
		}
		if (fCluHits.size() != cluListHandle->size())
		{
			mf::LogError("PMAlgTrackMaker") << "Hit-cluster map incorrect, better skip this event.";
			return false;
		}

		mf::LogVerbatim("PMAlgTrackMaker") << "...done, " << fCluHits.size() << " clusters for 3D tracking.";
		return true;
	}
	else return false;
}
// ------------------------------------------------------

bool PMAlgTrackMaker::sortHitsPfp(const art::Event& evt)
{
	fHitMap.clear(); fCluHits.clear(); fPfpClusters.clear(); fPfpPdgCodes.clear();
	fPfpVtx.clear();

	art::Handle< std::vector<recob::Hit> > allHitListHandle;
	art::Handle< std::vector<recob::Cluster> > cluListHandle;
	art::Handle< std::vector<recob::PFParticle> > pfparticleHandle;
	std::vector< art::Ptr<recob::Hit> > allhitlist;
	if (evt.getByLabel(fHitModuleLabel, allHitListHandle) && // all hits used to make clusters and PFParticles
	    evt.getByLabel(fCluModuleLabel, cluListHandle) &&    // clusters associated to PFParticles
	    evt.getByLabel(fCluModuleLabel, pfparticleHandle))   // and finally PFParticles
	{
		art::fill_ptr_vector(allhitlist, allHitListHandle);

		mf::LogVerbatim("PMAlgTrackMaker") << "Sort all hits for validation...";
		unsigned int cryo, tpc, view;
		for (auto const& h : allhitlist) // all hits used for validation
		{
			cryo = h->WireID().Cryostat;
			tpc = h->WireID().TPC;
			view = h->WireID().Plane;

			fHitMap[cryo][tpc][view].push_back(h);
		}
		mf::LogVerbatim("PMAlgTrackMaker") << "...done, " << allhitlist.size() << "hits.";

		mf::LogVerbatim("PMAlgTrackMaker") << "Sort hits by clusters assigned to PFParticles...";
		fCluHits.reserve(cluListHandle->size());
		art::FindManyP< recob::Hit > fbp(cluListHandle, evt, fCluModuleLabel);
		art::FindManyP< recob::Cluster > fpf(pfparticleHandle, evt, fCluModuleLabel);	
		art::FindManyP< recob::Vertex > fvf(pfparticleHandle, evt, fCluModuleLabel);

		for (size_t i = 0; i < cluListHandle->size(); ++i)
		{
			fCluHits.push_back(std::vector< art::Ptr<recob::Hit> >());
		}
		for (size_t i = 0; i < pfparticleHandle->size(); ++i)
		{
			fPfpPdgCodes[i] = pfparticleHandle->at(i).PdgCode();

			auto cv = fpf.at(i);
			for (const auto & c : cv)
			{
				fPfpClusters[i].push_back(c);

				if (fCluHits[c.key()].empty())
				{
					auto hv = fbp.at(c.key());
					fCluHits[c.key()].reserve(hv.size());
					for (auto const & h : hv) fCluHits[c.key()].push_back(h);
				}
			}

			if (fvf.at(i).size())
				fPfpVtx[i] = fvf.at(i).front();
		}

		mf::LogVerbatim("PMAlgTrackMaker") << "...done, "
			<< fCluHits.size() << " clusters from "
			<< fPfpClusters.size() << " pfparticles for 3D tracking.";
		return true;
	}
	else return false;
}
// ------------------------------------------------------
#include "PMAlgTrackMaker_module_continued.h"

DEFINE_ART_MODULE(PMAlgTrackMaker)

} // namespace trkf

