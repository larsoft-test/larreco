//
// Name: SpacePointCheater.cxx
//
// Purpose: Implementation file for module SpacePointCheater.
//
// Created: 15-Dec-2011  H. Greenlee
//

#include <vector>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TrackFinder/SpacePointCheater.h"
#include "Geometry/geo.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Prong.h"
#include "Utilities/AssociationUtil.h"

namespace trkf {

  SpacePointCheater::SpacePointCheater(const fhicl::ParameterSet& pset) :
    //
    // Purpose: Constructor.
    //
    // Arguments: pset - Module parameters.
    //
    fSptalg(pset.get<fhicl::ParameterSet>("SpacePointAlg")),
    fMinHits(0),
    fFilter(true),
    fMerge(false),
    fNumEvent(0),
    fNumProng2(0),
    fNumProng3(0),
    fNumSpt2(0),
    fNumSpt3(0)
  {
    reconfigure(pset);
    produces<std::vector<recob::Prong>                 >();
    produces<art::Assns<recob::Prong, recob::Cluster>  >();
    produces<art::Assns<recob::Prong, recob::Hit>      >();
    produces<std::vector<recob::SpacePoint>            >();
    produces<art::Assns<recob::SpacePoint, recob::Hit> >();

    // Report.

    mf::LogInfo("SpacePointCheater") 
      << "SpacePointCheater configured with the following parameters:\n"
      << "  ClusterModuleLabel = " << fClusterModuleLabel << "\n"
      << "  G4ModuleLabel = " << fG4ModuleLabel << "\n"
      << "  Minimum Hits per Cluster = " << fMinHits << "\n"
      << "  Filter = " << fFilter << "\n"
      << "  Merge = " << fMerge;
  }

  SpacePointCheater::~SpacePointCheater()
  //
  // Purpose: Destructor.
  //
  {}

  void SpacePointCheater::reconfigure(fhicl::ParameterSet const& pset)
  //
  // Purpose: Reconfigure method.
  //
  // Arguments: pset - Configuration parameters.
  //
  {
    fSptalg.reconfigure(pset.get<fhicl::ParameterSet>("SpacePointAlg"));
    fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
    fG4ModuleLabel = pset.get<std::string>("G4ModuleLabel");
    fMinHits = pset.get<unsigned int>("MinHits");
    fFilter = pset.get<bool>("Filter");
    fMerge = pset.get<bool>("Merge");
  }

  void SpacePointCheater::beginJob()
  {}

  void SpacePointCheater::produce(art::Event& evt)
  //
  // Purpose: Produce method.
  //
  // Arguments: event - Art event.
  //
  {
    ++fNumEvent;

    // Get Services.

    art::ServiceHandle<geo::Geometry> geom;

    // Get SimChannels.
    // Make a vector where each channel in the detector is an entry
    std::vector<const sim::SimChannel*> simchanh;
    std::vector<const sim::SimChannel*> simchanv(geom->Nchannels(),0);
    evt.getView(fG4ModuleLabel, simchanh);
    for(size_t i = 0; i < simchanh.size(); ++i)
      simchanv[simchanh[i]->Channel()] = simchanh[i];

    // Get clusters.

    art::Handle< std::vector<recob::Cluster> > clusterh;
    evt.getByLabel(fClusterModuleLabel, clusterh);

    // Make a double or triple loop over clusters in distinct views
    // (depending on minimum number of views configured in SpacePointAlg).

    if(clusterh.isValid()) {

      // Make a collection of prongs that will eventually be inserted into the event.

      std::auto_ptr<std::vector<recob::Prong> > prongs(new std::vector<recob::Prong>);
      std::auto_ptr< art::Assns<recob::Prong, recob::Cluster> > assn(new art::Assns<recob::Prong, recob::Cluster>);
      std::auto_ptr< art::Assns<recob::Prong, recob::Hit> > hassn(new art::Assns<recob::Prong, recob::Hit>);

      // Make a collection of space points that will be inserted into the event.

      std::auto_ptr<std::vector<recob::SpacePoint> > spts(new std::vector<recob::SpacePoint>);
      std::auto_ptr< art::Assns<recob::SpacePoint, recob::Hit> > spassn(new art::Assns<recob::SpacePoint, recob::Hit>);

      // Make a hit vector which will be used to store hits to be passed
      // to SpacePointAlg.

      art::PtrVector<recob::Hit> hits;      

      // Loop over first cluster.

      int nclus = clusterh->size();
      for(int iclus = 0; iclus < nclus; ++iclus) {
	art::Ptr<recob::Cluster> piclus(clusterh, iclus);
	geo::View_t iview = piclus->View();

	// Test first view.

	if(piclus->Hits().size() >= fMinHits &&
	   ((iview == geo::kU && fSptalg.enableU()) ||
	    (iview == geo::kV && fSptalg.enableV()) ||
	    (iview == geo::kW && fSptalg.enableW()))) {

	  // Store hits from first view into hit vector.

	  art::PtrVector<recob::Hit> ihits = piclus->Hits();
	  unsigned int nihits = ihits.size();
	  hits.clear();
	  hits.reserve(nihits);
	  for(art::PtrVector<recob::Hit>::const_iterator i = ihits.begin();
	      i != ihits.end(); ++i)
	    hits.push_back(*i);
	  
	  // Loop over second cluster.

	  for(int jclus = 0; jclus < iclus; ++jclus) {
	    art::Ptr<recob::Cluster> pjclus(clusterh, jclus);
	    geo::View_t jview = pjclus->View();

	    // Test second view.

	    if(pjclus->Hits().size() >= fMinHits &&
	       ((jview == geo::kU && fSptalg.enableU()) ||
		(jview == geo::kV && fSptalg.enableV()) ||
		(jview == geo::kW && fSptalg.enableW()))
	       && jview != iview) {

	      // Store hits from second view into hit vector.

	      art::PtrVector<recob::Hit> jhits = pjclus->Hits();
	      unsigned int njhits = jhits.size();
	      assert(hits.size() >= nihits);
	      //hits.resize(nihits);
	      while(hits.size() > nihits)
		hits.pop_back();
	      assert(hits.size() == nihits);
	      hits.reserve(nihits + njhits);
	      for(art::PtrVector<recob::Hit>::const_iterator j = jhits.begin();
		  j != jhits.end(); ++j)
		hits.push_back(*j);
	  
	      // If two-view space points are allowed, make them here.

	      if(fSptalg.minViews() <= 2) {
		std::vector<recob::SpacePoint> new_spts;
		fSptalg.makeMCTruthSpacePoints(hits, new_spts, simchanv,
 					       fFilter, fMerge, 0., 0.);

		// If we found some space points, make a prong.

		if(new_spts.size() > 0) {
		  fNumSpt2 += new_spts.size();
		  art::PtrVector<recob::Cluster> clusters;
		  clusters.reserve(2);
		  clusters.push_back(piclus);
		  clusters.push_back(pjclus);
		  prongs->push_back(recob::Prong(clusters, new_spts));
		  util::CreateAssn(*this, evt, *(prongs.get()), clusters,*(assn.get()));

		  // associate the cluster hits with this prong as well
		  for(size_t c = 0; c < clusters.size(); ++c){
		    art::PtrVector<recob::Hit> hits = util::FindManyP<recob::Hit>(clusters, evt, 
										  fClusterModuleLabel, c);
		    util::CreateAssn(*this, evt, *(prongs.get()), hits, *(hassn.get()));
		  }

		  // Insert newly found space points into event collection.

		  int nspt = spts->size();
		  spts->insert(spts->end(), new_spts.begin(), new_spts.end());

		  // Associate space points with hits.

		  for(unsigned int ispt = nspt; ispt < spts->size(); ++ispt) {
		    const recob::SpacePoint& spt = (*spts)[ispt];
		    const art::PtrVector<recob::Hit>& hits = fSptalg.getAssociatedHits(spt);
		    util::CreateAssn(*this, evt, *spts, hits, *spassn, ispt);
		  }

		  ++fNumProng2;
		}
	      }

	      // Loop over third cluster.

	      for(int kclus = 0; kclus < jclus; ++kclus) {
		art::Ptr<recob::Cluster> pkclus(clusterh, kclus);
		geo::View_t kview = pkclus->View();

		// Test third view.

		if(pkclus->Hits().size() >= fMinHits &&
		   ((kview == geo::kU && fSptalg.enableU()) ||
		    (kview == geo::kV && fSptalg.enableV()) ||
		    (kview == geo::kW && fSptalg.enableW()))
		   && kview != iview && kview != jview) {

		  // Store hits from third view into hit vector.

		  art::PtrVector<recob::Hit> khits = pkclus->Hits();
		  unsigned int nkhits = khits.size();
		  assert(hits.size() >= nihits + njhits);
		  //hits.resize(nihits + njhits);
		  while(hits.size() > nihits + njhits)
		    hits.pop_back();
		  assert(hits.size() == nihits + njhits);
		  hits.reserve(nihits + njhits + nkhits);
		  for(art::PtrVector<recob::Hit>::const_iterator k = khits.begin();
		      k != khits.end(); ++k)
		    hits.push_back(*k);

		  // Make three-view space points.

		  std::vector<recob::SpacePoint> new_spts;
		  fSptalg.makeMCTruthSpacePoints(hits, new_spts, simchanv,
						 fFilter, fMerge, 0., 0.);

		  // If we found some space points, make a prong.

		  if(new_spts.size() > 0) {
		    fNumSpt3 += new_spts.size();
		    art::PtrVector<recob::Cluster> clusters;
		    clusters.reserve(3);
		    clusters.push_back(piclus);
		    clusters.push_back(pjclus);
		    clusters.push_back(pkclus);
		    prongs->push_back(recob::Prong(clusters, new_spts));
		    util::CreateAssn(*this, evt, *(prongs.get()), clusters,*(assn.get()));
		    
		    // associate the cluster hits with this prong as well
		    for(size_t c = 0; c < clusters.size(); ++c){
		      art::PtrVector<recob::Hit> hits = util::FindManyP<recob::Hit>(clusters, evt, 
										    fClusterModuleLabel, c);
		      util::CreateAssn(*this, evt, *(prongs.get()), hits, *(hassn.get()));
		    }

		    // Insert newly found space points into event collection.

		    int nspt = spts->size();
		    spts->insert(spts->end(), new_spts.begin(), new_spts.end());

		    // Associate space points with hits.

		    for(unsigned int ispt = nspt; ispt < spts->size(); ++ispt) {
		      const recob::SpacePoint& spt = (*spts)[ispt];
		      const art::PtrVector<recob::Hit>& hits = fSptalg.getAssociatedHits(spt);
		      util::CreateAssn(*this, evt, *spts, hits, *spassn, ispt);
		    }

		    ++fNumProng3;
		  }
		}
	      }
	    }
	  }
	}
      }

      // Add prongs and associations to event.

      evt.put(prongs);
      evt.put(assn);
      evt.put(hassn);
      evt.put(spts);
      evt.put(spassn);
    }
  }

  void SpacePointCheater::endJob()
  //
  // Purpose: Print summary.
  //
  {
    mf::LogInfo("SpacePointCheater") 
      << "SpacePointCheater statistics:\n"
      << "  Number of events = " << fNumEvent << "\n"
      << "  Number of 2-view prongs created = " << fNumProng2 << "\n"
      << "  Number of 3-view prongs created = " << fNumProng3 << "\n"
      << "  Number of 2-view space points created = " << fNumSpt2 << "\n"
      << "  Number of 3-view space points created = " << fNumSpt3;
  }
}
