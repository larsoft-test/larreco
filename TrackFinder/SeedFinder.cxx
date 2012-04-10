//
// Name: SeedFinder.cxx
//
// Purpose: Implementation file for module SeedFinder.
//
// Ben Jones, MIT
//

#include <vector>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TrackFinder/SeedFinder.h"
#include "TrackFinder/SeedFinderService.h"
#include "TrackFinder/SpacePointService.h"
#include "Geometry/geo.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "RecoBase/Hit.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/Prong.h"
#include "Utilities/AssociationUtil.h"

namespace trkf {

  SeedFinder::SeedFinder(const fhicl::ParameterSet& pset) :
    fFilter(true),
    fMerge(false)
  {
    reconfigure(pset);
    produces<std::vector<recob::Seed> >();
  
  
    mf::LogInfo("SeedFinder") 
      << "SeedFinder configured with the following parameters:\n"
      << "  ClusterModuleLabel = " << fClusterModuleLabel << "\n"
      << "  Filter = " << fFilter << "\n"
      << "  Merge = " << fMerge;
  }

  SeedFinder::~SeedFinder()
  {
  }

  void SeedFinder::reconfigure(fhicl::ParameterSet const& pset)
  {
    fClusterModuleLabel = pset.get<std::string>("ClusterModuleLabel");
    fFilter             = pset.get<bool>("Filter");
    fMerge              = pset.get<bool>("Merge");
    fSeedMode           = pset.get<int>("SeedMode");
  }

  void SeedFinder::beginJob()
  {}

  void SeedFinder::produce(art::Event& evt)
  {
    art::ServiceHandle<trkf::SeedFinderService> sfsvc;  
 
    int EventNumber = evt.id().event();
    sfsvc->SetEventID(EventNumber);

    art::Handle< std::vector<recob::Cluster> > clusterh;
    evt.getByLabel(fClusterModuleLabel, clusterh);

    std::auto_ptr<std::vector<recob::Seed> > seeds(new std::vector<recob::Seed>);

    std::vector<std::vector<recob::SpacePoint> > SpacePointVectors = GetSpacePoints(clusterh);
    if(SpacePointVectors.size() > 0)
      {
	for(std::vector<std::vector<recob::SpacePoint> >::const_iterator it=SpacePointVectors.begin(); 
	    it!=SpacePointVectors.end();
	    it++)
	  {
	    std::vector<recob::Seed*> SeedFinderOutput;
	    if(fSeedMode==0)
	      SeedFinderOutput = sfsvc->FindSeedExhaustively(*it);
	    else if(fSeedMode==1)
	      SeedFinderOutput = sfsvc->FindAsManySeedsAsPossible(*it);
	    else 
	      throw cet::exception("SeedFinder") << 
		"Unkown seed mode " << fSeedMode<<"\n";

	    if(SeedFinderOutput.size()>0)
	      for(int i=0; i!=SeedFinderOutput.size(); i++)
		if(SeedFinderOutput.at(i)->IsValid())
		  seeds->push_back(*(SeedFinderOutput.at(i)));
	    
	  }			  
      }
    else
      std::cout<<"Seed finder made no seeds : no space points in event"<<std::endl;
    
    evt.put(seeds);
    
  }

  std::vector<std::vector<recob::SpacePoint> > SeedFinder::GetSpacePoints(art::Handle<std::vector<recob::Cluster> > clusterh)
  {
    // Get Services.

    art::ServiceHandle<trkf::SpacePointService> sptsvc;
    art::ServiceHandle<geo::Geometry> geom;

    // Get clusters.

    std::vector<std::vector<recob::SpacePoint> > SpacePointVectors;

    // Make a double or triple loop over clusters in distinct views
    // (depending on minimum number of views configured in SpacePointService).

    if(clusterh.isValid()) {



      art::PtrVector<recob::Hit> hits;      

      // Loop over first cluster.

      int nclus = clusterh->size();
      for(int iclus = 0; iclus < nclus; ++iclus) {
	art::Ptr<recob::Cluster> piclus(clusterh, iclus);
	geo::View_t iview = piclus->View();

	// Test first view.

	if((iview == geo::kU && sptsvc->enableU()) ||
	   (iview == geo::kV && sptsvc->enableV()) ||
	   (iview == geo::kW && sptsvc->enableW())) {

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

	    if(((jview == geo::kU && sptsvc->enableU()) ||
		(jview == geo::kV && sptsvc->enableV()) ||
		(jview == geo::kW && sptsvc->enableW()))
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
	  

	      // Loop over third cluster.

	      for(int kclus = 0; kclus < jclus; ++kclus) {
		art::Ptr<recob::Cluster> pkclus(clusterh, kclus);
		geo::View_t kview = pkclus->View();

		// Test third view.

		if(((kview == geo::kU && sptsvc->enableU()) ||
		    (kview == geo::kV && sptsvc->enableV()) ||
		    (kview == geo::kW && sptsvc->enableW()))
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

		  std::vector<recob::SpacePoint> spts;
		  sptsvc->makeSpacePoints(hits, spts,
					  fFilter, fMerge, 0., 0.);

		  if(spts.size() > 0) {
		    SpacePointVectors.push_back(spts);
		  }
		}
	      }
	    }
	  }
	}
      }

    }
    return SpacePointVectors;


  }

  void SeedFinder::endJob()
  {

  }
}
