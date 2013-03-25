#ifndef SEEDFINDER_H
#define SEEDFINDER_H

//
// Name: SeedFinderModule.h
//
//
// Ben Jones, MIT, April 2012
//   bjpjones@mit.edu
//

#include "art/Framework/Core/EDProducer.h"
#include "TrackFinder/SpacePointAlg.h"
#include "TrackFinder/SeedFinderAlgorithm.h"

namespace recob
{
  class SpacePoint;
  class Cluster;
  class Seed;
  class Hit;
}

namespace trkf {

  class SeedFinderModule : public art::EDProducer
  {
  public:
 
    // Constructors, destructor

    explicit SeedFinderModule(fhicl::ParameterSet const& pset);
    virtual ~SeedFinderModule();

    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void produce(art::Event& evt);
    
    void endJob();
  
    
    std::vector<std::vector<recob::SpacePoint> > GetSpacePointsFromClusters(std::string fInputModuleLabel, art::Event& evt);

    art::PtrVector<recob::Hit>       GetHitsFromEvent(std::string HitModuleLabel, art::Event & evt);



  private:

   
    // Fcl Attributes.

    SeedFinderAlgorithm    fSeedAlg;                  // Algorithm object
    SpacePointAlg          fSptalg;
    std::string            fInputModuleLabel;         // Where to find hits, if we need them
    int                    fInputSource;              // 1: Use Clusters
                                                      // 2: Use Hits

  };
  
}

#endif // SEEDFINDER_H



#include "art/Framework/Core/ModuleMacros.h" 


namespace trkf {
  DEFINE_ART_MODULE(SeedFinderModule);
}

//
// Name: SeedFinderModule.cxx
//
// Purpose: Implementation file for module SeedFinderModule.
//
// Ben Jones, MIT
//

#include <vector>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Geometry/Geometry.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "TMatrixD.h"
#include "TVectorD.h"

namespace trkf {

  //----------------------------------------------------------------------------
  SeedFinderModule::SeedFinderModule(const fhicl::ParameterSet& pset) :
    fSeedAlg(pset.get<fhicl::ParameterSet>("SeedAlg")),
    fSptalg(pset.get<fhicl::ParameterSet>("SpacePointAlg"))
  {
    reconfigure(pset);
    produces<std::vector<recob::Seed> >();
  
  
  }

  //----------------------------------------------------------------------------
  SeedFinderModule::~SeedFinderModule()
  {
  }

  //----------------------------------------------------------------------------
  void SeedFinderModule::reconfigure(fhicl::ParameterSet const& pset)
  {
    fSeedAlg               = SeedFinderAlgorithm ( pset.get<fhicl::ParameterSet>("SeedAlg") );
    fSptalg                = SpacePointAlg       ( pset.get<fhicl::ParameterSet>("SpacePointAlg") );
    fInputModuleLabel      = pset.get<std::string>("InputModuleLabel");
    fInputSource           = pset.get<int>("InputSource");
  
  }

  //----------------------------------------------------------------------------
  void SeedFinderModule::beginJob()
  {}

  //----------------------------------------------------------------------------
  void SeedFinderModule::produce(art::Event& evt)
  {
    
    std::unique_ptr<std::vector<recob::Seed> > seeds(new std::vector<recob::Seed>);

    std::vector<std::vector<recob::SpacePoint> > SpacePointsWithSeeds;
    std::vector<recob::Seed> SeedVector;
    
    if(fInputSource==1)
      {
	std::cout<<"SeedFinder: Getting space points from clusters"<<std::endl;
	
	std::vector<std::vector<recob::SpacePoint> > SpacePointVectors;
	SpacePointVectors = GetSpacePointsFromClusters(fInputModuleLabel, evt);
	
	std::vector<recob::Seed> SeedsToAdd;
	std::vector<std::vector<recob::SpacePoint> > ReturnedSPs;

	for(size_t i=0; i!=SpacePointVectors.size(); ++i)
	  {
	    SeedsToAdd = fSeedAlg.FindSeeds(SpacePointVectors.at(i),ReturnedSPs);
	    for(size_t j=0; j!=SeedsToAdd.size(); j++)
	      {
		SeedVector.push_back(SeedsToAdd.at(j));
		SpacePointsWithSeeds.push_back(ReturnedSPs.at(j));
	      }
	  }      

      }else if(fInputSource==0)
      {
	
	std::cout<<"SeedFinder: Getting space points from hits"<<std::endl;
	art::PtrVector<recob::Hit> Hits = GetHitsFromEvent(fInputModuleLabel, evt);
	std::cout<<"Hits extracted from event : " << Hits.size()<<std::endl;
       	std::vector<recob::SpacePoint> SPsFromHits = fSeedAlg.GetSpacePointsFromHitVector(Hits);
	std::cout<<"SPs extracted from hits : " << SPsFromHits.size()<<std::endl;
 	SeedVector = fSeedAlg.FindSeeds(SPsFromHits, SpacePointsWithSeeds);
      }
    else	  
      {
	throw cet::exception("SeedFinderModule") << 
	  "Unkown source mode " << fInputSource<<"\n";
      }

    if(SeedVector.size()>0)
      {
	for(size_t i=0; i!=SeedVector.size(); ++i)
	  {
	    seeds->push_back(SeedVector.at(i));
	  }
	
      }
    else
      std::cout<<"Seed finder made no seeds : no space points in event"<<std::endl;
    
    evt.put(std::move(seeds));
    


  }






  
  //----------------------------------------------------------------------------

  std::vector<std::vector<recob::SpacePoint> > SeedFinderModule::GetSpacePointsFromClusters(std::string ClusterModuleLabel, art::Event& evt)
  {
    // Get Services.

    art::ServiceHandle<geo::Geometry> geom;

    std::vector<art::Ptr<recob::Cluster> > Clusters;
   
    art::Handle< std::vector<recob::Cluster> > clusterh;
    evt.getByLabel(ClusterModuleLabel, clusterh);

    if(clusterh.isValid()) {
      art::fill_ptr_vector(Clusters, clusterh);
    }
    art::FindManyP<recob::Hit> fm(clusterh, evt, ClusterModuleLabel);


    std::vector<std::vector<recob::SpacePoint> > SpacePointVectors;

    // Make a double or triple loop over clusters in distinct views
    // (depending on minimum number of views configured in SpacePointAlg).

    art::PtrVector<recob::Hit> hits;      
    
    // Loop over first cluster.
    
    int nclus = Clusters.size();
    for(int iclus = 0; iclus < nclus; ++iclus) {
      art::Ptr<recob::Cluster> piclus = Clusters.at(iclus);
      geo::View_t iview = piclus->View();
      
      // Test first view.
      
      if((iview == geo::kU && fSptalg.enableU()) ||
	 (iview == geo::kV && fSptalg.enableV()) ||
	 (iview == geo::kW && fSptalg.enableW())) {
	
	// Store hits from first view into hit vector.
	
	std::vector< art::Ptr<recob::Hit> > ihits = fm.at(iclus);
	unsigned int nihits = ihits.size();
	hits.clear();
	hits.reserve(nihits);
	for(std::vector< art::Ptr<recob::Hit> >::const_iterator i = ihits.begin();
	    i != ihits.end(); ++i)
	  hits.push_back(*i);
	
	// Loop over second cluster.

	for(int jclus = 0; jclus < iclus; ++jclus) {
	  art::Ptr<recob::Cluster> pjclus = Clusters.at(jclus);
	  geo::View_t jview = pjclus->View();
	  
	    // Test second view.

	    if(((jview == geo::kU && fSptalg.enableU()) ||
		(jview == geo::kV && fSptalg.enableV()) ||
		(jview == geo::kW && fSptalg.enableW()))
	       && jview != iview) {

	      // Store hits from second view into hit vector.

	      std::vector< art::Ptr<recob::Hit> > jhits = fm.at(jclus);
	      unsigned int njhits = jhits.size();
	      assert(hits.size() >= nihits);
	      //hits.resize(nihits);
	      while(hits.size() > nihits)
		hits.pop_back();
	      assert(hits.size() == nihits);
	      hits.reserve(nihits + njhits);
	      for(std::vector< art::Ptr<recob::Hit> >::const_iterator j = jhits.begin();
		  j != jhits.end(); ++j)
		hits.push_back(*j);
	  

	      // Loop over third cluster.

	      for(int kclus = 0; kclus < jclus; ++kclus) {
		art::Ptr<recob::Cluster> pkclus = Clusters.at(kclus);
		geo::View_t kview = pkclus->View();

		// Test third view.

		if(((kview == geo::kU && fSptalg.enableU()) ||
		    (kview == geo::kV && fSptalg.enableV()) ||
		    (kview == geo::kW && fSptalg.enableW()))
		   && kview != iview && kview != jview) {

		  // Store hits from third view into hit vector.

		  std::vector< art::Ptr<recob::Hit> > khits = fm.at(kclus);
		  unsigned int nkhits = khits.size();
		  assert(hits.size() >= nihits + njhits);
		  //hits.resize(nihits + njhits);
		  while(hits.size() > nihits + njhits)
		    hits.pop_back();
		  assert(hits.size() == nihits + njhits);
		  hits.reserve(nihits + njhits + nkhits);
		  for(std::vector< art::Ptr<recob::Hit> >::const_iterator k = khits.begin();
		      k != khits.end(); ++k)
		    hits.push_back(*k);

		  // Make three-view space points.

		  std::vector<recob::SpacePoint> spts;
		  fSptalg.makeSpacePoints(hits, spts);

		  if(spts.size() > 0) {
		    SpacePointVectors.push_back(spts);
		  }
		}
	      }
	    }
	  }
	}
      }

  
    return SpacePointVectors;
  }


  // Extract vector of hits from event

  //----------------------------------------------------------------------------
  art::PtrVector<recob::Hit> SeedFinderModule::GetHitsFromEvent(std::string HitModuleLabel, art::Event & evt)
  {
    art::PtrVector<recob::Hit> TheHits;
    art::Handle< std::vector<recob::Hit> > hith;
    evt.getByLabel(HitModuleLabel, hith);
    for(unsigned int i=0; i<hith->size(); ++i)
      {
	art::Ptr<recob::Hit> hit(hith, i);
	TheHits.push_back(hit);
      }
    
    return TheHits;
  }

  




  //----------------------------------------------------------------------------
  void SeedFinderModule::endJob()
  {

  }
}
