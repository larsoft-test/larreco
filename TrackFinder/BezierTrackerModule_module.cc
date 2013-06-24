#include "art/Persistency/Common/PtrVector.h"

#ifndef BEZIERTRACKERMOD_H
#define BEZIERTRACKERMOD_H

//
// Name: BezierTrackerModule.h
//
// Purpose: Header file for module BezierTrackerModule.  This modules makes
//          bezier tracks out of seed collections
//
// Configuration parameters.
//
// SeedModuleLabel;     // Cluster module label (e.g. "dbcluster").
// HitModuleLabel;      // Hit module label (e.g. "FFTHitFinder")
//
// Ben Jones, MIT
//

#include "art/Framework/Core/EDProducer.h"
#include "RecoAlg/SeedFinderAlgorithm.h"
#include "Geometry/Geometry.h"

namespace recob
{
  class Seed;
  class Track;
  class Hit;
}


namespace trkf {

  class BezierTrack;
  class BezierTrackerAlgorithm;
  class SpacePointAlg;

  class BezierTrackerModule : public art::EDProducer
  {
  public:
 
    // Constructors, destructor

    explicit BezierTrackerModule(fhicl::ParameterSet const& pset);
    virtual ~BezierTrackerModule();

    
    // Overrides.

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void produce(art::Event& evt);
    void endJob();
    
    
    

  private:

    // Fcl Attributes.

    std::string fSeedModuleLabel;
    std::string fHitModuleLabel;
    std::string fClusterModuleLabel;

    int fTrackMode;

    trkf::BezierTrackerAlgorithm * fBTrackAlg;
    
    std::vector<std::vector<recob::SpacePoint> > GetSpacePointsFromClusters(std::string ClusterModuleLabel, art::Event& evt);

    std::vector<std::vector<recob::Seed> > GetSeedsFromClusterHits(std::map<geo::View_t, std::vector<art::PtrVector<recob::Hit> > > &);


    std::vector<std::vector<recob::Seed> > GetSeedsFromClusters(std::string ClusterModuleLabel, art::Event& evt);
    
    std::map<geo::View_t, std::vector<art::PtrVector<recob::Hit> > > GetHitsFromClusters(std::string ClusterModuleLabel, art::Event& evt);
    
    void GetHitsToAssn(std::map<geo::View_t, std::vector<art::PtrVector<recob::Hit> > >& SortedHits, int ClusterNum,  BezierTrack const& BTrack,  std::vector<double>& ReturnSValues, art::PtrVector<recob::Hit>& ReturnVec);
    
  };
}

#endif 

#include "art/Framework/Core/ModuleMacros.h" 


namespace trkf {
  DEFINE_ART_MODULE(BezierTrackerModule);
}

//
// Name: BezierTrackerModule.cxx
//
// Purpose: Implementation file for module BezierTrackerModule.
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
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoObjects/BezierTrack.h"
#include "Utilities/AssociationUtil.h"
#include "RecoAlg/BezierTrackerAlgorithm.h"

namespace trkf {

  BezierTrackerModule::BezierTrackerModule(const fhicl::ParameterSet& pset)
  {
    reconfigure(pset);
    produces< std::vector<recob::Track> >();
    produces< std::vector<recob::Seed> >();
    produces< art::Assns<recob::Track, recob::Hit> >();
      
  }

  BezierTrackerModule::~BezierTrackerModule()
  {
  }

  void BezierTrackerModule::reconfigure(fhicl::ParameterSet const& pset)
  {
    fSeedModuleLabel   = pset.get<std::string>("SeedModuleLabel");
    fClusterModuleLabel= pset.get<std::string>("ClusterModuleLabel");
    fHitModuleLabel    = pset.get<std::string>("HitModuleLabel");
    fTrackMode         = pset.get<double>("TrackMode");
      
    fBTrackAlg = new trkf::BezierTrackerAlgorithm(pset.get<fhicl::ParameterSet>("BezierTrackerAlgorithm"));
      
}

  void BezierTrackerModule::beginJob()
  {}


  void BezierTrackerModule::produce(art::Event& evt)
  {
 
    
    // Extract hits PtrVector from event

    art::Handle< std::vector<recob::Hit> > hith;
    evt.getByLabel(fHitModuleLabel, hith);

    std::vector<art::Ptr<recob::Hit> > HitVec;
    for(unsigned int i=0; i < hith->size(); ++i)
      {
	art::Ptr<recob::Hit> hit(hith,i);
	HitVec.push_back(hit);
      }
    

    // Declare products to store

    std::unique_ptr< std::vector<recob::Track > > btracks ( new std::vector<recob::Track>);
    std::unique_ptr< std::vector<recob::Seed > > seeds ( new std::vector<recob::Seed>);
    std::unique_ptr< art::Assns<recob::Track, recob::Hit > > assn( new art::Assns<recob::Track, recob::Hit>);
   

    std::vector<trkf::BezierTrack >           BTracks;
    
    std::vector<art::PtrVector<recob::Hit> >  HitsForAssns;
    std::vector<std::vector<double> >         SValuesOfHits;
    
   
    if(fTrackMode==1)
      {
	// Look for track-like features in seed collections
	art::Handle< std::vector<recob::Seed> > seedh;
        evt.getByLabel(fSeedModuleLabel, seedh);

	std::vector<recob::Seed> TrackSeeds(*seedh);

        fBTrackAlg->MakeBezierTracksFromSeeds(BTracks, TrackSeeds);
	// Insert hit collecting code here
      }

   
    else if(fTrackMode==2)
      {
	// Find tracks in amorphous hit collection
        fBTrackAlg->MakeBezierTracksFromHits(BTracks, HitVec, HitsForAssns);
      }

    else if(fTrackMode==3)
      {
	// Find tracks from cluster combinations
	//	mf::LogVerbatim("BezierTrackerModule")<<"Bezier tracker configured in mode 3, building tracks from cluster combinations"<<std::endl;

	std::map<geo::View_t, std::vector<art::PtrVector<recob::Hit> > > SortedHits = 
	  GetHitsFromClusters(fClusterModuleLabel, evt);
	
	std::vector<std::vector<recob::Seed> > Seeds = GetSeedsFromClusterHits(SortedHits);
	
	for(size_t i=0; i!=Seeds.size(); ++i)
	  {		  
	    for(size_t j=0; j!=Seeds.at(i).size(); ++j)
	      {
		seeds->push_back(Seeds.at(i).at(j));
	      }
	    std::vector<trkf::BezierTrack> BTracksThisCombo;
	    if(Seeds.at(i).size()>0)
	      fBTrackAlg->MakeBezierTracksFromSeeds(BTracksThisCombo, Seeds.at(i));
	
	    for(size_t j=0; j!=BTracksThisCombo.size(); ++j)
	      {

		BTracks.push_back(BTracksThisCombo.at(j));
		size_t ThisIndex = BTracks.size()-1;
		
		
		SValuesOfHits.push_back(std::vector<double>());
		HitsForAssns.push_back(art::PtrVector<recob::Hit>());

	
		GetHitsToAssn(SortedHits,i,BTracksThisCombo.at(j), SValuesOfHits.at(ThisIndex), HitsForAssns.at(ThisIndex));
		
		mf::LogVerbatim("BezierTracker") << " Bezier track " << ThisIndex<<" has hit collection size " << HitsForAssns.at(ThisIndex).size();
			
	      }
	    BTracksThisCombo.clear();	    
	  }
	for(size_t i=0; i!=Seeds.size(); ++i)
	  Seeds.at(i).clear();
	Seeds.clear();
	for(std::map<geo::View_t, std::vector<art::PtrVector<recob::Hit> > >::iterator it=SortedHits.begin();
	    it!=SortedHits.end();++it)
	  {
	    for(size_t i=0; i!=it->second.size(); ++i)
	      it->second.at(i).clear();
	    it->second.clear();
	  }
	SortedHits.clear();
	  
       }

    fBTrackAlg->MakeDirectJoins(BTracks, HitsForAssns);
    

    for(size_t i=0; i!=BTracks.size(); ++i)
      {
	
	std::unique_ptr<recob::Track>  ToStore = BTracks.at(i).GetBaseTrack();
	btracks->push_back(*ToStore);
	util::CreateAssn(*this, evt, *(btracks.get()), HitsForAssns.at(i), *(assn.get()));
      }
   
    mf::LogVerbatim("BezierTrackerAlgorithm")<<"Storing in evt - check"<<std::endl;
    evt.put(std::move(btracks));
    evt.put(std::move(seeds));
    evt.put(std::move(assn));

    // Now tidy up
    trkf::SpacePointAlg *Sptalg = fBTrackAlg->GetSeedFinderAlgorithm()->GetSpacePointAlg();
    Sptalg->clearHitMap();
    for(size_t i=0; i!=HitsForAssns.size(); ++i)
      HitsForAssns.at(i).clear();
    HitsForAssns.clear();
    

  }
  //-----------------------------------------

  std::map<geo::View_t, std::vector<art::PtrVector<recob::Hit> > > BezierTrackerModule::GetHitsFromClusters(std::string ClusterModuleLabel, art::Event& evt)
  {
    
    std::map<geo::View_t, std::vector<art::PtrVector<recob::Hit> > > ReturnVec;
    
    std::vector<art::Ptr<recob::Cluster> > Clusters;
    
    art::Handle< std::vector<recob::Cluster> > clusterh;
    evt.getByLabel(ClusterModuleLabel, clusterh);
    
    if(clusterh.isValid()) {
      art::fill_ptr_vector(Clusters, clusterh);
    }
    
    art::FindManyP<recob::Hit> fm(clusterh, evt, ClusterModuleLabel);
      
    for(size_t iclus = 0; iclus < Clusters.size(); ++iclus) {
      art::Ptr<recob::Cluster> ThisCluster = Clusters.at(iclus);
      
      std::vector< art::Ptr<recob::Hit> > ihits = fm.at(iclus);
      
      art::PtrVector<recob::Hit> HitsThisCluster;    
      for(std::vector< art::Ptr<recob::Hit> >::const_iterator i = ihits.begin();
	  i != ihits.end(); ++i)
	HitsThisCluster.push_back(*i);
      
      ReturnVec[ThisCluster->View()].push_back(HitsThisCluster);
    }
    return ReturnVec;
  }
  
  //----------------------------------------------


  std::vector<std::vector<recob::Seed> > BezierTrackerModule::GetSeedsFromClusterHits(std::map<geo::View_t, std::vector<art::PtrVector<recob::Hit> > >& SortedHits)
    {
      trkf::SpacePointAlg *Sptalg = fBTrackAlg->GetSeedFinderAlgorithm()->GetSpacePointAlg();
      
      std::vector<std::vector<recob::Seed> > ReturnVec;
      
      // This piece of code looks detector specific, but its not - 
      //   it also works for 2 planes, but one vector is empty.
      
      if(!(Sptalg->enableU()&&Sptalg->enableV()&&Sptalg->enableW()))
	mf::LogWarning("BezierTrackerModule")<<"Warning: SpacePointAlg is does not have three views enabled. This may cause unexpected behaviour in the bezier tracker.";
      
      for(std::vector<art::PtrVector<recob::Hit> >::iterator itU = SortedHits[geo::kU].begin();
	  itU !=SortedHits[geo::kU].end(); ++itU)
	for(std::vector<art::PtrVector<recob::Hit> >::iterator itV = SortedHits[geo::kV].begin();
	    itV !=SortedHits[geo::kV].end(); ++itV)
	  for(std::vector<art::PtrVector<recob::Hit> >::iterator itW = SortedHits[geo::kW].begin();
	      itW !=SortedHits[geo::kW].end(); ++itW)
	    {
	      art::PtrVector<recob::Hit> HitsFromThisCombo;
	    
	      if(Sptalg->enableU()) 
		for(size_t i=0; i!=itU->size(); ++i) 
		  HitsFromThisCombo.push_back(itU->at(i));
	
	      if(Sptalg->enableV()) 
		for(size_t i=0; i!=itV->size(); ++i) 
		  HitsFromThisCombo.push_back(itV->at(i));
	     
	      if(Sptalg->enableW())
		for(size_t i=0; i!=itW->size(); ++i) 
		  HitsFromThisCombo.push_back(itW->at(i));
 	    
	      std::vector<recob::SpacePoint> spts;
	      Sptalg->makeSpacePoints(HitsFromThisCombo, spts);
	    
	      
	      if(spts.size()>0)
		{
		  std::vector<std::vector<recob::SpacePoint> > CataloguedSPs;
		  
		  std::vector<recob::Seed> Seeds 
		    = fBTrackAlg->GetSeedFinderAlgorithm()->FindSeeds(spts,CataloguedSPs);
		  
		  ReturnVec.push_back(Seeds);
		  
		  spts.clear();
		}
	      
	      else
		{
		  ReturnVec.push_back(std::vector<recob::Seed>());
		}
	    }
      
      return ReturnVec;
      
    }


  //---------------------------------------------

  std::vector<std::vector<recob::Seed> > BezierTrackerModule::GetSeedsFromClusters(std::string ClusterModuleLabel, art::Event& evt)
  {
    // Get Services.
    trkf::SpacePointAlg *Sptalg = fBTrackAlg->GetSeedFinderAlgorithm()->GetSpacePointAlg();

    std::vector<std::vector<recob::Seed> > ReturnVec;
    
    art::ServiceHandle<geo::Geometry> geom;

    std::vector<art::Ptr<recob::Cluster> > Clusters;

    art::Handle< std::vector<recob::Cluster> > clusterh;
    evt.getByLabel(ClusterModuleLabel, clusterh);

    if(clusterh.isValid()) {
      art::fill_ptr_vector(Clusters, clusterh);
    }
    art::FindManyP<recob::Hit> fm(clusterh, evt, ClusterModuleLabel);



    // Make a double or triple loop over clusters in distinct views
    // (depending on minimum number of views configured in SpacePointAlg).
    art::PtrVector<recob::Hit> hits;

    // Loop over first cluster.

    int nclus = Clusters.size();
    mf::LogVerbatim("BezierTrackerModule")<< "There are " << nclus<< " clusters in the event";
    for(int iclus = 0; iclus < nclus; ++iclus) {
      art::Ptr<recob::Cluster> piclus = Clusters.at(iclus);
      geo::View_t iview = piclus->View();


      // Test first view.

      //    mf::LogVerbatim("BezierTrackerModule") << "View check: " << iview << " " << Sptalg->enableU()<< " " << Sptalg->enableV() << " " << Sptalg->enableW()<<std::endl;

      if((iview == geo::kU && Sptalg->enableU()) ||
         (iview == geo::kV && Sptalg->enableV()) ||
         (iview == geo::kZ && Sptalg->enableW())) {

        // Store hits from first view into hit vector.

	std::vector< art::Ptr<recob::Hit> > ihits = fm.at(iclus);
        unsigned int nihits = ihits.size();
	//	mf::LogVerbatim("BezierTrackerModule")<<"Cluster " << iclus<< " has " <<nihits<< " hits " <<std::endl;
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

	  if(((jview == geo::kU && Sptalg->enableU()) ||
	      (jview == geo::kV && Sptalg->enableV()) ||
	      (jview == geo::kZ && Sptalg->enableW()))
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

	      if(((kview == geo::kU && Sptalg->enableU()) ||
		  (kview == geo::kV && Sptalg->enableV()) ||
		  (kview == geo::kZ && Sptalg->enableW()))
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
		Sptalg->makeSpacePoints(hits, spts);

		if(spts.size() > 0) 
		  {
		    std::vector<std::vector<recob::SpacePoint> > CataloguedSPs;
		   
		    //  mf::LogVerbatim("BezierTrackerModule") 
		    //   << "Cluster combo found with " << spts.size() 
		    //  << " sps" << std::endl;
		    
		    std::vector<recob::Seed> Seeds 
		      = fBTrackAlg->GetSeedFinderAlgorithm()->FindSeeds(spts,CataloguedSPs);
		    mf::LogVerbatim("BezierTrackerModule") << "Found " << 
		      Seeds.size() << " seeds" << std::endl;
		    
		    if(Seeds.size()>0)
		      {
			ReturnVec.push_back(Seeds);

		      }
		    spts.clear();
		  }
	      
	      }
	    }
	  }
	}
      }
    }


    return ReturnVec;
  }






  std::vector<std::vector<recob::SpacePoint> > BezierTrackerModule::GetSpacePointsFromClusters(std::string ClusterModuleLabel, art::Event& evt)
  {
    // Get Services.
    trkf::SpacePointAlg *Sptalg = fBTrackAlg->GetSeedFinderAlgorithm()->GetSpacePointAlg();

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
    mf::LogVerbatim("BezierTrackerModule")<< "There are " << nclus<< " clusters in the event";
    for(int iclus = 0; iclus < nclus; ++iclus) {
      art::Ptr<recob::Cluster> piclus = Clusters.at(iclus);
      geo::View_t iview = piclus->View();

      // Test first view.

      //      mf::LogVerbatim("BezierTrackerModule") << "View check: " << iview << " " << Sptalg->enableU()<< " " << Sptalg->enableV() << " " << Sptalg->enableW()<<std::endl;

      if((iview == geo::kU && Sptalg->enableU()) ||
         (iview == geo::kV && Sptalg->enableV()) ||
         (iview == geo::kZ && Sptalg->enableW())) {

        // Store hits from first view into hit vector.

	std::vector< art::Ptr<recob::Hit> > ihits = fm.at(iclus);
        unsigned int nihits = ihits.size();
	//	mf::LogVerbatim("BezierTrackerModule")<<"Cluster " << iclus<< " has " <<nihits<< " hits " <<std::endl;
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

	  if(((jview == geo::kU && Sptalg->enableU()) ||
	      (jview == geo::kV && Sptalg->enableV()) ||
	      (jview == geo::kZ && Sptalg->enableW()))
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

	      if(((kview == geo::kU && Sptalg->enableU()) ||
		  (kview == geo::kV && Sptalg->enableV()) ||
		  (kview == geo::kZ && Sptalg->enableW()))
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
		Sptalg->makeSpacePoints(hits, spts);

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



  //---------------------------------------------------------------------
 
  void BezierTrackerModule::GetHitsToAssn(std::map<geo::View_t, std::vector<art::PtrVector<recob::Hit> > >& SortedHits, int ClusterNum,  BezierTrack const& BTrack, std::vector<double>& ReturnSValues, art::PtrVector<recob::Hit>& ReturnVec)
  {
    
    trkf::SpacePointAlg *Sptalg = fBTrackAlg->GetSeedFinderAlgorithm()->GetSpacePointAlg();

    int ComboNumber=0;
    
    if(!(Sptalg->enableU()&&Sptalg->enableV()&&Sptalg->enableW()))
      mf::LogWarning("BezierTrackerModule")<<"Warning: SpacePointAlg is does not have three views enabled. This may cause unexpected behaviour in the bezier tracker.";
    
    for(std::vector<art::PtrVector<recob::Hit> >::iterator itU = SortedHits[geo::kU].begin();
	    itU !=SortedHits[geo::kU].end(); ++itU)
      for(std::vector<art::PtrVector<recob::Hit> >::iterator itV = SortedHits[geo::kV].begin();
	  itV !=SortedHits[geo::kV].end(); ++itV)
	for(std::vector<art::PtrVector<recob::Hit> >::iterator itW = SortedHits[geo::kW].begin();
	    itW !=SortedHits[geo::kW].end(); ++itW)
	  {
	    if(ComboNumber==ClusterNum)
	      {
		art::PtrVector<recob::Hit> HitsFromThisCombo;
		
		if(Sptalg->enableU()) 
		  for(size_t i=0; i!=itU->size(); ++i) 
		    HitsFromThisCombo.push_back(itU->at(i));
		
		if(Sptalg->enableV()) 
		  for(size_t i=0; i!=itV->size(); ++i) 
		    HitsFromThisCombo.push_back(itV->at(i));
		
		if(Sptalg->enableW())
		  for(size_t i=0; i!=itW->size(); ++i) 
		    HitsFromThisCombo.push_back(itW->at(i));
		
		std::vector<int> HitIndices = fBTrackAlg->DetermineNearbyHits(HitsFromThisCombo, BTrack, ReturnSValues);
		
		for(size_t i=0; i!=HitIndices.size(); ++i)
		  {
		    ReturnVec.push_back(HitsFromThisCombo.at(HitIndices.at(i)));
		  }
		
		return;
	      }
	    
	    
	    ComboNumber++;	    
	    
	  }
  }




  //----------------------------------------------------------------------
  void BezierTrackerModule::endJob()
  {

  }
}
