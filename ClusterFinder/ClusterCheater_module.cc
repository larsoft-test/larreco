////////////////////////////////////////////////////////////////////////
// $Id: ClusterCheater_module.cc Exp $
//
// ClusterCheater module
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef CLUSTER_CLUSTERCHEATER_H
#define CLUSTER_CLUSTERCHEATER_H
#include <string>
#include <vector>

// ROOT includes

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "MCCheater/BackTracker.h"
#include "ClusterFinder/ClusterCheater.h"
#include "RecoBase/recobase.h"
#include "Utilities/AssociationUtil.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"
#include "ClusterFinder/HoughBaseAlg.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"

namespace cluster {
  class ClusterCheater : public art::EDProducer {
  public:
    explicit ClusterCheater(fhicl::ParameterSet const& pset);
    virtual ~ClusterCheater();

    void produce(art::Event& evt);

    void reconfigure(fhicl::ParameterSet const& pset);

 private:

    std::string fHitModuleLabel;    ///< label for module creating recob::Hit objects	   
    std::string fG4ModuleLabel;     ///< label for module running G4 and making particles, etc

    HoughBaseAlg fHLAlg;            ///< object holding algorithm to do Hough transforms

  };
}

namespace cluster{

  //--------------------------------------------------------------------
  ClusterCheater::ClusterCheater(fhicl::ParameterSet const& pset)
    : fHLAlg(pset.get< fhicl::ParameterSet >("HoughLineAlg"))
  {
    this->reconfigure(pset);

    produces< std::vector<recob::Cluster> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
  }

  //--------------------------------------------------------------------
  ClusterCheater::~ClusterCheater()
  {
  }

  //--------------------------------------------------------------------
  void ClusterCheater::reconfigure(fhicl::ParameterSet const& pset)
  {
    fHitModuleLabel    = pset.get< std::string >("HitModuleLabel",    "hit"     );
    fG4ModuleLabel     = pset.get< std::string >("G4ModuleLabel",     "largeant");
    fHLAlg.reconfigure(pset.get< fhicl::ParameterSet >("HoughLineAlg"));

    return;
  }

  //--------------------------------------------------------------------
  void ClusterCheater::produce(art::Event& evt)
  {
    art::ServiceHandle<geo::Geometry>      geo;
    art::ServiceHandle<cheat::BackTracker> bt;

    // grab the hits that have been reconstructed
    art::Handle< std::vector<recob::Hit> > hitcol;
    evt.getByLabel(fHitModuleLabel, hitcol);

    // make a vector of them - we aren't writing anything out to a file
    // so no need for a art::PtrVector here
    std::vector< art::Ptr<recob::Hit> > hits;
    art::fill_ptr_vector(hits, hitcol);
    
    // loop over the hits and figure out which particle contributed to each one
    std::vector< art::Ptr<recob::Hit> >::iterator itr = hits.begin();

    // adopt an EmEveIdCalculator to find the eve ID.  
    // will return a primary particle if it doesn't find 
    // a responsible particle for an EM process
    bt->SetEveIdCalculator(new sim::EmEveIdCalculator);

    mf::LogVerbatim("ClusterCheater") << bt->ParticleList();

    // make a map of vectors of art::Ptrs keyed by eveID values
    std::map< int, std::vector< art::Ptr<recob::Hit> > > eveHitMap;
    std::map< int, std::vector< art::Ptr<recob::Hit> > >::iterator hitMapItr = eveHitMap.begin();

    // loop over all hits and fill in the map
    while( itr != hits.end() ){

      std::vector<cheat::TrackIDE> eveides = bt->HitToEveID(*itr);

      // loop over all eveides for this hit
      for(size_t e = 0; e < eveides.size(); ++e){

	// don't worry about eve particles that contribute less than 10% of the
	// energy in the current hit
	if( eveides[e].energyFrac < 0.1) continue;

	hitMapItr = eveHitMap.find( eveides[e].trackID );
	
	// is this id already in the map, if so extend the collection 
	// by one hit, otherwise make a new collection and put it in
	// the map
	if( hitMapItr != eveHitMap.end() ){
	  ((*hitMapItr).second).push_back((*itr));
	}
	else{
	  std::vector< art::Ptr<recob::Hit> > hitvec;
	  hitvec.push_back(*itr);
	  eveHitMap[eveides[e].trackID] = hitvec;
	}

      } // end loop over eve IDs for this hit

      itr++;
    }// end loop over hits

    // loop over the map and make clusters
    std::unique_ptr< std::vector<recob::Cluster> > clustercol(new std::vector<recob::Cluster>);
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);

    unsigned int plane = 0;
    unsigned int wire  = 0;
    unsigned int tpc   = 0;
    unsigned int cstat = 0;
    for(hitMapItr = eveHitMap.begin(); hitMapItr != eveHitMap.end(); hitMapItr++){

      // separate out the hits for each particle into the different views
      std::vector< art::Ptr<recob::Hit> > eveHits( (*hitMapItr).second );

      for(size_t c = 0; c < geo->Ncryostats(); ++c){
	for(size_t t = 0; t < geo->Cryostat(c).NTPC(); ++t){
	  for(size_t pl = 0; pl < geo->Cryostat(c).TPC(t).Nplanes(); ++pl){
	    art::PtrVector<recob::Hit> ptrvs;
	    double startWire = 1.e6;
	    double startTime = 1.e6;
	    double endWire   = -1.e6;
	    double endTime   = -1.e6;
	    double dTdW      = 0.;
	    double dQdW      = 0.;
	    double totalQ    = 0.;
	    geo::View_t view = geo->Cryostat(c).TPC(t).Plane(pl).View();

	    for(size_t h = 0; h < eveHits.size(); ++h){
	      
	      geo->ChannelToWire(eveHits[h]->Channel(), cstat, tpc, plane, wire);
	      
	      if(plane != pl || tpc != t || cstat != c) continue;
	      
	      ptrvs.push_back(eveHits[h]);
	      totalQ += eveHits[h]->Charge();
	      
	      if(wire < startWire){
		startWire = wire;
		startTime = 1.e6;
	      }
	      if(wire > endWire  ){
		endWire = wire;
		endTime = -1.e-6;
	      }
	      
	      if(wire == startWire && eveHits[h]->StartTime() < startTime) 
		startTime = eveHits[h]->StartTime();
	      
	      if(wire == endWire   && eveHits[h]->EndTime()   > endTime  ) 
		endTime   = eveHits[h]->EndTime();
	      
	    } // end loop over hits for this particle	

	    // do not create clusters with zero size hit arrays, Andrzej
	    if(ptrvs.size()==0)
	      continue;
	    
	    // figure out the rest of the cluster information using these hits
	    
	    // use the HoughLineAlg to get dTdW for these hits
	    double intercept = 0.;
	    fHLAlg.Transform(eveHits, dTdW, intercept);
	    
	    ///\todo now figure out the dQdW
	    
	    // add a cluster to the collection.  Make the ID be the eve particle
	    // trackID*1000 + plane number*100 + tpc that the current hits are from
	    
	    clustercol->push_back(recob::Cluster(startWire, 0.,
						 startTime, 0.,
						 endWire,   0.,
						 endTime,   0.,
						 dTdW,      0.,
						 dQdW,      0.,
						 totalQ,
						 view,
						 ((*hitMapItr).first * 1000) + pl*100 + tpc*10 + cstat));
	    
	    // association the hits to this cluster
	    util::CreateAssn(*this, evt, *(clustercol.get()), ptrvs, *(assn.get()));
	    
	    mf::LogInfo("ClusterCheater") << "adding cluster: \n" 
					  << clustercol->back()
					  << "\nto collection.";

	  } // end loop over the number of planes
	} // end loop over the tpcs
      } // end loop over the cryostats
    } // end loop over the map

    evt.put(std::move(clustercol));
    evt.put(std::move(assn));

    return;

  } // end produce

} // end namespace

namespace cluster{

  DEFINE_ART_MODULE(ClusterCheater);

}

#endif


