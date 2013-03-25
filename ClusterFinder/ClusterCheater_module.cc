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
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "Utilities/AssociationUtil.h"
#include "Simulation/EmEveIdCalculator.h"
#include "ClusterFinder/HoughBaseAlg.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Utilities/GeometryUtilities.h"


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
    
    std::string fMCGeneratorLabel;  ///< label for module to get MC truth information
    std::string fHitModuleLabel;    ///< label for module creating recob::Hit objects	   
    std::string fG4ModuleLabel;     ///< label for module running G4 and making particles, etc

    HoughBaseAlg fHLAlg;            ///< object holding algorithm to do Hough transforms

  };
}

namespace cluster{

  //--------------------------------------------------------------------
  ClusterCheater::ClusterCheater(fhicl::ParameterSet const& pset)
    : fHLAlg(pset.get< fhicl::ParameterSet >("HoughBaseAlg"))
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
    fMCGeneratorLabel  = pset.get< std::string >("MCGeneratorLabel",  "generator");
    fHitModuleLabel    = pset.get< std::string >("HitModuleLabel",    "hit"     );
    fG4ModuleLabel     = pset.get< std::string >("G4ModuleLabel",     "largeant");
    fHLAlg.reconfigure(pset.get< fhicl::ParameterSet >("HoughBaseAlg"));

    return;
  }

  //--------------------------------------------------------------------
  void ClusterCheater::produce(art::Event& evt)
  {
    art::ServiceHandle<geo::Geometry>      geo;
    art::ServiceHandle<cheat::BackTracker> bt;
    
    // ###################################
    // ### Getting Detector Properties ###
    // ###################################
    art::ServiceHandle<util::DetectorProperties> detp;
    
    // #######################################
    // ### Getting MC Truth Info from simb ###
    // #######################################
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    evt.getByLabel(fMCGeneratorLabel,mctruthListHandle);
    
    // ##################################
    // ### Getting MC Truth simb Info ###
    // ##################################
    ///\todo: Fix the following loop to work properly in case of multiple MCTruths per Event 
    ///\todo: and in case of multiple particles per interaction
    float vertex[5] = {0.};
    for (size_t ii = 0; ii <  mctruthListHandle->size(); ++ii){
      for(int i = 0; i < mctruthListHandle->at(ii).NParticles(); ++i){
	simb::MCParticle const& p = mctruthListHandle->at(ii).GetParticle(i); 
	
	vertex[0] = p.Vx();
	vertex[1] = p.Vy();
	vertex[2] = p.Vz();
	std::cout<<"Vertex X = "<<vertex[0]<<std::endl;
	std::cout<<"Vertex Y = "<<vertex[1]<<std::endl;
	std::cout<<"Vertex Z = "<<vertex[2]<<std::endl;
      }
    }//<---end MC Truth Vertex
    
    // grab the hits that have been reconstructed
    art::Handle< std::vector<recob::Hit> > hitcol;
    evt.getByLabel(fHitModuleLabel, hitcol);

    // make a vector of them - we aren't writing anything out to a file
    // so no need for a art::PtrVector here
    std::vector< art::Ptr<recob::Hit> > hits;
    art::fill_ptr_vector(hits, hitcol);
    
    // adopt an EmEveIdCalculator to find the eve ID.  
    // will return a primary particle if it doesn't find 
    // a responsible particle for an EM process
    bt->SetEveIdCalculator(new sim::EmEveIdCalculator);

    LOG_DEBUG("ClusterCheater") << bt->ParticleList();

    // make a map of vectors of art::Ptrs keyed by eveID values
    std::map< int, std::vector< art::Ptr<recob::Hit> > > eveHitMap;

    // loop over all hits and fill in the map
    for( auto const& itr : hits ){

      std::vector<cheat::TrackIDE> eveides = bt->HitToEveID(itr);

      // loop over all eveides for this hit
      for(size_t e = 0; e < eveides.size(); ++e){

	// don't worry about eve particles that contribute less than 10% of the
	// energy in the current hit
	if( eveides[e].energyFrac < 0.1) continue;

	eveHitMap[eveides[e].trackID].push_back(itr);

      } // end loop over eve IDs for this hit

    }// end loop over hits

    // loop over the map and make clusters
    std::unique_ptr< std::vector<recob::Cluster> > clustercol(new std::vector<recob::Cluster>);
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);

    for(auto const& hitMapItr : eveHitMap){
      unsigned int vtx_wire = 0;

      // separate out the hits for each particle into the different views
      std::vector< art::Ptr<recob::Hit> > eveHits( hitMapItr.second );

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
            	    
	    for(auto const& h : eveHits){
	      
	      geo::WireID wid = h->WireID();

	      if(wid.Plane    != pl || 
		 wid.TPC      != t  || 
		 wid.Cryostat != c) continue;
	      // ================================================================================
	      // === Only keeping clusters with 5 or more hits to help cut down on the number ===
	      // ===                   of clusters found by the algorithm                     ===
	      // ================================================================================
	      if(eveHits.size() < 5){continue;}
	      
	      
	      ptrvs.push_back(h);
	      totalQ += h->Charge();
	      
	      if(wid.Wire < startWire){
		startWire = wid.Wire;
		startTime = 1.e6;
	      }
	      if(wid.Wire > endWire  ){
		endWire = wid.Wire;
		endTime = -1.e-6;
	      }
	      
	      if(wid.Wire == startWire && h->StartTime() < startTime) 
		startTime = h->StartTime();
	      
	      if(wid.Wire == endWire   && h->EndTime()   > endTime  ) 
		endTime   = h->EndTime();
	      
	    } // end loop over hits for this particle	

	    // do not create clusters with zero size hit arrays, Andrzej
	    if(ptrvs.size() == 0) continue;
	    
	    // figure out the rest of the cluster information using these hits
	    
	    // use the HoughLineAlg to get dTdW for these hits
	    double intercept = 0.;
	    fHLAlg.Transform(eveHits, dTdW, intercept);
	    
	     // ==============================================================================
    	    // Translating the truth vertex xyz information into truth wire and time position
            // ===============================================================================
	    unsigned int vtx_channel = geo->NearestChannel(vertex, pl, t, c);
     	    std::vector<geo::WireID> vtxwire = geo->ChannelToWire(vtx_channel);
	    geo::WireID wgeom = vtxwire[0];
				
	    vtx_wire = wgeom.Wire;
	     
	    // === If the cluster is to the downstream of the vertex the positions are switched ===
	    if(startWire - vtx_wire < 0 && endWire - vtx_wire < 0){
	      float tempstartWire = endWire;
	      float tempendWire = startWire;
	      
	      startWire = tempstartWire;
	      endWire = tempendWire;
	      
	      float tempstartTime = endTime;
	      float tempendTime = startTime;
	      
	      startTime = tempstartTime;
	      endTime = tempendTime;
	    }
	    
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
						 (hitMapItr.first * 1000) + pl*100 + t*10 + c));
	    
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


