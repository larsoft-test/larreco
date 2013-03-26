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
#include "TStopwatch.h"

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

  struct eveLoc{

    eveLoc(int id, size_t c, size_t t, size_t p)
      : eveID(id)
      , cryostat(c)
      , tpc(t)
      , plane(p)
    {}

    friend bool operator < (eveLoc const& a, eveLoc const& b) { return a.eveID < b.eveID; }
    
    int    eveID;
    size_t cryostat;
    size_t tpc;
    size_t plane;
  };

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

    // make a map of vectors of art::Ptrs keyed by eveID values and 
    // location in cryostat, TPC, plane coordinates of where the hit originated
    std::map< eveLoc, std::vector< art::Ptr<recob::Hit> > > eveHitMap;

    TStopwatch ts;

    // loop over all hits and fill in the map
    for( auto const& itr : hits ){

      std::vector<cheat::TrackIDE> eveides = bt->HitToEveID(itr);

      // loop over all eveides for this hit
      for(size_t e = 0; e < eveides.size(); ++e){

	// don't worry about eve particles that contribute less than 10% of the
	// energy in the current hit
	if( eveides[e].energyFrac < 0.1) continue;

	eveLoc el(eveides[e].trackID, 
		  itr->WireID().Cryostat,
		  itr->WireID().TPC,
		  itr->WireID().Plane);

	eveHitMap[el].push_back(itr);

      } // end loop over eve IDs for this hit

    }// end loop over hits

    // loop over the map and make clusters
    std::unique_ptr< std::vector<recob::Cluster> > clustercol(new std::vector<recob::Cluster>);
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);

    for(auto const& hitMapItr : eveHitMap){

      LOG_DEBUG("ClusterCheater") << "make cluster for eveID: " << hitMapItr.first.eveID
				  << " in cryostat: "           << hitMapItr.first.cryostat
				  << " tpc: "         	        << hitMapItr.first.tpc     
				  << " plane: "       	        << hitMapItr.first.plane;   

      double startWire =  1.e6;
      double startTime =  1.e6;
      double endWire   = -1.e6;
      double endTime   = -1.e6;
      double dTdW      =  0.;
      double dQdW      =  0.;
      double totalQ    =  0.;

      // ==============================================================================
      // Translating the truth vertex xyz information into truth wire and time position
      // ===============================================================================
      unsigned int vtx_channel = geo->NearestChannel(vertex, 
						     hitMapItr.first.plane, 
						     hitMapItr.first.tpc, 
						     hitMapItr.first.cryostat);
      unsigned int vtx_wire = geo->ChannelToWire(vtx_channel)[0].Wire;
      ts.Stop();
      ts.Print();

      // ================================================================================
      // === Only keeping clusters with 5 or more hits to help cut down on the number ===
      // ===                   of clusters found by the algorithm                     ===
      // ================================================================================
      if(hitMapItr.second.size() < 5) continue;
	            
      for(auto const& h : hitMapItr.second){
	      
	geo::WireID wid = h->WireID();

	totalQ += h->Charge();
	      
	if(wid.Wire < startWire){
	  startWire = wid.Wire;
	  startTime = 1.e6;
	}
	if(wid.Wire > endWire  ){
	  endWire = wid.Wire;
	  endTime = -1.e-6;
	}
	
	if(wid.Wire == startWire && h->StartTime() < startTime) startTime = h->StartTime();
	
	if(wid.Wire == endWire   && h->EndTime()   > endTime  ) endTime   = h->EndTime();

      } // end loop over hits for this particle	

      mf::LogVerbatim("ClusterCheater") << "run the hough transform";
      ts.Reset();
      ts.Start();

      // remove need to use hough algorithm and instead use MCTruth info
      // to convert from world coordinates to time and wire number in order
      // to get dTdW
      double intercept = 0.;
      fHLAlg.Transform(hitMapItr.second, dTdW, intercept);
      ts.Stop();
      ts.Print();
      
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
					   hitMapItr.second.at(0)->View(),
					   (hitMapItr.first.eveID*1000 + 
					    hitMapItr.first.plane*100  + 
					    hitMapItr.first.tpc*10     + 
					    hitMapItr.first.cryostat)
					   )
			    );
      
      // association the hits to this cluster
      util::CreateAssn(*this, evt, *clustercol, hitMapItr.second, *assn);
      
      mf::LogInfo("ClusterCheater") << "adding cluster: \n" 
				    << clustercol->back()
				    << "\nto collection.";
      
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


