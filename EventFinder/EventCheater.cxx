////////////////////////////////////////////////////////////////////////
/// \file    EventCheater.cxx
/// \brief   create perfectly reconstructed events
/// \version $Id: GeometryTest.cxx,v 1.1 2011/02/17 01:45:48 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#include <vector>

// ROOT includes

// LArSoft includes
#include "MCCheater/BackTracker.h"
#include "EventFinder/EventCheater.h"
#include "Utilities/AssociationUtil.h"
#include "RecoBase/recobase.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/FindOne.h"

namespace event{

  //--------------------------------------------------------------------
  EventCheater::EventCheater(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    produces< std::vector<recob::Event> >();
    produces< art::Assns<recob::Event, recob::Vertex> >();
    produces< art::Assns<recob::Event, recob::Hit> >();
  }

  //--------------------------------------------------------------------
  EventCheater::~EventCheater()
  {
  }

  //--------------------------------------------------------------------
  void EventCheater::reconfigure(fhicl::ParameterSet const& pset)
  {
    fCheatedVertexLabel = pset.get< std::string >("CheatedVertexLabel", "prong" );
    fG4ModuleLabel      = pset.get< std::string >("G4ModuleLabel",      "largeant");

    return;
  }

  //--------------------------------------------------------------------
  void EventCheater::produce(art::Event& evt)
  {

    mf::LogError("EventCheater") << "the EventCheater is temporarily broken, please don't use it";

    // grab the sim::ParticleList
    sim::ParticleList plist = sim::SimListUtils::GetParticleList(evt, fG4ModuleLabel);

    art::Handle< std::vector<sim::Particle> > pcol;
    evt.getByLabel(fG4ModuleLabel, pcol);

    // grab the vertices that have been reconstructed
    art::Handle< std::vector<recob::Vertex> > vertexcol;
    evt.getByLabel(fCheatedVertexLabel, vertexcol);

    // make a vector of them - we aren't writing anything out to a file
    // so no need for a art::PtrVector here
    std::vector< art::Ptr<recob::Vertex> > vertices;
    art::fill_ptr_vector(vertices, vertexcol);

    // loop over the vertices and figure out which primaries they are associated with
    std::vector< art::Ptr<recob::Vertex> >::iterator vertexitr = vertices.begin();

    // make a map of primary product id's to collections of vertices
    std::map<art::ProductID, std::vector< art::Ptr<recob::Vertex> > > vertexMap;
    std::map<art::ProductID, std::vector< art::Ptr<recob::Vertex> > >::iterator vertexMapItr = vertexMap.begin();

    // loop over all prongs
    while( vertexitr != vertices.end() ){

      // the ID of the vertex should be a primary particle
      // from the event generator so the mother has to be 0
      art::FindOne<simb::MCTruth> fomct(pcol, evt, fG4ModuleLabel);
      //art::Ptr<simb::MCTruth> mctp;
      for(size_t p = 0; p < pcol->size(); ++p){
	art::Ptr<sim::Particle> pptr(pcol, p);
	if( pptr->TrackId() == (*vertexitr)->ID() ){
	  // do something here with the FindOne object to get the
	  // product ID of the MCTruth object.  Probably requires
	  // waiting for an ART version bump from v1.00.05
	}
      }
      
      art::ProductID primary;// = mctp.id();

      if(vertexMap.find(primary) != vertexMap.end()){
	  ((*vertexMapItr).second).push_back((*vertexitr));
      }
      else{
	std::vector< art::Ptr<recob::Vertex> > vertexvec;
	vertexvec.push_back(*vertexitr);
	vertexMap[primary] = vertexvec;
      }

      vertexitr++;
    }// end loop over vertices

    std::auto_ptr< std::vector<recob::Event> > eventcol(new std::vector<recob::Event>);
    std::auto_ptr< art::Assns<recob::Event, recob::Vertex> > evassn(new art::Assns<recob::Event, recob::Vertex>);
    std::auto_ptr< art::Assns<recob::Event, recob::Hit> > ehassn(new art::Assns<recob::Event, recob::Hit>);

    // loop over the map and associate all vertex objects with an event
    for(vertexMapItr = vertexMap.begin(); vertexMapItr != vertexMap.end(); vertexMapItr++){

      // Vertex objects require PtrVectors of showers and tracks as well
      // as a vertex position for their constructor
      art::PtrVector<recob::Vertex> ptrvs;

      std::vector< art::Ptr<recob::Vertex> > verts( (*vertexMapItr).second );

      for(size_t v = 0; v < verts.size(); ++v)
	ptrvs.push_back(verts[v]);
	
      // add an event to the collection.  
      eventcol->push_back(recob::Event(ptrvs, (*vertexMapItr).first.productIndex()));

      // associate the event with its vertices
      util::CreateAssn(*this, evt, *(eventcol.get()), ptrvs, *(evassn.get()));
      
      // get the hits associated with each vertex and associate those with the event
      for(size_t p = 0; p < ptrvs.size(); ++p){
	art::PtrVector<recob::Hit> hits = util::FindManyP<recob::Hit>(ptrvs, evt, fCheatedVertexLabel, p);
	util::CreateAssn(*this, evt, *(eventcol.get()), hits, *(ehassn.get()));
      }


      mf::LogInfo("EventCheater") << "adding event: \n" 
				  << eventcol->back()
				  << "\nto collection";

    } // end loop over the map

    evt.put(eventcol);
    evt.put(evassn);
    evt.put(ehassn);

    return;

  } // end produce

} // end namespace
