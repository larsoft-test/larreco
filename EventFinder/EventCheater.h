///////////////////////////////////////////////////////////////////////
/// \file    EventCheater.h
/// \brief   make events using MC truth information
/// \author  brebel@fnal.gov
/// \version $Id: GeometryTest.h,v 1.1 2011/02/17 01:45:48 brebel Exp $
///////////////////////////////////////////////////////////////////////
#ifndef EVENT_EVENTCHEATER_H
#define EVENT_EVENTCHEATER_H
#include <string>

#include "art/Framework/Core/EDProducer.h"

namespace event {
  class EventCheater : public art::EDProducer {
  public:
    explicit EventCheater(fhicl::ParameterSet const& pset);
    virtual ~EventCheater();

    void produce(art::Event& evt);

    void reconfigure(fhicl::ParameterSet pset);

  private:

    std::string fCheatedVertexLabel; ///< label for module creating recob::Vertex objects	   
    std::string fG4ModuleLabel;      ///< label for module running G4 and making particles, etc
    std::string fDetSimModuleLabel;  ///< label for module creating sim::SimChannel objects

  };
}
#endif
