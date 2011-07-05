////////////////////////////////////////////////////////////////////////
/// \file  AggregateVertex.h
/// \brief Module to find vertices based on 2-d clusters
///
/// echurch@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef AGGREGATEVERTEX_H
#define AGGREGATEVERTEX_H

// Framework includes
#include "art/Framework/Core/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Handle.h" 
#include "art/Persistency/Common/View.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Core/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h" 
#include "art/Framework/Core/EDProducer.h" 

// LArSoft includes
#include "RecoBase/recobase.h"

#include <vector>
#include <string>

namespace vertex {


  class AggregateVertex : public art::EDProducer
  {

  public:

    explicit AggregateVertex(fhicl::ParameterSet const& pset);
    virtual ~AggregateVertex();

    void produce(art::Event& evt); 
    void beginJob(); 

    std::auto_ptr< std::vector<recob::Vertex> >  MatchV2T();

  private:

    std::string fDBScanModuleLabel;
    std::string fHoughModuleLabel;
    std::string fTrack3DModuleLabel;
    std::string fEndPointModuleLabel;

    art::PtrVector<recob::Vertex>  feplist;
    art::PtrVector<recob::Track>   ftracklist;
    art::PtrVector<recob::Vertex>  feplistStrong;

  }; // class AggregateVertex

}  // Namespace vertex

#endif // AGGREGATEVERTEX_H
