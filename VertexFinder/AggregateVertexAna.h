
#ifndef AGGREGATEVTXANA_H
#define AGGREGATEVTXANA_H

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
#include "art/Framework/Core/EDAnalyzer.h"

// LArSoft includes
#include "RecoBase/recobase.h"

#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include <vector>
#include <string>


namespace vertex {


  class AggregateVertexAna : art::EDAnalyzer 
  {

  public:

    explicit AggregateVertexAna(fhicl::ParameterSet const& pset);
    ~AggregateVertexAna();

    void analyze (const art::Event& evt);
    void beginJob(); 

  private:

    TH1F* HnTrksVtx;
    TH1F* HnVtxes;
    TH1F* HVtxSep;
    TH2F* HVtxRZ;

    std::string fDBScanModuleLabel;
    std::string fHoughModuleLabel;
    std::string fHitModuleLabel;
    std::string fTrack3DModuleLabel;
    std::string fEndPointModuleLabel;
    std::string fVertexModuleLabel;

    art::PtrVector<recob::Hit>        fhitlist;
    art::PtrVector<recob::EndPoint2D> feplist;
    art::PtrVector<recob::Track>      ftracklist;
    art::PtrVector<recob::Vertex>     fVertexlist;


  }; // class AggregateVertexAna

}  // Namespace aggr

#endif // AGGREGATEVTXANA_H
