#ifndef VERTEXMATCH_H
#define VERTEXMATCH_H

#include "TMath.h"
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include "FWCore/Framework/interface/EDProducer.h"

namespace vertex {
   
  class VertexMatch : public edm::EDProducer {
    
  public:
    
    explicit VertexMatch(edm::ParameterSet const& pset); 
    ~VertexMatch();
    void produce(edm::Event& evt, edm::EventSetup const&);
    void beginJob(edm::EventSetup const&);
    
  private:
    std::string fVertexModuleLabel; 
    std::string fHoughClusterLabel;
    double fMaxDistance;
  protected:

    
    
  };


    
}

#endif // VERTEXMATCH_H
