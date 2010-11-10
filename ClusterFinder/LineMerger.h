#ifndef LINEMERGER_H
#define LINEMERGER_H

#include "FWCore/Framework/interface/EDProducer.h"
#include <string>

namespace cluster {
   
  class LineMerger : public edm::EDProducer {
    
  public:
    
    explicit LineMerger(edm::ParameterSet const& pset); 
    ~LineMerger();
    
    void produce(edm::Event& evt, edm::EventSetup const&);
    void beginJob(edm::EventSetup const&);
    
  private:
        
    std::string     fClusterModuleLabel;
    double          fSlope; // tolerance for matching Slopes 
    double          fIntercept; // tolerance for matching Intercepts (in time samples) 
   
    bool SlopeCompatibility(double slope1,double slope2);
    bool InterceptCompatibility(double intercept1,double intercept2);
    
  protected: 
    
  }; // class LineMerger

}

#endif // LINEMERGER_H
