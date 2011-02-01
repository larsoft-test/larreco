////////////////////////////////////////////////////////////////////////
// $Id: HoughLineFinderAna.cxx,v 1.36 2010/09/15  bpage Exp $
//
// HoughLineFinderAna class
//
// josh
//
////////////////////////////////////////////////////////////////////////
#ifndef LINEMERGER_H
#define LINEMERGER_H

#include "art/Framework/Core/EDProducer.h"
#include <string>

namespace cluster {
   
  class LineMerger : public art::EDProducer {
    
  public:
    
    explicit LineMerger(fhicl::ParameterSet const& pset); 
    ~LineMerger();
    
    void produce(art::Event& evt);
    void beginJob();
    
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
