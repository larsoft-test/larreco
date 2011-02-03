////////////////////////////////////////////////////////////////////////
// $Id: HoughLineFinderAna.cxx,v 1.36 2010/09/15  bpage Exp $
//
// HoughLineFinderAna class
//
// josh
//
////////////////////////////////////////////////////////////////////////
#ifndef HOUGHLINEMODULE_H
#define HOUGHLINEMODULE_H

#include "TMath.h"
#include <vector>
#include <string>

#include "art/Framework/Core/EDProducer.h"

class TH1F;
class TTree;

namespace cluster {
   
  class HoughLineModule : public art::EDProducer {
    
  public:
    
    explicit HoughLineModule(fhicl::ParameterSet const& pset); 
    virtual ~HoughLineModule();
         
    void produce(art::Event& evt);
     
    
  private:
    std::string fDBScanModuleLabel;    
    
      
  
  };
  
  
}



#endif // HOUGHLineMODULE_H
