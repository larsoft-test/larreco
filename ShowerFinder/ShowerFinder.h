#ifndef SHOWERFINDER_H
#define SHOWERFINDER_H

#include <vector>
#include <string>

#include "FWCore/Framework/interface/EDProducer.h" 


namespace shwf {
   
  class ShowerFinder : public edm::EDProducer {
    
  public:
    
    explicit ShowerFinder(edm::ParameterSet const&); 
    virtual ~ShowerFinder();
         
    void produce(edm::Event& evt, edm::EventSetup const&); 
    void beginJob(const edm::EventSetup&); 
    void endJob(); 
 

  private:

    std::string     fvertices;  
    std::string     fclusters;  
    std::string     fhoughlines;  
    std::string     fhits;  
    double     fRcone;  
    double     fLcone;  

  protected: 
    
    
  }; // class ShowerFinder
  
  
}

#endif // SHOWERFINDER_H
