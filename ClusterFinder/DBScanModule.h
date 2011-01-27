/////////////////////////////////////////////////////////////////
//  DBSCANfinder.h
//  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////

#include <vector>
#include <cmath>
#include <iostream>
#include "FWCore/Framework/interface/EDProducer.h"

class TH1F;

namespace cluster{
   
  //--------------------------------------------------------------- 
  class DBScanModule : public edm::EDProducer
  {
  public:
    explicit DBScanModule(edm::ParameterSet const& pset); 
    ~DBScanModule();
    void produce(edm::Event& evt, edm::EventSetup const&);
    void beginJob(edm::EventSetup const&);
    
  private:
       
    TH1F *fhitwidth;
    TH1F *fhitwidth_ind_test;  
    TH1F *fhitwidth_coll_test;  
    
  protected:
    std::string fhitsModuleLabel;
    
  };

}
