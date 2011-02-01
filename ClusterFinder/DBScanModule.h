/////////////////////////////////////////////////////////////////
//  DBSCANfinder.h
//  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////

#include <vector>
#include <cmath>
#include <iostream>
#include "art/Framework/Core/EDProducer.h"

class TH1F;

namespace cluster{
   
  //--------------------------------------------------------------- 
  class DBScanModule : public art::EDProducer
  {
  public:
    explicit DBScanModule(fhicl::ParameterSet const& pset); 
    ~DBScanModule();
    void produce(art::Event& evt);
    void beginJob();
    
  private:
       
    TH1F *fhitwidth;
    TH1F *fhitwidth_ind_test;  
    TH1F *fhitwidth_coll_test;  
    
    std::string fhitsModuleLabel;
    
  };

}
