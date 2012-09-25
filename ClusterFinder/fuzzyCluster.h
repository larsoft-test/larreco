/////////////////////////////////////////////////////////////////
//  \file fuzzyCluster.h
//  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////

#include <vector>
#include <cmath>
#include <iostream>
#include "art/Framework/Core/EDProducer.h"
#include "ClusterFinder/fuzzyClusterAlg.h"

class TH1F;

namespace cluster{
   
  //--------------------------------------------------------------- 
  class fuzzyCluster : public art::EDProducer
  {
  public:
    explicit fuzzyCluster(fhicl::ParameterSet const& pset); 
    ~fuzzyCluster();
    void produce(art::Event& evt);
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);
    
  private:
       
    TH1F *fhitwidth;
    TH1F *fhitwidth_ind_test;  
    TH1F *fhitwidth_coll_test;  
  
    std::string fhitsModuleLabel;
   
    fuzzyClusterAlg ffuzzyCluster; ///< object that implements the fuzzy cluster algorithm
  };

}
