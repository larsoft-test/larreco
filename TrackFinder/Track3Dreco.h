#ifndef TRACK3DRECO_H
#define TRACK3DRECO_H

#include "FWCore/Framework/interface/EDProducer.h" //

#include <vector>
#include <string>

namespace trkf {
   
  class Track3Dreco : public edm::EDProducer {
    
  public:
    
    explicit Track3Dreco(edm::ParameterSet const& pset);
    ~Track3Dreco();
    
    //////////////////////////////////////////////////////////
    void produce(edm::Event& evt, edm::EventSetup const&); 
    void beginJob(const edm::EventSetup&);
    void endJob();

  private:
        
    int             ftmatch; // tolerance for time matching (in time samples) 
    double          fchi2dof;// tolerance for chi2/dof of cluster fit to function
    std::string     fClusterModuleLabel;// label for input cluster collection
  protected: 
    
  
  }; // class Track3Dreco

}

#endif // TRACK3DRECO_H
