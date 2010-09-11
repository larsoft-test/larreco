#ifndef FFTHITFINDER_H
#define FFTHITFINDER_H

#include <vector>
#include <string>

#include "FWCore/Framework/interface/EDProducer.h" 


#include "CalData/inc/CalWire.h"
#include "RecoBase/inc/recobase.h"
#include "Utilities/inc/LArFFT.h"

namespace hit {
   
  class FFTHitFinder : public edm::EDProducer {
    
  public:
    
    explicit FFTHitFinder(edm::ParameterSet const& ); 
    virtual ~FFTHitFinder();
         
    void produce(edm::Event& evt, edm::EventSetup const&); 
    void beginJob(const edm::EventSetup&); 
    void endJob();                 

  private:
        
    bool            fSpacer;
    std::string     fCalDataModuleLabel;
    double          fMinSigInd;     //Induction signal height threshold 
    double          fMinSigCol;     //Collection signal height threshold 
    double          fIndWidth;      //Initial width for induction fit
    double          fColWidth;     //Initial width for collection fit
    double          fDrift;    //Drift Velocity
    double          fPOffset;  //Time delay between planes
    double          fOOffset;  //Distance(cm) from induction plane to origin
    int             fMaxMultiHit;   //maximum hits for multi fit   
  protected: 
    
  
  }; // class FFTHitFinder


}

#endif // FFTHITFINDER_H
