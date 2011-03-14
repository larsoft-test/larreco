#ifndef FFTHITFINDER_H
#define FFTHITFINDER_H

#include <vector>
#include <string>

#include "art/Framework/Core/EDProducer.h" 

///localizations of energy depositions
namespace hit {
   
  class FFTHitFinder : public art::EDProducer {
    
  public:
    
    explicit FFTHitFinder(fhicl::ParameterSet const& ); 
    virtual ~FFTHitFinder();
         
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 

  private:
        
    std::string     fCalDataModuleLabel;
    double          fMinSigInd;     ///>Induction signal height threshold 
    double          fMinSigCol;     ///>Collection signal height threshold 
    double          fIndWidth;      ///>Initial width for induction fit
    double          fColWidth;      ///>Initial width for collection fit
    double          fPOffset;       ///>Time delay between planes
    double          fOOffset;       ///>Distance(cm) from induction plane to origin
    int             fMaxMultiHit;   ///>maximum hits for multi fit   
  protected: 
    
  
  }; // class FFTHitFinder


}

#endif // FFTHITFINDER_H
