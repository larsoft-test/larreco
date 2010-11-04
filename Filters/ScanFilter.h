////////////////////////////////////////////////////////////////////////
//
// ScanFilter class:
// Tells the downstream rreconstruction to not process events that
// do not pass certain hand scan criteria
//
// joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef SCANFILTER_H
#define SCANFILTER_H

#include "FWCore/Framework/interface/EDFilter.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TH2D.h"
namespace filt {

  class ScanFilter : public edm::EDFilter  {
    
  public:
    
    explicit ScanFilter(edm::ParameterSet const& ); 
    virtual ~ScanFilter();
         
    
    bool filter(edm::Event& evt, edm::EventSetup const&);
    void beginJob(const edm::EventSetup&);
   
  private: 
 
    std::string fScanModuleLabel;
    int fNeutrino_req, fNumShowers_req, fNumTracks_req;
  
  protected: 
    
  }; // class ScanFilter

}

#endif // SCANFILTER_H
