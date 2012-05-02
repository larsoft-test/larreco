////////////////////////////////////////////////////////////////////////
//
// GausHitFinder class designed to analyze signal on a wire in the TPC
//
// jaasaadi@syr.edu
//
// Note: This has been based (stolen) from the FFTHitFinderAna thus
// there is still some unneeded hold overs that will get cleaned up later
////////////////////////////////////////////////////////////////////////
#ifndef GAUSHITFINDERANA_H
#define GAUSHITFINDERANA_H


#include "RecoBase/recobase.h"
#include "Utilities/LArFFT.h"

#include "TComplex.h"
#include "TString.h"
#include "TGraph.h"
#include "TH2.h"
#include "TTree.h"

#include "art/Framework/Core/EDAnalyzer.h"

#include <vector>
#include <string>


namespace geo { class Geometry; }

///Detector simulation of raw signals on wires
namespace hit {

  /// Base class for creation of raw signals on wires. 
  class GausHitFinderAna : public art::EDAnalyzer {
    
  public:
        
    explicit GausHitFinderAna(fhicl::ParameterSet const& pset); 
    virtual ~GausHitFinderAna();
    
    /// read/write access to event
    void analyze (const art::Event& evt);
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    std::string            fGausHitFinderModuleLabel;
    std::string            fLArG4ModuleLabel;
    
      TTree* fHTree;
      Int_t fRun;
      Int_t fEvt;
      Int_t fnhits; //<---Number of Hits in the Event
      Int_t fnOnePulseHits; //<---Number of One Pulse hit events
      Int_t fmulitPulseHits; //<---Number of Multi pulse hit events
      
      Int_t fSingleHit; //<<---Set a indicator to know if this is a single pulse hit or multihit
      Int_t fMultiHit;  //<<---Set a indicator to know if this is a single pulse hit or multihit
      
      Float_t fWiren1; //<---Wire number of hit with multiplicity of 1
      Float_t fgoodoffitn1; //<---Goodness of hit with mulitplicity of 1
      Float_t fChargen1; //<---Charge of hit with multiplicity of 1
      Float_t fSigmaChargen1; //<---Uncertainty of charge of hit with multiplicity of 1
      Float_t fWidthn1; //<---Width of the Hit (Endtime - Peaktime) with multiplicity of 1
      Float_t fPeakn1; //<---Peak time of hit with multiplicity of 1
      Float_t fPeakUncertn1; //<---Uncertainty in the peak position of the hit with multiplicity of 1
      Float_t fStartTimen1; //<---Start Position of the hit with multiplicity of 1
      Float_t fStartTimeUncertn1; //<---Start Time Uncertainty of the hit with multiplicity of 1
      Float_t fEndTimen1; //<---End Position of the hit with multiplicity of 1
      Float_t fEndTimeUncertn1; //<---End Time Uncertainty of the hit with multiplicity of 1
      
      
      Float_t fWirenGT1; //<---Wire number of hit with multiplicity of 1
      Float_t fWidthnGT1; //<---Width of the Hit (Endtime - Peaktime) with multiplicity > 1
      Float_t fPeaknGT1; //<---Peak time of hit with multiplicity > 1
      Float_t fPeakUncertnGT1; //<---Uncertainty in the peak position of the hit with multiplicity > 1
      Float_t fSigmaChargenGT1; //<---Uncertainty of charge of hit with multiplicity > 1
      Float_t fgoodoffitnGT1; //<---Goodness of hit with mulitplicity > 1
      Float_t fChargenGT1; //<---Charge of hit with multiplicity > 1
      Float_t fStartTimenGT1; //<---Start Position of the hit with multiplicity > 1
      Float_t fStartTimeUncertnGT1; //<---Start Time Uncertainty of the hit with multiplicity > 1
      Float_t fEndTimenGT1; //<---End Position of the hit with multiplicity > 1
      Float_t fEndTimeUncertnGT1; //<---End Time Uncertainty of the hit with multiplicity > 1
      
      Int_t fNp0;
      Int_t fNp1;
      Int_t fNp2;
      Int_t fN3p0;
      Int_t fN3p1;
      Int_t fN3p2;
      Float_t* fPeakTime0;
      Float_t* fPeakTime1;
      Float_t* fPeakTime2;
      Int_t* fWirep0;
      Int_t* fWirep1;
      Int_t* fWirep2;
      Float_t* fChgp0;
      Float_t* fChgp1;
      Float_t* fChgp2;
      Float_t* fXYZp0;
      Float_t* fXYZp1;
      Float_t* fXYZp2;

      Int_t*  fMCPdg0;
      Int_t*  fMCTId0;
      Float_t*  fMCE0;
      Int_t*  fMCPdg1;
      Int_t*  fMCTId1;
      Float_t*  fMCE1;
      Int_t*  fMCPdg2;
      Int_t*  fMCTId2;
      Float_t*  fMCE2;

  }; // class GausHitFinderAna

}

#endif // GAUSHITFINDERANA_H
