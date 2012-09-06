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
#include "Simulation/sim.h"
#include "SimulationBase/simbase.h"

#include "TComplex.h"
#include "TString.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"

#include "art/Framework/Core/EDAnalyzer.h"

#include <vector>
#include <string>
using namespace std;

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
    std::string		   fHitCheaterModuleLabel;
    std::vector<const sim::SimChannel*>    fSimChannels;           ///< all the SimChannels for the event
      
      TH1F* fRun;
      TH1F* fEvt;
      TH1F* fWireNumbMulti1; //<---Wire number of hit with multiplicity of 1
      TH1F* fFitGoodnessMulti1; //<---Goodness of hit with mulitplicity of 1
      TH1F* fChargeMulti1; //<---Charge of hit with multiplicity of 1
      
      TH1F* fTruthPeakPosition;//<---Peak time of hit from backtracker with multiplicity of 1
      TH1F* fTruthPeakPositionPlane0Multi1;
      TH1F* fTruthPeakPositionPlane1Multi1;
      TH1F* fTruthPeakPositionPlane2Multi1;
      
      TH1F* fRecoPeakPositionMulti1;//<---Peak time of hit with multiplicity of 1
      TH1F* fRecoPeakPositionPlane0Multi1;//<---Peak time of hit with multiplicity of 1 in plane 0
      TH1F* fRecoPeakPositionPlane1Multi1;//<---Peak time of hit with multiplicity of 1 in plane 1
      TH1F* fRecoPeakPositionPlane2Multi1;//<---Peak time of hit with multiplicity of 1 in plane 2
      
      TH1F* fRecoPeakPositionUncertMulti1;//<---Peak time of hit with multiplicity of 1
      TH1F* fRecoPeakPositionUncertPlane0Multi1;//<---Peak time of hit with multiplicity of 1 in plane 0
      TH1F* fRecoPeakPositionUncertPlane1Multi1;//<---Peak time of hit with multiplicity of 1 in plane 1
      TH1F* fRecoPeakPositionUncertPlane2Multi1;//<---Peak time of hit with multiplicity of 1 in plane 2
      
      TH1F* fRecoPeakPositionUncertMultiGT1;//<---Peak time of hit with multiplicity of > 1
      TH1F* fRecoPeakPositionUncertPlane0MultiGT1;//<---Peak time of hit with multiplicity of > 1 in plane 0
      TH1F* fRecoPeakPositionUncertPlane1MultiGT1;//<---Peak time of hit with multiplicity of > 1 in plane 1
      TH1F* fRecoPeakPositionUncertPlane2MultiGT1;//<---Peak time of hit with multiplicity of > 1 in plane 2
      
      TH1F* fHitResidualAll;
      TH1F* fHitResidualMulti1;
      TH1F* fHitResidualPlane0Multi1;
      TH1F* fHitResidualPlane1Multi1;
      TH1F* fHitResidualPlane2Multi1;
      
      TTree* fHTree;
      //Int_t fRun;
      //Int_t fEvt;
      Int_t fnhits; //<---Number of Hits in the Event
      Int_t fnOnePulseHits; //<---Number of One Pulse hit events
      Int_t fmulitPulseHits; //<---Number of Multi pulse hit events
      
      
      Int_t fSingleHit; //<<---Set a indicator to know if this is a single pulse hit or multihit
      Int_t fMultiHit;  //<<---Set a indicator to know if this is a single pulse hit or multihit
      
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
      
      
      Int_t fWirenGT1; //<---Wire number of hit with multiplicity of 1
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
/*      
/////////////////////////////////////////////////////////////////////////////////////
// 
// Declare and fill 1D histo with Title
//
void FillHisto( double var, const string& hname, 
				int nbins, float xmin, float xmax ) 
{
  // If histo not already booked, book it
  if ( !gDirectory->FindObject( hname.c_str() ) ) {
    std::vector<TH1D*> vectHist ;
    vectHist.push_back( new TH1D( hname.c_str(), "", nbins, xmin, xmax ) )  ;
  }

  // Fill histo
  ((TH1D*)gDirectory->Get(hname.c_str()))->Fill( var ) ;

}

/////////////////////////////////////////////////////////////////////////////////////
// 
// Book and fill 1D histo with Title
//
void FillHisto( double var, const string& hname, const string& htit, 
				int nbins, float xmin, float xmax ) 
{
  // If histo not already booked, book it
  if ( !gDirectory->FindObject( hname.c_str() ) ) {
    std::vector<TH1D*> vectHist ;
    vectHist.push_back( new TH1D( hname.c_str(), "", nbins, xmin, xmax ) )  ;
    vectHist[vectHist.size()-1]->SetXTitle( htit.c_str() ) ;
  }

  // Fill histo
  ((TH1D*)gDirectory->Get(hname.c_str()))->Fill( var ) ;

}


/////////////////////////////////////////////////////////////////////////////////////
//
// Book and fill 2D histo
//
void FillHisto( double var1, double var2, const string& hname, 
				int nbinsx, float xmin, float xmax, 
				int nbinsy, float ymin, float ymax ) 
{
  // If histo not already booked, book it
  if ( !gDirectory->FindObject( hname.c_str() ) ) {
    std::vector<TH2D*> vectHist ;
    vectHist.push_back( new TH2D( hname.c_str(), "", nbinsx, xmin, xmax,
                                  nbinsy, ymin, ymax ) )  ;
  }

  // Fill histo
  ((TH2D*)gDirectory->Get(hname.c_str()))->Fill( var1, var2 ) ;

}

/////////////////////////////////////////////////////////////////////////////////////
//
// Book and fill 2D histo with Title
//
void FillHisto( double var1, double var2, const string& hname, const string& htit,
				int nbinsx, float xmin, float xmax, 
				int nbinsy, float ymin, float ymax ) 
{
  // If histo not already booked, book it
  if ( !gDirectory->FindObject( hname.c_str() ) ) {
    std::vector<TH2D*> vectHist ;
    vectHist.push_back( new TH2D( hname.c_str(), "", nbinsx, xmin, xmax,
                                  nbinsy, ymin, ymax ) )  ;
    vectHist[vectHist.size()-1]->SetXTitle( htit.c_str() ) ;
  }

  // Fill histo
  ((TH2D*)gDirectory->Get(hname.c_str()))->Fill( var1, var2 ) ;

}


///////////////////////////////////////////////////////////////
//
// Book and fill 1D profile histo
//
void FillProf( double var1, double var2, const string& hname, 
			       int nbins, float xmin, float xmax ,float ymin, float ymax ) 
{
  // If histo not already booked, book it
  if ( !gDirectory->FindObject( hname.c_str() ) ) {
    std::vector<TProfile*> vectHist ;
    vectHist.push_back( new TProfile( hname.c_str(), "", nbins, xmin, xmax, ymin, ymax )
                        )  ;
  }

  // Fill histo
  ((TProfile*)gDirectory->Get(hname.c_str()))->Fill( var1, var2, 1. ) ;

}


///////////////////////////////////////////////////////////////
//
// Book and fill 2D profile histo
//
void FillProf( double var1, double var2, const string& hname, const string& htit,
			       int nbins, float xmin, float xmax ,float ymin, float ymax ) 
{
  // If histo not already booked, book it
  if ( !gDirectory->FindObject( hname.c_str() ) ) {
    std::vector<TProfile*> vectHist ;
    vectHist.push_back( new TProfile( hname.c_str(), "", nbins, xmin, xmax, ymin, ymax )
                        )  ;
    vectHist[vectHist.size()-1]->SetXTitle( htit.c_str() ) ;
  }

  // Fill histo
  ((TProfile*)gDirectory->Get(hname.c_str()))->Fill( var1, var2, 1. ) ;

}*/
  }; // class GausHitFinderAna

}

#endif // GAUSHITFINDERANA_H
