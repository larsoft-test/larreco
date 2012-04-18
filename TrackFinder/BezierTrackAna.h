#ifndef BEZIERTRACKANA_H
#define BEZIERTRACKANA_H

//
// Name: BezierTrackAna.h
//
// Purpose: Header file for module BezierTrackAna.
//
// Created: April 2012, bjpjones@mit.edu
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "TH1F.h"
#include "TH2D.h"
class TTree;

namespace trkf {

  class BezierTrackAna : public art::EDAnalyzer
  {
  public:
 
    // Constructors, destructor

    explicit BezierTrackAna(fhicl::ParameterSet const& pset);
    virtual ~BezierTrackAna();

    // Book histograms.

    void bookHistograms(bool mc);

    // Overrides.

    void analyze(const art::Event& evt);

  private:

    // Fcl Attributes.

    std::string fBezierTrackModuleLabel;
    bool fBooked;

    // Histograms
    
    //TH1D * fRMSCurvatures;
    //TH1D * fLengths;
    
    TH1D * fhHitDistance;
    TH1D * fhHitS;
    TH1D * fhdQdxU;
    TH1D * fhdQdxV;
    TH1D * fhdQdxW;
    TH2D * fhdQdxVW;
    TH1D * fhCurv;
    
    
   
    //Entries for TTree
    TTree*  fTree;
    Float_t fLength;
    Float_t fRMSCurvature;
    Int_t   fNumEvent;
    Int_t   fNHitsU;
    Int_t   fNHitsV;
    Int_t   fNHitsW;
    Int_t   fNHits;
    Float_t fAverageS;
    Float_t fAverageDistance;
    
    Float_t fdQdxU;
    Float_t fdQdxV;
    Float_t fdQdxW;
    
    Float_t fChargeU;
    Float_t fChargeV;
    Float_t fChargeW;
    
    Int_t fSegments;


  };
}

#endif 
