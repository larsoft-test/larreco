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
    
    TH1D * fRMSCurvatures;
    TH1D * fLengths;
    
   
    //Entries for TTree
    TTree*  fTree;
    Float_t fLength;
    Float_t fRMSCurvature;
    Int_t   fNumEvent;

  };
}

#endif 
