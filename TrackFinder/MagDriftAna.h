////////////////////////////////////////////////////////////////////////
// $Id$
//
// Measure the charge drift induced by a magnetic field
//
// dmckee@phys.ksu.edu
//
//
////////////////////////////////////////////////////////////////////////
#ifndef MAGDRIFTANA_H
#define MAGDRIFTANA_H


#include "TComplex.h"
#include "TString.h"
#include "TGraph.h"
class TH1D;
class TH2D;

#include "art/Framework/Core/EDAnalyzer.h"

#include <vector>
#include <string>


namespace geo { class Geometry; }

///Detector simulation of raw signals on wires
namespace hit {

  /// Base class for creation of raw signals on wires. 
  class MagDriftAna : public art::EDAnalyzer {
    
  public:
        
    explicit MagDriftAna(fhicl::ParameterSet const& pset); 
    virtual ~MagDriftAna();
    
    /// read/write access to event
    void analyze (const art::Event& evt);
    void beginJob(){};
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);

    // intilize the histograms
    //
    // Can't be done in Begin job because I want to use LArProperties
    // which used the database, so I test and run on each
    // event. Wasteful and silly, but at least it *works*.
    void ensureHists();

  private:

    std::string            fFFTHitFinderModuleLabel;
    std::string            fTrackFinderModuleLabel;
    std::string            fLArG4ModuleLabel;

    // Flag for initialization done, because we set up histograms the
    // first time through beginRun() so that we can use the
    // database...
    bool initDone;

    // Drift properties
    double fDirCosY;
    double fDirCosZ;

    TH1D * fChargeXpos;  // << position of the MC Truth charge deposition 
    TH1D * fChargeYpos;
    TH1D * fChargeZpos;
    TH1D * fHitZpos;     // << Z position of the recorded hit (from the
			 //    z-sensitive wire)

    TH1D * fDriftDeltaZ; // << Difference in MC charge Z and recorded hit Z
    TH1D * fDeltaZoverX; // << Delta Z as a function of drift distance
    TH2D * fDeltaZvsX;

    // Same as above, but only for long drift distances (greater than
    // 4/5 of the detector)
    TH1D * fDriftDeltaZAway; // << Difference in MC charge Z and recorded hit Z
    TH1D * fDeltaZoverXAway; // << Delta Z as a function of drift distance

  }; // class MagdriftAna

}

#endif // MAGDRIFTANA_H
