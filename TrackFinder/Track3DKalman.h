////////////////////////////////////////////////////////////////////////
//
//   \file Track3DKalman.h
//
//   soderber@fnal.gov
//   kinga.partyka@yale.edu
//   joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef TRACK3DRECO_H
#define TRACK3DRECO_H

#include "art/Framework/Core/EDProducer.h"
#include <TTree.h>
#include <TMatrixT.h>
#include <TRandom3.h>

#include "Genfit/GFAbsTrackRep.h"

#include <vector>
#include <string>

//#include "RecoBase/SpacePoint.h"

namespace trkf {

  class Track3DKalman : public art::EDProducer {
    
  public:
    
    explicit Track3DKalman(fhicl::ParameterSet const& pset);
    ~Track3DKalman();
    
    //////////////////////////////////////////////////////////
    void produce(art::Event& evt); 
    void beginJob();
    void endJob();

  private:
        
    int             ftmatch; // tolerance for time matching (in time samples) 
    double          fchi2dof;// tolerance for chi2/dof of cluster fit to function
    std::string     fTrackModuleLabel;// label for input collection
    std::string     fGenieGenModuleLabel;// label for input MC single particle generator
    TRandom3*              fRandom;           //< random number generator 
    bool fGenfPRINT;
      
    TFile *fileGENFIT;
    TTree *tree;

    TMatrixT<Double_t> *stMCT;
    TMatrixT<Double_t> *covMCT;
    TMatrixT<Double_t> *stREC;
    TMatrixT<Double_t> *covREC;
    Float_t chi2;
    Float_t chi2ndf;
    Double_t *pREC;
    Double_t *pMCT;
    int nfail;
    int ndf;

    genf::GFAbsTrackRep *repMC;
    genf::GFAbsTrackRep *rep;

  protected: 
    
  
  }; // class Track3DKalman

}

#endif // TRACK3DRECO_H
