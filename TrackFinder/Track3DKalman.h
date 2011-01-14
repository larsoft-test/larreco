#ifndef TRACK3DRECO_H
#define TRACK3DRECO_H

#include "FWCore/Framework/interface/EDProducer.h"
#include <TTree.h>
#include <TMatrixT.h>
#include <TRandom3.h>

#include <vector>
#include <string>

//#include "RecoBase/SpacePoint.h"

namespace trkf {


   
  class Track3DKalman : public edm::EDProducer {
    
  public:
    
    explicit Track3DKalman(edm::ParameterSet const& pset);
    ~Track3DKalman();
    
    //////////////////////////////////////////////////////////
    void produce(edm::Event& evt, edm::EventSetup const&); 
    void beginJob(const edm::EventSetup&);
    void endJob();

  private:
        
    int             ftmatch; // tolerance for time matching (in time samples) 
    double          fchi2dof;// tolerance for chi2/dof of cluster fit to function
    std::string     fHitModuleLabel;// label for input cluster collection
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


  protected: 
    
  
  }; // class Track3DKalman

}

#endif // TRACK3DRECO_H
