////////////////////////////////////////////////////////////////////////
//
//   \file Track3DKalmanSPS.h
//
//   soderber@fnal.gov
//   kinga.partyka@yale.edu
//   joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////
#ifndef TRACK3DKSPS_H
#define TRACK3DKSPS_H

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include <TTree.h>
#include <TMatrixT.h>

#include "Genfit/GFAbsTrackRep.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

#include <vector>
#include <string>

//#include "RecoBase/SpacePoint.h"

namespace trkf {

  class Track3DKalmanSPS : public art::EDProducer {
    
  public:
    
    explicit Track3DKalmanSPS(fhicl::ParameterSet const& pset);
    virtual ~Track3DKalmanSPS();
    
    //////////////////////////////////////////////////////////
    void produce(art::Event& evt); 
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);
    double energyLossBetheBloch(const double& mass,
				const double p
				);

  private:
        
    std::string     fClusterModuleLabel;// label for input collection
    std::string     fProngModuleLabel;// label for input collection
    std::string     fGenieGenModuleLabel;// label for input MC single particle generator
    std::string     fG4ModuleLabel;// label for input MC single particle generator
    bool fGenfPRINT;
      
    TFile *fileGENFIT;
    TTree *tree;

    TMatrixT<Double_t> *stMCT;
    TMatrixT<Double_t> *covMCT;
    TMatrixT<Double_t> *stREC;
    TMatrixT<Double_t> *covREC;
    Float_t chi2;
    Float_t chi2ndf;
    int     fcont;

    Float_t *fpRECt3D;    
    Float_t *fpRECLE;
    Float_t *fpRECL;
    Float_t *fpREC;
    Float_t *fpMCMom;
    Float_t *fpMCPos;
    Float_t *fState0;
    Float_t *fCov0;
    int nfail;
    int ndf;
    unsigned int evtt;
    unsigned int nTrks;
    unsigned int fptsNo;
    Float_t *fshx;
    Float_t *fshy;
    Float_t *fshz;
    Float_t *feshx;
    Float_t *feshy;
    Float_t *feshz;
    Float_t *feshyz;
    Float_t *fupdate;
    Float_t *fth;
    Float_t *feth;
    Float_t *fedudw;
    Float_t *fedvdw;
    Float_t *feu;
    Float_t *fev;
    Float_t *fsep;
    unsigned int fDimSize; // if necessary will get this from pset in constructor.
    Float_t *fPCmeans;
    Float_t *fPCevals;
    Float_t *fPCsigmas;
    Float_t *fPC1;
    Float_t *fPC2;
    Float_t *fPC3;

    std::vector<double> fPosErr;
    std::vector<double> fMomErr;
    std::vector<double> fMomStart;
    double fPerpLim;
    bool fDoFit;
    int fNumIt;
    double fErrScaleS;
    double fErrScaleM;
    int fDecimate;
    double fMaxUpdate;
    int fDecimateU;
    double fMaxUpdateU;
    double fMomLow;
    double fMomHigh;
    int fPdg;
    
    genf::GFAbsTrackRep *repMC;
    genf::GFAbsTrackRep *rep;

  protected: 
    
  
  }; // class Track3DKalmanSPS

}

#endif // TRACK3DKSPS_H
