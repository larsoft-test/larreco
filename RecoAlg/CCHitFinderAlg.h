////////////////////////////////////////////////////////////////////////
// ClusterCrawlerAlg.h
//
// ClusterCrawlerAlg class
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef CCHITFINDERALG_H
#define CCHITFINDERALG_H


// Insert code for studying hit fitting. This study is required to determine
// the best values of the fcl inputs - ChiNorms, MinSigInd, MinRMSInd, etc.
// This study only needs to be done once for a detector configuration using
// MC and real data. This study should be performed on a single event that 
// contains one or more shallow angle tracks with minimal activity elsewhere
#define STUDYHITS

// Insert code for printing out hit finding/fitting
// #define PRINTHITS

#include "TMath.h"

#include <vector>

#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 

#include "Geometry/Geometry.h"
#include "RecoBase/Hit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"

//namespace recob { class Hit; }

namespace cluster {

  class CCHitFinderAlg {
  
  public:
    
    struct CCHit {
      float Charge;
      float ChargeErr;
      float Amplitude;
      float AmplitudeErr;
      float Time;
      float TimeErr;
      float RMS;
      float RMSErr;
      float ChiDOF;
      art::Ptr<recob::Wire> Wire;
      unsigned short WireNum;
      unsigned short numHits;
      unsigned short LoHitID;
    };
    std::vector< CCHit > allhits;
    
    // struct for passing hit fitting cuts to ClusterCrawler
    struct HitCuts {
      float MinSigInd;
      float MinSigCol;
      float MinRMSInd;
      float MinRMSCol;
      float ChiSplit;
      std::vector<float> ChiNorms;
      std::vector<float> TimeOffsets;
      std::vector<float> ChgNorms;
    };
    HitCuts hitcuts;

    CCHitFinderAlg(fhicl::ParameterSet const& pset);
    virtual ~CCHitFinderAlg();

    void reconfigure(fhicl::ParameterSet const& pset);

    void RunCCHitFinder(art::Event & evt);
    
  private:
    
    std::string     fCalDataModuleLabel;
    float fMinSigInd;     ///<Induction signal height threshold 
    float fMinSigCol;     ///<Collection signal height threshold 
    float fMinRMSInd;      ///<Initial rms for induction fit
    float fMinRMSCol;      ///<Initial rms for collection fit
    unsigned short fMaxBumps; // make a crude hit if > MaxBumps are found in the RAT
    unsigned short fMaxXtraHits; // max num of hits in Region Above Threshold
    float fChiSplit;      ///<Estimated noise error on the Signal
    float ChgNorm;     // Area norm for the wire we are working on

    std::vector<float> fChiNorms;
    std::vector<float> fTimeOffsets;
    std::vector<float> fChgNorms;

    uint32_t theChannel;
    unsigned short theWireNum;
    unsigned short thePlane;
    float minRMS;
    float minSig;
    float chinorm;
    float timeoff;
    const float Sqrt2Pi = 2.5066;
    const float SqrtPi  = 1.7725;

#ifdef STUDYHITS
    std::vector<int> bumpCnt;
    std::vector<float> bumpChi;
    std::vector<float> bumpRMS;
    std::vector<int> hitCnt;
    std::vector<float> hitRMS;
#endif

#ifdef PRINTHITS
    bool prt;
#endif
    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;

    // fit n Gaussians possibly with bounds setting (parmin, parmax)
    void FitNG(unsigned short nGaus, unsigned short npt, float *ticks,
       float *signl);
    // parameters, errors, lower limit, upper limits for FitNG
    std::vector<double> par;
    std::vector<double> parerr;
    std::vector<double> parmin;
    std::vector<double> parmax;
    float chidof;
    std::vector<unsigned short> bumps;
    
    // make a cruddy hit if fitting fails
    void MakeCrudeHit(unsigned short npt, float *ticks, float *signl);
    // store the hits
    void StoreHits(unsigned short TimeOffset, art::Ptr<recob::Wire>& theWire);

  }; // class CCHitFinderAlg

} // cchit

#endif // ifndef CCHITFINDERALG_H
