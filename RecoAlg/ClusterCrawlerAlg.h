////////////////////////////////////////////////////////////////////////
// ClusterCrawlerAlg.h
//
// ClusterCrawlerAlg class
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef CLUSTERCRAWLERALG_H
#define CLUSTERCRAWLERALG_H

#include "TMath.h"

#include <vector>
#include "TLinearFitter.h"
#include "TF1.h"

#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 

#include "Geometry/Geometry.h"
#include "RecoBase/Hit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

namespace recob { class Hit; }

namespace cluster {

    // structure of temporary clusters

    struct ClusterStore {
      int ID;         // Cluster ID = 1000 * pass + number
      int StopCode;   // code for the reason for stopping cluster tracking
      float BeginSlp; // beginning slope (= DS end = high wire number)
      int   BeginWir; // begin wire
      float BeginTim; // begin time
      float BeginChg; // beginning average charge
      float EndSlp;   // end slope (= US end = low  wire number)
      int   EndWir;   // end wire
      float EndTim;   // ending time
      float EndChg;   // ending average charge
      std::vector<int> tclhits; // hits on the cluster
    };
    std::vector< ClusterStore > tcl;


  class ClusterCrawlerAlg {
    public:

    ClusterCrawlerAlg(fhicl::ParameterSet const& pset);
    virtual ~ClusterCrawlerAlg();

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void RunCrawler(art::PtrVector<recob::Hit>& plnhits,
                 int plane, std::vector< ClusterStore >& tcl);
    
    int fNumPass;       ///< number of passes over the hit collection
    std::vector<int> fMaxHitsFit; ///< Max number of hits fitted
    std::vector<int> fMinHits;    ///< Min number of hits to make a cluster
    std::vector<int> fNHitsAve;   ///< number of hits used to compute fAveChg and fAveWid
    std::vector<float> fChiCut;   ///< stop adding hits to clusters if chisq too high
    std::vector<float> fKinkChiRat;  ///< Max consecutive chisq increase for the last 
                                     ///< 3 hits on the cluster
    std::vector<float> fKinkAngCut; ///< kink angle cut (dT/dW) made after fKinkChiRat
    std::vector<float> fWidCut;   ///< chisq cut for adding a hit to a cluster
    std::vector<float> fChgCut;   ///< charge difference cut for adding a hit to a cluster
    std::vector<int> fMaxWirSkip; ///< max number of wires that can be skipped while following
                                  ///< a cluster
    std::vector<int> fMinWirAfterSkip; ///< minimum number of hits on consecutive wires
                                     ///< after skipping
    std::vector<bool> fDoMerge;     ///< try to merge clusters?
    std::vector<float> fTimDelta;   ///< max time difference for matching

    // global cuts
    float fHitErrFac;   ///< hit time error = fHitErrFac * (EndTime - PeakTime)
    float fBEChgRat;    ///< Begin/End Charge ratio to determine cluster direction
                        ///< Set to 0 to turn off
    float fFudgeBigHits;   ///< increase the width of hits with large charge to
                           ///< prevent incorrect tracking stops. This is a workaround
                           ///< that wont be necessary when/if the hit reconstruction
                           ///< improves
    
    private:
    
    bool prt;

    std::map<int, int> FirstWirHit; ///< map of the first hit on each wire
                                    ///< returns 0 if no hits on the wire
                                    ///< returns -1 if the wire is dead
    std::map<int, float> hiterr2; ///< map of hit error of each hit
    std::map<int, float> hitwid; ///< map of hit width of each hit
    
    TVectorD clpar;     ///< cluster parameters for the current fit with
                        ///< origin at wire0
    TVectorD clparerr;  ///< cluster parameter errors
    float clChisq;     ///< chisq of the current fit
    float clBeginSlp;  ///< begin slope (= DS end = high wire number)
    int   clBeginWir;  ///< begin wire
    float clBeginTim;  ///< begin time
    float clBeginChg;  ///< begin average charge
    float clEndSlp;    ///< slope at the end   (= US end = low  wire number)
    int   clEndWir;    ///< begin wire
    float clEndTim;    ///< begin time
    float clEndChg;    ///< end average charge
    int clStopCode;     ///< code for the reason for stopping cluster tracking
                        ///< 0 = no signal on the next wire
                        ///< 1 = skipped too many occupied/dead wires
                        ///< 2 = failed the fMinWirAfterSkip cut
                        ///< 3 = ended on a kink. Fails fKinkChiRat
                        ///< 4 = failed the fChiCut cut
    
    // TLinearFitter arrays
    TLinearFitter *lf;
    Double_t *xwir;
    Double_t *ytim;
    Double_t *ytimerr;
    
    int fFirstWire;    ///< the first wire with a hit
    float fAveWid;  ///< average hit width at leading edge of cluster
    float fAveChg;     ///< average pulse height at leading edge of cluster
    int wire0;          ///< wire number origin of the fit => usually end wire
    
    int pass;
    
    std::vector<int> fcl2hits;  ///< vector of hits used in the cluster
    

    std::string fhitsModuleLabel;
    
    // Finds a hit on wire kwire, adds it to the cluster and re-fits it
    void cl2AddHit(art::PtrVector<recob::Hit>& plnhits, int kwire, bool phchk,
      bool& HitOK, bool& SigOK);
    // Fits the cluster hits in fcl2hits to a straight line
    void cl2Fit(art::PtrVector<recob::Hit>& plnhits);
    // Follows a trail of hits UpStream
    void cl2FollowUS(art::PtrVector<recob::Hit>& plnhits);
    // Stores cluster information in a temporary vector
    void cl2TmpStore(art::PtrVector<recob::Hit>& plnhits, 
      std::vector<ClusterStore>& tcl,int pass, int nClusters);
    // Compares two cluster combinations to see if they should be merged
    void cl2ChkMerge(art::PtrVector<recob::Hit>& plnhits,
      std::vector<ClusterStore>& tcl);
    // Checks merge for cluster cl2 within the bounds of cl1
    void cl2ChkMerge12(art::PtrVector<recob::Hit>& plnhits,
      std::vector<ClusterStore>& tcl, unsigned int it1, unsigned int it2, bool& didit);
    // Merges clusters cl1 and cl2
    void cl2DoMerge(art::PtrVector<recob::Hit>& plnhits, 
      std::vector<ClusterStore>& tcl, unsigned int it1, unsigned int it2);
    // Prints cluster information to the screen
    void cl2Print(art::PtrVector<recob::Hit>& plnhits, 
     std::vector<ClusterStore>& tcl);
    // Cleans the current cluster = add missing hits, remove bad ones
    void cl2Clean(art::PtrVector<recob::Hit>& plnhits);   
    // Sets the beginning and end direction of the cluster. This should
    // be called last before returning to the calling routine
    void cl2SetBeginEnd(art::PtrVector<recob::Hit>& plnhits,
      std::vector<ClusterStore>& tcl);
    

  }; // class ClusterCrawlerAlg


} // namespace cluster

#endif // ifndef CLUSTERCRAWLERALG_H
