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

    // struct of temporary clusters

    struct ClusterStore {
      short ID;         // Cluster ID. ID < 0 = abandoned cluster
      short ProcCode;   // Processor code for debugging
      short Assn;       // coded pointer to associated clusters
      short StopCode;   // code for the reason for stopping cluster tracking
      float BeginSlp;   // beginning slope (= DS end = high wire number)
      float BeginSlpErr; // error
      short BeginWir;   // begin wire
      float BeginTim;   // begin time
      float BeginChg;   // beginning average charge
      float EndSlp;     // end slope (= US end = low  wire number)
      float EndSlpErr;  // error
      short EndWir;     // end wire
      float EndTim;     // ending time
      float EndChg;     // ending average charge
      short BeginVtx;   // ID of the begin vertex
      short EndVtx;     // ID of the end vertex
      std::vector<short> tclhits; // hits on the cluster
    };
    std::vector< ClusterStore > tcl;

    // struct of temporary vertices
    struct VtxStore {
      short Wire;
      float Time;
      float Wght;
      short Topo;
    };
    std::vector< VtxStore > vtx;

  class ClusterCrawlerAlg {
    public:

    ClusterCrawlerAlg(fhicl::ParameterSet const& pset);
    virtual ~ClusterCrawlerAlg();

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void RunCrawler(art::PtrVector<recob::Hit>& plnhits,
       int plane, std::vector< ClusterStore >& tcl, std::vector< VtxStore >& vtx);
    
    short fNumPass;                 ///< number of passes over the hit collection
    std::vector<short> fMaxHitsFit; ///< Max number of hits fitted
    std::vector<short> fMinHits;    ///< Min number of hits to make a cluster
    std::vector<short> fNHitsAve;   ///< number of US hits used to compute fAveChg and fAveWid
                                    ///< set to > 2 to do a charge fit using fNHitsAve hits
    std::vector<float> fChiCut;     ///< stop adding hits to clusters if chisq too high
    std::vector<float> fKinkChiRat; ///< Max consecutive chisq increase for the last 
                                    ///< 3 hits on the cluster
    std::vector<float> fKinkAngCut; ///< kink angle cut made after fKinkChiRat
    std::vector<float> fWidCut;     ///< chisq cut for adding a hit to a cluster
    std::vector<float> fChgCut;     ///< charge difference cut for adding a hit to a cluster
    std::vector<short> fMaxWirSkip; ///< max number of wires that can be skipped while following
                                    ///< a cluster
    std::vector<short> fMinWirAfterSkip; ///< minimum number of hits on consecutive wires
                                    ///< after skipping
    std::vector<bool> fDoMerge;     ///< try to merge clusters?
    std::vector<float> fTimeDelta;  ///< max time difference for matching
    std::vector<float> fTimeDeltaLA; ///< max time difference for matching
                                     ///< large angle clusters, abs(slope) > 20

    // global cuts and parameters
    float fHitErrFac;   ///< hit time error = fHitErrFac * (EndTime - PeakTime)
                        ///< used for cluster fit
    float fHitWidFac;   ///< hit width = fHitWidFac * (EndTime - PeakTime)
                        ///< used to decide if there is a signal near the projected
                        ///< cluster position
    float fBEChgRat;    ///< Begin/End Charge ratio to determine cluster direction
                        ///< Set to 0 to turn off
    float fFudgeBigHits;   ///< increase the width of hits with large charge to
                           ///< prevent incorrect tracking stops. This is a workaround
                           ///< that wont be necessary when/if the hit reconstruction
                           ///< improves
    float fPairAngCut;     ///< Close pair angle cut. Set <= 0 to turn off
                           ///< Applied to clusters with >= 10 hits
    float fCurlyMergeAngCut;  ///< Run final merge on short curly clusters?
    bool fDoVertex;        ///< run vertexing code
    short fDebugWire;  ///< set to the Begin Wire and Hit of a cluster to print
    short fDebugHit;   ///< out detailed information while crawling
    
    private:
    
    bool prt;
    short NClusters;
    bool ArgoFix;

    
    TVectorD clpar;     ///< cluster parameters for the current fit with
                        ///< origin at wire0
    TVectorD clparerr;  ///< cluster parameter errors
    TVectorD clchgpar;  ///< cluster charge fit parameters
    float clChisq;     ///< chisq of the current fit
    float clBeginSlp;  ///< begin slope (= DS end = high wire number)
    float clBeginSlpErr; 
    short clBeginWir;  ///< begin wire
    float clBeginTim;  ///< begin time
    float clBeginChg;  ///< begin average charge
    float clEndSlp;    ///< slope at the end   (= US end = low  wire number)
    float clEndSlpErr; 
    short clEndWir;    ///< begin wire
    float clEndTim;    ///< begin time
    float clEndChg;    ///< end average charge
    short clStopCode;     ///< code for the reason for stopping cluster tracking
                        ///< 0 = no signal on the next wire
                        ///< 1 = skipped too many occupied/dead wires
                        ///< 2 = failed the fMinWirAfterSkip cut
                        ///< 3 = ended on a kink. Fails fKinkChiRat
                        ///< 4 = failed the fChiCut cut
                        ///< 5 = cluster split by cl2ChkPair
    short clProcCode;     ///< Processor code = pass number
                        ///< +   10 cl2ChkMerge
                        ///< +  100 cl2ChkMerge12
                        ///< + 1000 cl2ChkPair
                        ///< + 2000 failed pass N cuts but passes pass N=1 cuts
                        ///< +10000 cl2CurlyMerge
    short clAssn;         ///< index of a parent cluster. -1 if no parent.
                        ///< Parent clusters are not associated with daughters
    
    // TLinearFitter arrays
    TLinearFitter *lf;
    Double_t *xwir;
    Double_t *ytim;
    Double_t *ytimerr;
    Double_t *ychg;
    Double_t *ychgerr;
    
    short fFirstWire;    ///< the first wire with a hit
    short fLastWire;      ///< the last wire with a hit
    float fAveWid;  ///< average hit width at leading edge of cluster
    float fAveChg;     ///< average pulse height at leading edge of cluster
    short wire0;          ///< wire number origin of the fit => usually end wire
    
    short pass;
    float ScaleF;     ///< scale factor from Tick/Wire to dx/du

    // vector of pairs of first (.first) and last+1 (.second) hit on each wire
    // in the range fFirstWire to fLastWire. A value of -2 indicates that there
    // are no hits on the wire. A value of -1 indicates that the wire is dead
    std::vector< std::pair<short, short> > WireHitRange;
    std::map<short, float> hiterr2; ///< map of hit error of each hit
    std::map<short, float> hitwid; ///< map of hit width of each hit
    
    std::vector<short> fcl2hits;  ///< vector of hits used in the cluster

    std::string fhitsModuleLabel;
    
/*
    // Sorts clusters in tcl by decreasing number of hits, while ignoring
    // abandoned clusters
    void cl2SortByLength(std::vector<ClusterStore>& tcl,
        std::vector<unsigned int>& sortindex);
*/
    // Finds a hit on wire kwire, adds it to the cluster and re-fits it
    void cl2AddHit(art::PtrVector<recob::Hit>& plnhits, short kwire, bool phchk,
      bool& HitOK, bool& SigOK);
    // Fits the cluster hits in fcl2hits to a straight line
    void cl2Fit(art::PtrVector<recob::Hit>& plnhits);
    // Fits the middle of a temporary cluster it1 using hits iht to iht + nhit
    void cl2FitMid(art::PtrVector<recob::Hit>& plnhits, 
      std::vector<ClusterStore>& tcl, unsigned int it1, short iht, short nhit);
    // Follows a trail of hits UpStream
    void cl2FollowUS(art::PtrVector<recob::Hit>& plnhits);
    // Stores cluster information in a temporary vector
    void cl2TmpStore(art::PtrVector<recob::Hit>& plnhits, 
      std::vector<ClusterStore>& tcl);
    // attaches hits to the ends of tcl clusters
    void cl2Clean(art::PtrVector<recob::Hit>& plnhits, 
      std::vector<ClusterStore>& tcl);
    // Prepares close pair clusters
    void cl2ChkPair(art::PtrVector<recob::Hit>& plnhits,
      std::vector<ClusterStore>& tcl);
    // Splits close pair clusters
    void cl2DoSplit(art::PtrVector<recob::Hit>& plnhits,
      std::vector<ClusterStore>& tcl, unsigned int it1, unsigned int it2);
    // Compares two cluster combinations to see if they should be merged
    void cl2ChkMerge(art::PtrVector<recob::Hit>& plnhits,
      std::vector<ClusterStore>& tcl);
    // Checks merge for cluster cl2 within the bounds of cl1
    void cl2ChkMerge12(art::PtrVector<recob::Hit>& plnhits,
      std::vector<ClusterStore>& tcl, unsigned int it1, unsigned int it2, bool& didit);
    // Merges clusters cl1 and cl2
    void cl2DoMerge(art::PtrVector<recob::Hit>& plnhits, 
      std::vector<ClusterStore>& tcl, unsigned int it1, unsigned int it2, short ProcCode);
    void cl2CurlyMerge(art::PtrVector<recob::Hit>& plnhits,
      std::vector<ClusterStore>& tcl);
    // Prints cluster information to the screen
    void cl2Print(art::PtrVector<recob::Hit>& plnhits, 
     std::vector<ClusterStore>& tcl);
    // Sets the beginning and end direction of the cluster. This should
    // be called last before returning to the calling routine
    void cl2SetBeginEnd(art::PtrVector<recob::Hit>& plnhits,
      std::vector<ClusterStore>& tcl);
    // make 2D vertices
    void cl2DoVertex(art::PtrVector<recob::Hit>& plnhits,
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx);
    // check for a signal on all wires between two points
    void cl2ChkSignal(art::PtrVector<recob::Hit>& plnhits,
      short wire1, float time1, short wire2, float time2, bool& SigOK);
    // check a vertex (vw, fvt) made with clusters it1, and it2 against the
    // vector of existing clusters
    void cl2ChkVertex(art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
        short vw, float fvt, unsigned int it1, unsigned int it2, short topo);
    // try to attach a cluster to an existing vertex
    void cl2ClsVertex(art::PtrVector<recob::Hit>& plnhits, 
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
        unsigned int it2);

  }; // class ClusterCrawlerAlg


} // namespace cluster

#endif // ifndef CLUSTERCRAWLERALG_H
