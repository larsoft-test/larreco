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
      int ID;                  // Cluster ID = 1000 * pass + number
      double slpstart;  // slope at the start (= DS end = high wire number)
      double slpend;    // slope at the end   (= US end = low  wire number)
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
    
    int fNumPass;     ///< number of passes over the hit collection
    int fMaxHitsFit[5];  ///< Max number of hits fitted
    int fMinHits[5];  ///< Min number of hits to make a cluster
    double fHitErrFac;///< hit time error = fHitErrFac * (EndTime - PeakTime)
    double fChiCut;   ///< stop adding hits to clusters if chisq too high
    double fSigCut;   ///< chisq cut for adding a hit to a cluster
    double fChgCut;   ///< charge difference cut for adding a hit to a cluster
    int fMaxWirSkip;  ///< max number of wires that can be skipped while following
                      ///< a cluster
    bool fDoMerge;    ///< try to merge clusters?
    int fWirDelta;    ///< max wire difference between start/end for merging
    int fSlpDelta;    ///< max slope difference for merging
    double fChgDelta; ///< max charge ratio difference for merging
    int fTimDelta;    ///< max time difference for matching
    
    private:
    
    bool prt;

    std::map<int, int> FirstWirHit; ///< map of the first hit on each wire
                                    ///< returns 0 if no hits on the wire
                                    ///< returns -1 if the wire is dead
    std::map<int, double> hiterr2; ///< map of hit error for each hit
    
    TVectorD clpar;     ///< cluster parameters for the current fit with
                        ///< origin at wire0
    TVectorD clparerr;  ///< cluster parameter errors
    int wire0;          ///< wire number origin of the fit
    double clchisq;     ///< chisq of the current fit
    double clslpstart;  ///< slope at the start (= DS end = high wire number)
    double clslpend;    ///< slope at the end   (= US end = low  wire number)
    
    // TLinearFitter arrays
    TLinearFitter *lf;
    Double_t *xwir;
    Double_t *ytim;
    Double_t *ytimerr;
    
    int fFirstWire;    ///< the first wire with a hit
    double fAveWid;  ///< average hit width at leading edge of cluster
    double fAveChg;     ///< average pulse height at leading edge of cluster
    
    int pass;
    
    std::vector<int> fcl2hits;  ///< vector of hits used in the cluster
    

    std::string fhitsModuleLabel;
    
    // Finds a hit on wire kwire, adds it to the cluster and re-fits it
    int cl2AddHit(art::PtrVector<recob::Hit>& plnhits, int kwire, bool phchk,
      bool SigOK);
    // Fits the cluster hits in fcl2hits to a straight line
    void cl2Fit(art::PtrVector<recob::Hit>& plnhits);
    // Fits hits on wires lowire to hiwire to a straight line w origin at wire0
//    void cl2Fit(art::PtrVector<recob::Hit>& plnhits, int lowire, int hiwire);
    // Kills a cluster
    void cl2Kill(art::PtrVector<recob::Hit>& plnhits);
    // Follows a trail of hits UpStream
    void cl2FollowUS(art::PtrVector<recob::Hit>& plnhits);
    // Re-follows the current cluster to pick up missing hits
    void cl2ReFollow(art::PtrVector<recob::Hit>& plnhits);
    // Stores cluster information in a temporary vector
    void cl2TmpStore(art::PtrVector<recob::Hit>& plnhits, 
      std::vector<ClusterStore>& tcl,int pass, int nClusters);
    // Compares two cluster combinations to see if they should be merged
    void cl2ChkMerge(art::PtrVector<recob::Hit>& plnhits,
      std::vector<ClusterStore>& tcl);
    // Merges clusters cl1 and cl2
    void cl2DoMerge(art::PtrVector<recob::Hit>& plnhits, 
      std::vector<ClusterStore>& tcl, unsigned int it1, unsigned int it2);
    // Prints cluster information to the screen
    void cl2Print(art::PtrVector<recob::Hit>& plnhits, 
     std::vector<ClusterStore>& tcl);
    // Cleans the current cluster = add missing hits, remove bad ones
    void cl2Clean(art::PtrVector<recob::Hit>& plnhits);    

  }; // class ClusterCrawlerAlg


} // namespace cluster

#endif // ifndef CLUSTERCRAWLERALG_H
