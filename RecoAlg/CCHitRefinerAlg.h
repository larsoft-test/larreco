////////////////////////////////////////////////////////////////////////
// CCHitRefiner.h
//
// Cluster Crawler Hit Refiner class
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef CCHITREFINERALG_H
#define CCHITREFINERALG_H


#include "TMath.h"

#include <vector>

#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 

#include "Geometry/Geometry.h"
#include "RecoBase/Hit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "RecoAlg/ClusterCrawlerAlg.h"
#include "RecoAlg/CCHitFinderAlg.h"

namespace cluster {

  class CCHitRefinerAlg {
  
  public:

    CCHitRefinerAlg(fhicl::ParameterSet const& pset);
    virtual ~CCHitRefinerAlg();

    void reconfigure(fhicl::ParameterSet const& pset);

    void RunCCHitRefiner(std::vector<CCHitFinderAlg::CCHit>& allhits,
      CCHitFinderAlg::HitCuts& hitcuts,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl,
      std::vector<ClusterCrawlerAlg::VtxStore>& vtx, ClusterCrawlerAlg& fCCAlg);

    bool fRefineHits;   ///< refine hits?
    float fBEChgRat;    ///< Begin/End Charge ratio to determine cluster direction
                        ///< Set to 0 to turn off
    
  private:
    
    
    bool prt;
    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;
    
    unsigned short plane;
    std::vector< std::pair<short, short> > WireHitRange;
    unsigned short fFirstWire;
    unsigned short fLastWire;

    // fit n Gaussians possibly with bounds setting (parmin, parmax)
    void FitNG(unsigned short nGaus, unsigned short tstart, 
      unsigned short tend, std::vector<float>& signal, 
      CCHitFinderAlg::HitCuts& hitcuts);
    // parameters, errors, lower limit, upper limits for FitNG
    std::vector<double> par;
    std::vector<double> parerr;
    std::vector<double> parlo;
    std::vector<double> parhi;
    float chidof;
    float fitChiDOF;

    void ChkTopo1Vtx(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      CCHitFinderAlg::HitCuts& hitcuts,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl, 
      std::vector<ClusterCrawlerAlg::VtxStore>& vtx, unsigned int iv,
      std::vector<unsigned short>& clusters);

    // Gets the signal Region Above Threshold
    void GetRAT(std::vector<CCHitFinderAlg::CCHit>& allhits,
      unsigned short inwire, unsigned short intime, std::vector<float>& rat,
      unsigned short& firsttick);
    // Sets the beginning and end direction of the cluster
    void cl2SetBeginEnd(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl);

  }; // class CCHitRefinerAlg

} // cchit

#endif // ifndef CCHITREFINERALG_H
