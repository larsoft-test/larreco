////////////////////////////////////////////////////////////////////////
// $Id: HoughLineAlg.h,v 1.36 2010/09/15  bpage Exp $
//
// HoughLineAlg class
//
// josh
//
////////////////////////////////////////////////////////////////////////
#ifndef HOUGHLINEALG_H
#define HOUGHLINEALG_H

#include "TMath.h"
#include <vector>

#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 

#include "ClusterFinder/HoughBaseAlg.h"

namespace recob { 
  class Hit;
  class Cluster; 
}

namespace cluster {
   
  class HoughLineAlg:public HoughBaseAlg {
    
  public:
    
    HoughLineAlg(fhicl::ParameterSet const& pset); 
    virtual ~HoughLineAlg();
         
    size_t Transform(art::PtrVector<recob::Cluster>                 & clusIn,
     	             std::vector<recob::Cluster>               	    & ccol,  
		     std::vector< art::PtrVector<recob::Hit> >      & clusHitsOut,
		     art::Event                                const& evt,
		     std::string                               const& label);

    size_t Transform(std::vector< art::Ptr<recob::Hit> >& hits,
		     double                             & slope,
		     double                             & intercept);

    void reconfigure(fhicl::ParameterSet const& pset);
          
  private:
  
    int    fMaxLines;         ///< Max number of lines that can be found 
    int    fMinHits;          ///< Min number of hits in the accumulator to consider 
                              ///< (number of hits required to be considered a line).
    int    fSaveAccumulator;  ///< Save bitmap image of accumulator for debugging?
    int    fNumAngleCells;    ///< Number of angle cells in the accumulator 
                              ///< (a measure of the angular resolution of the line finder). 
                              ///< If this number is too large than the number of votes 
                              ///< that fall into the "correct" bin will be small and consistent with noise.
    double fMaxDistance;
    double fMaxSlope;
    int    fRhoZeroOutRange;
    int    fThetaZeroOutRange;
    int    fRhoResolutionFactor;
    int    fPerCluster;
    int    fMissedHits;
      
  protected:

    friend class HoughTransform;
  };
  
  
}// namespace

#endif // HOUGHLINEALG_H
