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

namespace recob { 
  class Hit;
  class Cluster; 
}

namespace cluster {
   
  class HoughTransform {
  public:
    
    HoughTransform();
    ~HoughTransform();
     
    void Init(int dx, int dy, int rhoresfact, int numACells);
    bool AddPoint(int x, int y);
    int  GetCell(int row, int col)            { return m_accum[row][col]; }
    void SetCell(int row, int col, int value) { m_accum[row][col] = value; }
    void IncrementCell(int row, int col)      { m_accum[row][col]++;}
    void GetAccumSize(int &numRows, int &numCols) 
    { 
      numRows = m_accum.size();
      numCols  = m_rowLength;
    }
    int NumAccumulated()                      { return m_numAccumulated; }
    void GetEquation(double row, double col, double &rho, double &theta)
    {
      theta = (TMath::Pi()*row)/m_numAngleCells;
      rho   = (col - (m_rowLength/2.))/m_rhoResolutionFactor;
    }
    int GetMax(int & xmax, int & ymax);
      
    private:
         
    int m_dx;
    int m_dy;
    std::vector<std::map<int,int> > m_accum;  // column=rho, row=theta
    //int distCenter;// \todo Why is this here? Only used locally by DoAddPoint
    //int lastDist;
    //int dist;
    //int stepDir; 
    //int cell;
    int m_rowLength;
    int m_numAccumulated;
    int m_rhoResolutionFactor;
    int m_numAngleCells;
    std::vector<double> m_cosTable;
    std::vector<double> m_sinTable;
    bool DoAddPoint(int x, int y);
  }; // class HoughTransform  

  class HoughLineAlg {
    
  public:
    
    HoughLineAlg(fhicl::ParameterSet const& pset); 
    virtual ~HoughLineAlg();
         
    size_t Transform(art::PtrVector<recob::Cluster>& clusIn,
     	             std::vector<recob::Cluster>& ccol);

    size_t Transform(std::vector< art::Ptr<recob::Hit> >& hits,
		     double &slope,
		     double &intercept);

    void reconfigure(fhicl::ParameterSet const& pset);
          
  private:
  
    void HLSSaveBMPFile(char const*, unsigned char*, int, int);

  
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
