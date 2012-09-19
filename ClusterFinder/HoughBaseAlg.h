////////////////////////////////////////////////////////////////////////
// HoughBaseAlg.h
//
// HoughBaseAlg class
//
// Ben Carls (bcarls@fnal.gov)
//
////////////////////////////////////////////////////////////////////////
#ifndef HOUGHBASEALG_H
#define HOUGHBASEALG_H

#include "TMath.h"
#include <vector>

#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 

namespace recob { 
  class Hit;
  class Cluster; 
}

    struct lineSlope
    {
      unsigned int clusterNumber;
      unsigned int oldClusterNumber;
      double clusterSlope;
      double clusterIntercept;
      double pMin0;
      double pMin1;
      double pMax0;
      double pMax1;
      double isolation;
      bool merged;
      lineSlope(unsigned int num, 
          double slope, 
          double intercept, 
          double Min0, 
          double Min1, 
          double Max0, 
          double Max1,
          double signalToBkg)
      {
        clusterNumber = num;
        oldClusterNumber = num;
        clusterSlope = slope;
        clusterIntercept = intercept;
        pMin0 = Min0;
        pMin1 = Min1;
        pMax0 = Max0;
        pMax1 = Max1;
        merged = false;
        isolation = signalToBkg;
      }
    };



namespace cluster {
   
  class HoughTransform {
  public:
    
    HoughTransform();
    ~HoughTransform();
     
    void Init(int dx, int dy, int rhoresfact, int numACells);
    bool AddPoint(int x, int y);
    bool AddPointReturnMax(int x, int y);
    int  AddPointReturnMax(int x, int y, int *yMax, int *xMax, int minHits);
    bool SubtractPoint(int x, int y);
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

    void reconfigure(fhicl::ParameterSet const& pset);

    private:
         
    int m_dx;
    int m_dy;
    // Note, m_accum is a vector of associative containers, the vector elements are called by rho, theta is the container key, the number of hits is the value corresponding to the key
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
    int  DoAddPointReturnMax(int x, int y, int *yMax, int *xMax, int minHits);
    bool DoSubtractPoint(int x, int y);


  }; // class HoughTransform  





  class HoughBaseAlg {
    
  public:
    
    HoughBaseAlg(fhicl::ParameterSet const& pset); 
    virtual ~HoughBaseAlg();
        

    // Copy the functions over entirely to the base class
    size_t Transform(std::vector<art::Ptr<recob::Cluster> >           & clusIn,
     	             std::vector<recob::Cluster>                      & ccol,  
		     std::vector< art::PtrVector<recob::Hit> >      & clusHitsOut,
		     art::Event                                const& evt,
		     std::string                               const& label);

    size_t Transform(std::vector<art::Ptr<recob::Hit> >& hits,
                     std::vector<unsigned int>     *fpointId_to_clusterId,
                     unsigned int clusterId, // The id of the cluster we are examining
                     int *nClusters,
                     std::vector<unsigned int> corners);

    size_t Transform(std::vector<art::Ptr<recob::Hit> >& hits);

    size_t Transform(std::vector< art::Ptr<recob::Hit> >& hits,
		     double                             & slope,
		     double                             & intercept);

    virtual void reconfigure(fhicl::ParameterSet const& pset);
         
  protected:

    void HLSSaveBMPFile(char const*, unsigned char*, int, int);

  private:

    int    fMaxLines;                      ///< Max number of lines that can be found 
    int    fMinHits;                       ///< Min number of hits in the accumulator to consider 
                                           ///< (number of hits required to be considered a line).
    int    fSaveAccumulator;               ///< Save bitmap image of accumulator for debugging?
    int    fNumAngleCells;                 ///< Number of angle cells in the accumulator 
                                           ///< (a measure of the angular resolution of the line finder). 
                                           ///< If this number is too large than the number of votes 
                                           ///< that fall into the "correct" bin will be small and consistent with noise.
    double fMaxDistance;                   ///< Max distance that a hit can be from a line to be considered part of that line
    double fMaxSlope;                      ///< Max slope a line can have
    int    fRhoZeroOutRange;               ///< Range in rho over which to zero out area around previously found lines in the accumulator
    int    fThetaZeroOutRange;             ///< Range in theta over which to zero out area around previously found lines in the accumulator
    int    fRhoResolutionFactor;           ///< Factor determining the resolution in rho
    int    fPerCluster;                    ///< Tells the original Hough algorithm to look at clusters individually, or all hits
                                           ///< at once
    int    fMissedHits;                    ///< Number of wires that are allowed to be missed before a line is broken up into
                                           ///< segments
    double fEndPointCutoff;                ///< Max distance from a Hough line's end point and an end point from the end point
                                           ///< finder, this is obsolete
    double fHoughLineMergeAngle;           ///< Max angle between slopes before two lines are merged (muon tracks), only for fuzzy clustering
    double fParaHoughLineMergeAngle;       ///< Max angle between slopes before two lines are merged, they should 
                                           ///< be close to parallel (electron showers), only for fuzzy clustering
    double fLineIsolationCut;              ///< Cut on the Hough line isolation, only for fuzzy clustering
    double fHoughLineMergeCutoff;          ///< Max distance between Hough lines before two lines are merged (muon tracks), only for fuzzy clustering
    double fParaHoughLineMergeCutoff;      ///< Max distance between Hough lines before two lines are merged (electron showers),
                                           ///< they are generally farther apart from each other, only for fuzzy clustering
     
    void mergeHoughLines(unsigned int k,
        std::vector<lineSlope> *linesFound, 
        std::vector<int> *newClusNum, 
        std::vector<double> *newClusDist, 
        double xyScale);
    void mergeHoughLinesBySegment(unsigned int k,
        std::vector<lineSlope> *linesFound, 
        std::vector<int> *newClusNum, 
        std::vector<double> *newClusDist, 
        double xyScale);
    int mergeParaHoughLines(unsigned int k,
        std::vector<lineSlope> *linesFound, 
        std::vector<int> *newClusNum, 
        std::vector<double> *newClusDist, 
        double xyScale);
    void mergeParaHoughLinesBySegment(unsigned int k,
        std::vector<lineSlope> *linesFound, 
        std::vector<int> *newClusNum, 
        std::vector<double> *newClusDist, 
        double xyScale);
     

    //std::vector<lineSlope> linesFound;
    double HoughLineDistance(double p0MinLine1, 
        double p1MinLine1, 
        double p0MaxLine1, 
        double p1MaxLine1, 
        double p0MinLine2, 
        double p1MinLine2, 
        double p0MaxLine2, 
        double p1MaxLine2);
    bool   HoughLineIntersect(double x11,
        double  y11,
        double  x12,
        double  y12,
        double  x21,
        double  y21,
        double  x22,
        double  y22);
    double PointSegmentDistance(double px,
        double  py,
        double  x1,
        double  y1,
        double  x2,
        double  y2);




    


  protected:

    friend class HoughTransformClus;
  };
  
  
}// namespace

#endif // HOUGHBASEALG_H
