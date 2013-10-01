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

//#include "RecoBase/Hit.h"

namespace recob { 
  class Hit;
  class Cluster; 
}

    struct houghCorner
    {
      double strength=0;
      double p0=0;
      double p1=0;
      houghCorner(double strengthTemp=0,
          double p0Temp=0,
          double p1Temp=0)
      {
        strength=strengthTemp;
        p0=p0Temp;
        p1=p1Temp;
      }

      bool operator < (const houghCorner& houghCornerComp) const
      {
        return (strength < houghCornerComp.strength);
      }
    };


    // This stores information about merged lines
    struct mergedLines
    {
      double totalQ=0;
      double pMin0=0;
      double pMin1=0;
      double pMax0=0;
      double pMax1=0;
      int clusterNumber=-999999;
      double showerLikeness=0;
      mergedLines (double totalQTemp=0,
          double pMin0Temp=0,
          double pMin1Temp=0,
          double pMax0Temp=0,
          double pMax1Temp=0,
          double clusterNumberTemp=-999999,
          double showerLikenessTemp=0)
      {
        totalQ=totalQTemp;
        pMin0=pMin0Temp;
        pMin1=pMin1Temp;
        pMax0=pMax0Temp;
        pMax1=pMax1Temp;
        clusterNumber=clusterNumberTemp;
        showerLikeness=showerLikenessTemp;
      }
    };


    // This stores information about merged lines
    struct showerLine
    {
      double slope=0;
      double intercept=0;
      int clusterNumber=-999999;
      int oldClusterNumber=-999999;
      int lineSize=0;
      int iMinWire=0;
      int iMaxWire=0;
      bool showerMerged=false;
      showerLine (double slopeTemp=0,
          double interceptTemp=0,
          int clusterNumberTemp=-999999,
          int oldClusterNumberTemp=-999999,
          int lineSizeTemp=0,
          int iMinWireTemp=0,
          int iMaxWireTemp=0)
      {
        slope=slopeTemp;
        intercept=interceptTemp;
        clusterNumber=clusterNumberTemp;
        oldClusterNumber=oldClusterNumberTemp;
        lineSize=lineSizeTemp;
        iMinWire=iMinWireTemp;
        iMaxWire=iMaxWireTemp;
      }
      bool operator < (const showerLine& showerLineComp) const
      {
        return (lineSize < showerLineComp.lineSize);
      }
    };



    struct lineSlope
    {
      int clusterNumber=999999;
      int oldClusterNumber=999999;
      double clusterSlope=999999;
      double clusterIntercept=999999;
      double totalQ=-999999;
      double pMin0=999999;
      double pMin1=999999;
      double pMax0=-999999;
      double pMax1=-999999;
      double iMinWire=999999;
      double iMaxWire=-999999;
      double minWire=999999;
      double maxWire=-999999;
      double isolation=-999999;
      double showerLikeness=-999999;
      bool merged=false;
      bool showerMerged=false;
      std::vector<art::Ptr<recob::Hit>> hits;
      lineSlope(unsigned int num=999999, 
          double slope=999999, 
          double intercept=999999,
          double totalQTemp=-999999,
          double Min0=999999, 
          double Min1=999999, 
          double Max0=-999999, 
          double Max1=-999999,
          int    iMinWireTemp=999999,
          int    iMaxWireTemp=-999999,
          int    minWireTemp=999999,
          int    maxWireTemp=-999999,
          std::vector<art::Ptr<recob::Hit>> hitsTemp=std::vector<art::Ptr<recob::Hit>>())
      {
        clusterNumber = num;
        oldClusterNumber = num;
        clusterSlope = slope;
        clusterIntercept = intercept;
        totalQ=totalQTemp;
        pMin0 = Min0;
        pMin1 = Min1;
        pMax0 = Max0;
        pMax1 = Max1;
        iMinWire = iMinWireTemp;
        iMaxWire = iMaxWireTemp;
        minWire = minWireTemp;
        maxWire = maxWireTemp;
        merged = false;
        showerMerged = false;
        showerLikeness = 0;
        hits.swap(hitsTemp);
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
    void GetEquation( float row, float col, float &rho, float &theta)
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
    int m_rowLength;
    int m_numAccumulated;
    int m_rhoResolutionFactor;
    int m_numAngleCells;
    std::vector<float> m_cosTable;
    std::vector<float> m_sinTable;
    bool DoAddPoint(int x, int y);
    int  DoAddPointReturnMax(int x, int y, int *yMax, int *xMax, int minHits);
    bool DoSubtractPoint(int x, int y);


  }; // class HoughTransform  





  class HoughBaseAlg {
    
  public:
    
    HoughBaseAlg(fhicl::ParameterSet const& pset); 
    virtual ~HoughBaseAlg();

    size_t FastTransform(std::vector<art::Ptr<recob::Cluster> >         & clusIn,
			 std::vector<recob::Cluster>                    & ccol,  
			 std::vector< art::PtrVector<recob::Hit> >      & clusHitsOut,
			 art::Event                                const& evt,
			 std::string                               const& label);

    size_t Transform(std::vector<art::Ptr<recob::Hit> > const& hits,
                     std::vector<unsigned int>     *fpointId_to_clusterId,
                     unsigned int clusterId, // The id of the cluster we are examining
                     int *nClusters,
                     std::vector<unsigned int> corners,
                     std::vector<lineSlope> *linesFound);
    
    
    // interface to look for lines only on a set of hits,without slope and totalQ arrays
    size_t FastTransform(std::vector < art::Ptr < recob::Hit > >                 & clusIn,
     	             std::vector< art::PtrVector<recob::Hit> >      & clusHitsOut
		     );
    
    // interface to look for lines only on a set of hits
    size_t FastTransform(std::vector < art::Ptr < recob::Hit > >                 & clusIn,
     	             std::vector< art::PtrVector<recob::Hit> >      & clusHitsOut,
		     std::vector<double> &slope, std::vector<double> &totalQ
		     );
    

    size_t Transform(std::vector<art::Ptr<recob::Hit> > const& hits);

    size_t Transform(std::vector< art::Ptr<recob::Hit> > const& hits,
		     double                                   & slope,
		     double                                   & intercept);

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
    double fMissedHitsDistance;            ///< Distance between hits in a hough line before a hit is considered missed
    double fMissedHitsToLineSize;          ///< Ratio of missed hits to line size for a line to be considered a fake




    


  protected:

    friend class HoughTransformClus;
  };
  
  
}// namespace

#endif // HOUGHBASEALG_H
