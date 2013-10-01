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
      std::vector<std::pair<double,double> > pHit;
      std::vector<std::pair<double,double> > pHitChargeSigma;
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
    // Note, m_accum is a vector of associative containers, 
    // the vector elements are called by rho, 
    // theta is the container key, 
    // the number of hits is the value corresponding to the key
    std::vector<std::map<int,int> > m_accum;  // column=rho, row=theta
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

    size_t FastTransform(std::vector<art::Ptr<recob::Cluster> >         & clusIn,
			 std::vector<recob::Cluster>                    & ccol,  
			 std::vector< art::PtrVector<recob::Hit> >      & clusHitsOut,
			 art::Event                                const& evt,
			 std::string                               const& label);

    // clusterId is the id of the cluster we are examining
    size_t Transform(std::vector<art::Ptr<recob::Hit> > const& hits,
                     std::vector<unsigned int>                *fpointId_to_clusterId,
                     unsigned int                              clusterId, 
                     int                                      *nClusters,
                     std::vector<unsigned int>                 corners);
    
    // interface to look for lines only on a set of hits,without slope and totalQ arrays
    size_t FastTransform(std::vector < art::Ptr < recob::Hit > >                 & clusIn,
     	             std::vector< art::PtrVector<recob::Hit> >                   & clusHitsOut);
    
    // interface to look for lines only on a set of hits
    size_t FastTransform(std::vector < art::Ptr < recob::Hit > >                 & clusIn,
			 std::vector< art::PtrVector<recob::Hit> >               & clusHitsOut,
			 std::vector<double>                                     & slope, 
			 std::vector<double>                                     & totalQ);
    

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
                                           ///< that fall into the "correct" bin will be 
                                           ///< small and consistent with noise.
    double fMaxDistance;                   ///< Max distance that a hit can be from a line to be 
                                           ///< considered part of that line
    double fMaxSlope;                      ///< Max slope a line can have
    int    fRhoZeroOutRange;               ///< Range in rho over which to zero out area around 
                                           ///< previously found lines in the accumulator
    int    fThetaZeroOutRange;             ///< Range in theta over which to zero out area around 
                                           ///< previously found lines in the accumulator
    int    fRhoResolutionFactor;           ///< Factor determining the resolution in rho
    int    fPerCluster;                    ///< Tells the original Hough algorithm to look at 
                                           ///< clusters individually, or all hits at once
    int    fMissedHits;                    ///< Number of wires that are allowed to be missed before 
                                           ///< a line is broken up into segments
    double fMissedHitsDistance;            ///< Distance between hits in a hough line before a 
                                           ///< hit is considered missed
    double fMissedHitsToLineSize;          ///< Ratio of missed hits to line size for a line to 
                                           ///< be considered a fake
    int    fDoFuzzyRemnantMerge;           ///< Tell the algorithm to merge fuzzy cluster remnants
                                           ///< into showers or tracks (0-off, 1-on)
    int    fDoHoughLineMerge;              ///< Turns on Hough line merging (0-off, 1-on)
    double fHoughLineMergeAngle;           ///< Max angle between slopes before two lines are merged
                                           ///<  (muon tracks), only for fuzzy clustering
    double fShowerHoughLineMergeAngle;     ///< Max angle between slopes before two lines are merged,
                                           ///<  for lines in shower line regions
    int    fDoShowerHoughLineMerge;        ///< Turns on shower Hough line merging (0-off, 1-on)
                                           ///< for (electron showers), only for fuzzy clustering
    int    fDoChargeAsymAngleMerge;        ///< Turn on cut on product of charge asymmetry and sin of
                                           ///<  angle between slopes of lines
    double fChargeAsymAngleCut;            ///< Cut on product of charge asymmetry and sin of angle
                                           ///<  between slopes of lines
    double fSigmaChargeAsymAngleCut;       ///< Cut on product of charge asymmetry and sin of angle
                                           ///<  between slopes of lines
    double fChargeAsymAngleCutoff;         ///< Distance between lines before cut on product of 
                                           ///< charge asymmetry and sin of angle between slopes of lines
                                           ///< is applied
    double fShowerWidthAngle;              ///< Half of the angle defining how wide a shower is 
    double fHoughLineMergeCutoff;          ///< Max distance between Hough lines before two lines 
                                           ///< are merged (muon tracks), only for fuzzy clustering
    double fShowerHoughLineMergeCutoff;    ///< Max distance between Hough lines before two lines 
                                           ///< are merged (electron showers),
    int    fDoShowerHoughLineInterceptMerge;///< Turns on Hough line merging for shower like lines, 
                                           ///< merging if they intercept (0-off, 1-on) 
                                           ///< line closest to each other
    double fShowerLikenessCut;             ///< Cut on shower likeness (larger the more shower like, 
                                           ///< smaller the less shower like)


    void mergeHoughLinesBySegment(unsigned int k,
				  std::vector<lineSlope> *linesFound, 
				  double xyScale,
				  int mergeStyle,
				  double wire_dist,
				  double tickToDist);
    
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
    
    double DistanceBetweenHits(art::Ptr<recob::Hit> hit0,
			       art::Ptr<recob::Hit> hit1,
			       double wire_dist,
			       double tickToDist);

  protected:

    friend class HoughTransformClus;
  };
  
  
}// namespace

#endif // HOUGHBASEALG_H
