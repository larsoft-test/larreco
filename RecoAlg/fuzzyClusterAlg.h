/////////////////////////////////////////////////////////////////
//  \filefuzzyClusterAlg.h
//  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////
#ifndef fuzzyClusterALG_H
#define fuzzyClusterALG_H
#include <vector>
#include <cmath>
#include <iostream>
#include <stdint.h>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

#include "RecoAlg/HoughBaseAlg.h"
#include "RecoAlg/DBScanAlg.h"
#include "Geometry/Geometry.h"

#include "TMatrixF.h"
#include "TMatrixDEigen.h"
#include "TVectorF.h"
#include "TVector.h"


class TH1F;

namespace recob { class Hit; }

namespace cluster{

  //--------------------------------------------------------------- 
  class fuzzyClusterAlg {
  public:
    
    
    fuzzyClusterAlg(fhicl::ParameterSet const& pset);
    virtual ~fuzzyClusterAlg();
    
    void reconfigure(fhicl::ParameterSet const& p);
    void InitFuzzy(std::vector<art::Ptr<recob::Hit> >& allhits, std::set<uint32_t> badChannels);
    // Three differnt version of the clustering code
    void run_fuzzy_cluster(std::vector<art::Ptr<recob::Hit> >& allhits);     
   
    //Functions for fuzzy clustering
    void computeCentroids(int k);
    void computeCentroids2(int k);
    bool mergeClusters();
    bool updateMembership(int k);
    inline bool canStop(){
      double epsilon = 0.01;
      TMatrixF diffMatrix = fpsMembership - fpsNewMembership;
      double difference = diffMatrix.Norm1();
      return difference < epsilon;
    }


    std::vector<std::vector<unsigned int> > fclusters;               ///< collection of something
    std::vector<std::vector<float> >       fps;                     ///< the collection of points we are working on     
    std::vector<unsigned int>               fpointId_to_clusterId;   ///< mapping point_id -> clusterId     
    std::vector<std::vector<float> >       fsim;                    ///<
    std::vector<std::vector<float> >       fsim2;            	     ///<
    std::vector<std::vector<float> >       fsim3;            	     ///<
    float fMaxWidth;

    //Matrices for Ben's fuzzy cluster
    TMatrixF                         fpsMat;

    // Get functions and structures from HoughBaseAlg
    //friend class HoughBaseAlg;

  private:
    
    // The fuzzyness factor needed for fuzzy clustering, commonly known as just m
    float fFuzzifier;
    // The maximum number of clusters to try, needed for fuzzy clustering
    int fMaxNumClusters;
    // The maximum number of iterations to try, needed for fuzzy clustering
    int nIterations;
    // The parameter beta used in determining the radius of a cluster
    float fBeta;
    int    fDoFuzzyRemnantMerge;           ///< Tell the algorithm to merge fuzzy cluster remnants into showers or tracks (0-off, 1-on)
    int    fDoHoughLineMerge;              ///< Turns on Hough line merging (0-off, 1-on)
    double fHoughLineMergeAngle;           ///< Max angle between slopes before two lines are merged (muon tracks), only for fuzzy clustering
    double fShowerHoughLineMergeAngle;     ///< Max angle between slopes before two lines are merged, for lines in shower line regions
    int    fDoShowerHoughLineMerge;        ///< Turns on shower Hough line merging (0-off, 1-on)
                                           ///< for (electron showers), only for fuzzy clustering
    int    fDoChargeAsymAngleMerge;        ///< Turn on cut on product of charge asymmetry and sin of angle between slopes of lines
    double fChargeAsymAngleCut;            ///< Cut on product of charge asymmetry and sin of angle between slopes of lines
    double fSigmaChargeAsymAngleCut;       ///< Cut on product of charge asymmetry and sin of angle between slopes of lines
    double fChargeAsymAngleCutoff;         ///< Distance between lines before cut on product of charge asymmetry and sin of angle between slopes of lines
                                           ///< is applied
    double fFuzzyRemnantMergeCutoff;       ///< cut off on merging the fuzzy cluster remnants into the nearest shower or track 
    double fShowerWidthAngle;              ///< Half of the angle defining how wide a shower is 
    double fHoughLineMergeCutoff;          ///< Max distance between Hough lines before two lines are merged (muon tracks), 
                                           ///< only for fuzzy clustering
    double fShowerHoughLineMergeCutoff;    ///< Max distance between Hough lines before two lines are merged (electron showers),
    int    fDoShowerHoughLineInterceptMerge;///< Turns on Hough line merging for shower like lines, merging if they intercept (0-off, 1-on)
                                           ///< line closest to each other
    double fShowerLikenessCut;             ///< Cut on shower likeness (larger the more shower like, smaller the less shower like)

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

    double DistanceBetweenHits(
        art::Ptr<recob::Hit> hit0,
        art::Ptr<recob::Hit> hit1,
        double wire_dist,
        double tickToDist);



    // noise vector
    std::vector<bool>      fnoise;	
    std::vector<bool>      fvisited;					     
    std::vector<float>    fWirePitch;     ///< the pitch of the wires in each plane
    std::set<uint32_t>     fBadChannels;   ///< set of bad channels in this detector
    std::vector<uint32_t>  fBadWireSum;    ///< running total of bad channels. Used for fast intervening 
                                           ///< dead wire counting ala fBadChannelSum[m]-fBadChannelSum[n]. 

    //Matrices for Ben's fuzzy cluster
    TMatrixF                         fpsMembership;
    TMatrixF                         fpsNewMembership;
    TMatrixF                         fpsMembershipFinal;
    TMatrixF                         fpsCentroids;

    // Object used for Hough transforms
    HoughBaseAlg fHBAlg;        

    // Object used for DBScan
    DBScanAlg fDBScan;

    art::ServiceHandle<geo::Geometry> fGeom; ///< handle to geometry service

  }; // class fuzzyClusterAlg

} // namespace

#endif // ifndef fuzzyClusterALG_H
