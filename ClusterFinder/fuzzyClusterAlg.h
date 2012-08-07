/////////////////////////////////////////////////////////////////
//  \filefuzzyClusterAlg.h
//  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////
#ifndef fuzzyClusterALG_H
#define fuzzyClusterALG_H
#include <vector>
#include <cmath>
#include <iostream>

#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "ClusterFinder/RStarTree/RStarTree.h"
#include "ClusterFinder/HoughClusAlg.h"
#include "ClusterFinder/EndPointClusAlg.h"

#include "Geometry/geo.h"
#include "TMatrixTBase.h"


class TH1F;

namespace recob { class Hit; }

// RStarTree related infrastructure
//
// Our core objects have a physical extent (I.e. there are not
// points), but a R*-tree should be able to deal with that.
//class dbsPoint; // forward declaration
typedef RStarTree< unsigned int, 2, 32, 64 > RTree; // payload is just an index
typedef RTree::BoundingBox BoundingBox;

namespace cluster{

  //--------------------------------------------------------------- 
  class fuzzyClusterAlg {
  public:
    
    
    fuzzyClusterAlg(fhicl::ParameterSet const& pset);
    virtual ~fuzzyClusterAlg();
    
    void reconfigure(fhicl::ParameterSet const& p);
    void InitFuzzy(art::PtrVector<recob::Hit>& allhits, std::set<unsigned int> badChannels);
    // Three differnt version of the clustering code
    void run_fuzzy_cluster(art::PtrVector<recob::Hit>& allhits);     
    double getSimilarity(const std::vector<double> v1, const std::vector<double> v2); 
    std::vector<unsigned int> findNeighbors( unsigned int pid, double threshold, double threshold2);
    void computeSimilarity();
    void run_cluster();     
    double getSimilarity2(const std::vector<double> v1, const std::vector<double> v2); 
    void computeSimilarity2();
    double getWidthFactor(const std::vector<double> v1, const std::vector<double> v2); 
    void computeWidthFactor();
   
    //Functions for fuzzy clustering
    void computeCentroids(int k);
    void computeCentroids2(int k);
    void mergeClusters(int k);
    bool updateMembership(int k);
    inline bool canStop(){
      double epsilon = 0.01;
      TMatrixT<double> diffMatrix = fpsMembership - fpsNewMembership;
      double difference = diffMatrix.Norm1();
      return difference < epsilon;
    }


    std::vector<std::vector<unsigned int> > fclusters;               ///< collection of something
    std::vector<std::vector<double> >       fps;                     ///< the collection of points we are working on     
    std::vector<unsigned int>               fpointId_to_clusterId;   ///< mapping point_id -> clusterId     
    std::vector<std::vector<double> >       fsim;                    ///<
    std::vector<std::vector<double> >       fsim2;            	     ///<
    std::vector<std::vector<double> >       fsim3;            	     ///<
    double fMaxWidth;

    RTree fRTree;
    //std::vector< dbsPoint > fRect;
   
    //Matrices for Ben's fuzzy cluster
    TMatrixT<double>                         fpsMat;




  private:
    
    // eps radius
    // Two points are neighbors if the distance 
    // between them does not exceed threshold value.
    double fEps;
    double fEps2;
    // The fuzzyness factor needed for fuzzy clustering, commonly known as just m
    double fFuzzifier;
    // The maximum number of clusters to try, needed for fuzzy clustering
    int nMaxClusters;
    // The maximum number of iterations to try, needed for fuzzy clustering
    int nIterations;
    // The maximum distance between clusters prior to merging, needed for fuzzy clustering
    int fMergeCutoff;
    //minimum number of points
    unsigned int fMinPts;
    // Which clustering to run
    unsigned int fClusterMethod;  ///< Which clustering method to use
    unsigned int fDistanceMetric; ///< Which distance metric to use
      
    // noise vector
    static const unsigned int kNO_CLUSTER    = UINT_MAX;
    static const unsigned int kNOISE_CLUSTER = UINT_MAX-1;
    std::vector<bool>      fnoise;	
    std::vector<bool>      fvisited;					     
    std::vector<double>    fWirePitch;     ///< the pitch of the wires in each plane
    std::set<unsigned int> fBadChannels;   ///< set of bad channels in this detector
    std::vector<unsigned int> fBadWireSum; ///< running total of bad channels. Used for fast intervening 
                                           ///< dead wire counting ala fBadChannelSum[m]-fBadChannelSum[n]. 

    //Matrices for Ben's fuzzy cluster
    TMatrixT<double>                         fpsMembership;
    TMatrixT<double>                         fpsNewMembership;
    TMatrixT<double>                         fpsMembershipFinal;
    TMatrixT<double>                         fpsCentroids;

    // Object used for Hough transforms
    HoughClusAlg fHCAlg;        
    EndPointClusAlg fEPCAlg;        ///< object that does the Hough Transform


    // Helper routined for run_dbscan_cluster() names and
    // responsibilities taken directly from the paper
    bool ExpandCluster(unsigned int point /* to be added */, 
		       unsigned int clusterID /* which is being expanded */);
    std::set<unsigned int> RegionQuery(unsigned int point);
    // Helper for the accelerated run_FN_cluster()
    std::vector<unsigned int> RegionQuery_vector(unsigned int point);


  }; // class fuzzyClusterAlg
} // namespace

#endif // ifndef fuzzyClusterALG_H