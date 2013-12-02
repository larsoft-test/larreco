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
#include "TVectorF.h"
#include "TVector.h"
#include "TH1.h"


class TH1F;

namespace recob { class Hit; }

namespace cluster{


  // This stores information about a showerlike cluster
  class showerCluster
    {
      public:
        int clusterNumber=-999999;
        std::vector<protoTrack> clusterProtoTracks;
        showerCluster (protoTrack protoTrackTemp)
        {
          clusterNumber=protoTrackTemp.clusterNumber;
          clusterProtoTracks.push_back(protoTrackTemp);
        }

        void addProtoTracks(std::vector<protoTrack> tracksToAdd){
          
          for(auto tracksToAddItr = tracksToAdd.begin(); tracksToAddItr != tracksToAdd.end(); tracksToAddItr++)
            tracksToAddItr->clusterNumber = clusterNumber;
          clusterProtoTracks.insert(clusterProtoTracks.end(),tracksToAdd.begin(),tracksToAdd.end());
        }
        
        void clearProtoTracks(){
          clusterProtoTracks.clear();
        }

    };

  // This stores information about a tracklike cluster
  class trackCluster
    {
      public:
        int clusterNumber=-999999;
        std::vector<protoTrack> clusterProtoTracks;
        trackCluster (protoTrack protoTrackTemp)
        {
          clusterNumber=protoTrackTemp.clusterNumber;
          clusterProtoTracks.push_back(protoTrackTemp);
        }

        void addProtoTracks(std::vector<protoTrack> tracksToAdd){

          for(auto tracksToAddItr = tracksToAdd.begin(); tracksToAddItr != tracksToAdd.end(); tracksToAddItr++)
            tracksToAddItr->clusterNumber = clusterNumber;
          clusterProtoTracks.insert(clusterProtoTracks.end(),tracksToAdd.begin(),tracksToAdd.end());
        }
        
        void clearProtoTracks(){
          clusterProtoTracks.clear();
        }

    };


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
      float epsilon = 0.01;
      TMatrixF diffMatrix = fpsMembership - fpsNewMembership;
      float difference = diffMatrix.Norm1();
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
    float  fFuzzyRemnantMergeCutoff;       ///< cut off on merging the fuzzy cluster remnants into the nearest shower or track 

    int    fDoTrackClusterMerge;           ///< Turn on cut on product of charge asymmetry and sin of angle between slopes of lines
    float  fTrackClusterMergeCutoff;          ///< Max distance between Hough lines before two lines are merged (muon tracks), 
    float  fChargeAsymAngleCut;            ///< Cut on product of charge asymmetry and sin of angle between slopes of lines
    float  fSigmaChargeAsymAngleCut;       ///< Cut on product of charge asymmetry and sin of angle between slopes of lines
  
    int    fDoShowerClusterMerge;          ///< Turns on shower Hough line merging (0-off, 1-on)
    float  fShowerClusterMergeAngle;       ///< Max angle between slopes before two lines are merged, for lines in shower line regions
    float  fShowerClusterMergeCutoff;    ///< Max distance between Hough lines before two lines are merged (electron showers),

    int    fDoShowerTrackClusterMerge;     ///< Turn on cut on product of charge asymmetry and sin of angle between slopes of lines
    float  fShowerTrackClusterMergeCutoff;    ///< Max distance between Hough lines before two lines are merged (electron showers),
    float  fShowerTrackClusterMergeAngle;  ///< Max angle between slopes before two lines are merged, for lines in shower line regions
    
    float  fShowerLikenessCut;             ///< Cut on shower likeness (larger the more shower like, smaller the less shower like)

    int    fMaxVertexLines;                ///< Max number of line end points allowed in a Hough line merge region for a merge to happen







    void mergeHoughLinesBySegment(unsigned int k,
        std::vector<protoTrack> *protoTracks, 
        float xyScale,
        int mergeStyle,
        float wire_dist,
        float tickToDist);

    bool mergeTrackClusters(unsigned int k,
        std::vector<trackCluster> *trackClusters, 
        float xyScale,
        float wire_dist,
        float tickToDist);

    bool mergeShowerClusters(unsigned int k,
        std::vector<showerCluster> *showerClusters, 
        float xyScale,
        float wire_dist,
        float tickToDist);
    
    bool mergeShowerTrackClusters(showerCluster *showerClusterI, 
        trackCluster *trackClusterJ, 
        float xyScale,
        float wire_dist,
        float tickToDist);

    //std::vector<lineSlope> linesFound;
    float HoughLineDistance(float p0MinLine1, 
        float p1MinLine1, 
        float p0MaxLine1, 
        float p1MaxLine1, 
        float p0MinLine2, 
        float p1MinLine2, 
        float p0MaxLine2, 
        float p1MaxLine2);
    bool   HoughLineIntersect(float x11,
        float  y11,
        float  x12,
        float  y12,
        float  x21,
        float  y21,
        float  x22,
        float  y22);
    float PointSegmentDistance(float px,
        float  py,
        float  x1,
        float  y1,
        float  x2,
        float  y2);

    float DistanceBetweenHits(
        art::Ptr<recob::Hit> hit0,
        art::Ptr<recob::Hit> hit1,
        float wire_dist,
        float tickToDist);

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














