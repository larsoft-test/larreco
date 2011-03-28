////////////////////////////////////////////////////////////////////////
// $Id: KingaClusterAna.cxx,v 1.36 2010/09/15  kpartyka Exp $
//
// KingaCluster class
//
// 
//
////////////////////////////////////////////////////////////////////////
#ifndef KingaCluster_H
#define KingaCluster_H

#include "TMath.h"
#include <vector>
#include <string>
#include "art/Persistency/Common/PtrVector.h" 


#include "art/Framework/Core/EDProducer.h"

class TH1F;
class TH2F;
class TTree;
namespace recob { class Hit; }
namespace cluster {
   
  class KingaCluster : public art::EDProducer {
    
  public:
    
    explicit KingaCluster(fhicl::ParameterSet const& pset); 
    KingaCluster();
    ~KingaCluster();
         
    void produce(art::Event& evt);
    void beginJob();
    void AngularDistribution(int plane);
    void FindMax(int plane);
    //void FinalPeaks();
    void FitAngularDistributions();
    void FindClusters(int plane);
 
  
     
  private:
    std::string fDBScanModuleLabel;  
    std::string fGenieGenModuleLabel;  
    TH1F *fh_theta_ind;
    TH1F *fh_theta_coll;
    TH1F *fh_theta_ind_2D;
    TH1F *fh_theta_coll_2D;
    TH1F *fh_theta_coll_Area;
    TH1F *fh_theta_ind_Area;
     //std::vector<TH1F*> fh_theta;     /**Histo for the angular distribution theta of the shower*/
    art::PtrVector<recob::Hit> allhits;
    std::vector<int> maxBin;    //stores bin # of local maximum
    std::vector<int> MaxStartPoint; //bin no of the starting point of a peak
    std::vector<int> MaxEndPoint;  //bin no of the end point of a peak
    std::vector<int> MaxStartPointTheta; //theta value of the starting point of a peak
    std::vector<int> MaxEndPointTheta; //theta value of the end point of a peak
    std::vector<unsigned int> fwire_vertex;
    std::vector<unsigned int> ftime_vertex;
    std::vector<int> HitsWithClusterID;
    double ftimetick; //get from parameterset
double fdriftvelocity;  //get from paramtereset 9either k and V)
double fpi;
int fMC;
double MCvertex [3];

std::vector<double> maxBinValues;
std::vector<double> OriginalmaxBinValues;
std::vector<int> SortedMaxBin;
std::vector<int> FinalPeaks;
  protected:

    
  };
  
 // class AngCluster{
//   public:
//    AngCluster();
//    ~AngCluster();
//    void beginJob();
//    void AngularDistribution(art::PtrVector<recob::Hit>);
//   
//   private:
//   TH1F *fh_theta_ind;
//     TH1F *fh_theta_coll;
//     TH1F *fh_theta_ind_2D;
//     TH1F *fh_theta_coll_2D;
//   
//   
//   };
}



#endif // KingaCluster_H
