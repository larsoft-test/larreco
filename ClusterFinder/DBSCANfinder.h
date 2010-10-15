/////////////////////////////////////////////////////////////////
//  DBSCANfinder.h
//  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////
//#include "ClusterFinder.h"
#include <vector>
#include <cmath>
#include <iostream>
#include "FWCore/Framework/interface/EDProducer.h"
#include "RecoBase/Cluster.h"
#include "TH2F.h"

namespace edm { 
class Event;
class ParameterSet;}


class TH1F;
namespace cluster{
  
  //The line below is for my testing purposes only: 
  //void randomInit(std::vector<std::vector<double> > & ps, int dims = //5,unsigned int num_points = 10);
 
  //--------------------------------------------------------------- 
  class DBcluster : public edm::EDProducer
  {
  public:
    explicit DBcluster(edm::ParameterSet const& pset); 
    ~DBcluster();
    void produce(edm::Event& evt, edm::EventSetup const&);
   void beginJob(edm::EventSetup const&);
    
  private:
  
  edm::PtrVector<recob::Hit> allhits;
  edm::PtrVector<recob::Hit> clusterHits;
    // double fEps, fEps2;
//     int fMinPts;
     
 TH1F *fhitwidth;
 TH1F *fhitwidth_ind_test;  
 TH1F *fhitwidth_coll_test;  
    //the collection of clusters
    
     TH1F* fNoParticles_pdg;
    TH1F* fNoParticles_trackid; 
    TH1F* fNoParticles_trackid_mother;
    TH1F* fNoParticles_trackid_per_event;  
    TH1F* fNoParticles_pdg_per_event;
    TH1F* fCl_for_Muon;
   /*  TH1F* fCl_for_Electron; */
/*     TH1F* fCl_for_Positron; */
/*     TH1F* fCl_for_Pion_111; */
/*     TH1F* fCl_for_Pion_211; */
/*     TH1F* fCl_for_Pion_m211; */
/*     TH1F* fCl_for_Proton; */
    TH1F* fNoClustersInEvent;
    TH1F* fPercentNoise;
    TH1F* fno_of_clusters_per_track;
    TH1F* fPercent_lost_muon_hits;
    TH1F* fPercent_lost_electron_hits;
    TH1F* fPercent_lost_positron_hits;
    TH1F* fPercent_lost_111_hits;
    TH1F* fPercent_lost_211_hits;
    TH1F* fPercent_lost_m211_hits;
    TH1F* fPercent_lost_2212_hits;
    TH1F* fPercent_lost_2112_hits;

    TH1F* fPercent_lost_muon_energy;
    TH1F* fPercent_lost_electron_energy;
    TH1F* fPercent_lost_positron_energy;
    TH1F* fPercent_lost_111_energy;
    TH1F* fPercent_lost_211_energy;
    TH1F* fPercent_lost_m211_energy;
    TH1F* fPercent_lost_2212_energy;
    TH1F* fPercent_lost_2112_energy;
    TH1F* fEnergy;
    TH2F* fbrian_in;
    TH2F* fbrian_coll;
    TH1F* fhitwidth_;
    TH1F* fhitwidth_0;
    TH1F* fhitwidth_1;

  protected:
std::string fhitsModuleLabel;
    double fEps, fEps2;
    int fMinPts;
    
  };
  //--------------------------------------------------------------
  class DBScan {
  public:
    std::vector<std::vector<unsigned int> > _clusters;
    
    
    
    DBScan ( std::vector<std::vector<double> > & ps, double eps,double eps2, int minPts) : 
      _ps(ps), _eps(eps),_eps2(eps2), _minPts(minPts)
      {
	_pointId_to_clusterId.resize(_ps.size(), 0);
	_noise.resize(ps.size(), false);
	_visited.resize(ps.size(), false);
      };
      
      
      ~DBScan(){};
      

     double getSimilarity(const std::vector<double> v1, const std::vector<double> v2); 
     std::vector<unsigned int> findNeighbors( unsigned int pid, double threshold, double threshold2);
      void computeSimilarity();
      void run_cluster() ;
      
      double getSimilarity2(const std::vector<double> v1, const std::vector<double> v2); 
       void computeSimilarity2();
       double getWidthFactor(const std::vector<double> v1, const std::vector<double> v2); 
       void computeWidthFactor();
      // protected:
      
      
      // the collection of points we are working on
     
        std::vector<std::vector<double> > _ps;
      
      // mapping point_id -> clusterId
       std::vector<unsigned int> _pointId_to_clusterId; 
      
      // the collection of clusters
      // std::vector<std::vector<unsigned int>> _clusters;
      
      // simarity_matrix
      
      std::vector<std::vector<double> >_sim;
      std::vector<std::vector<double> >_sim2;
      std::vector<std::vector<double> >_sim3;
       	
      friend std::ostream& operator << (std::ostream& o, const DBScan& c);
      friend 	std::ostream& operator<<(std::ostream& o,const std::vector<double> & p);
      
  private:
      
      // eps radiuus
      // Two points are neighbors if the distance 
      // between them does not exceed threshold value.
      double _eps, _eps2;
      
      //minimum number of points
      unsigned int _minPts;
      
      // noise vector
      std::vector<bool> _noise;
      
      // noise vector
      std::vector<bool> _visited;
     
      
  };
}
