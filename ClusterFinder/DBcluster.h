/////////////////////////////////////////////////////////////////
//  DBSCANfinder.h
//  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////

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
    
     
 TH1F *fhitwidth;
 TH1F *fhitwidth_ind_test;  
 TH1F *fhitwidth_coll_test;  
  
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
      
      
      std::vector<std::vector<double> >_sim;
      std::vector<std::vector<double> >_sim2;
      std::vector<std::vector<double> >_sim3;
       	
      friend std::ostream& operator << (std::ostream& o, const DBScan& c);
      friend 	std::ostream& operator<<(std::ostream& o,const std::vector<double> & p);
      
  private:
      
      // eps radius
      // Two points are neighbors if the distance 
      // between them does not exceed threshold value.
      double _eps, _eps2;
      
      //minimum number of points
      unsigned int _minPts;
      
      // noise vector
      std::vector<bool> _noise;
      std::vector<bool> _visited;
     
      
  };
}
