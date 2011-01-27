/////////////////////////////////////////////////////////////////
//  DBSCANfinder.h
//  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////

#include <vector>
#include <cmath>
#include <iostream>
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"

class TH1F;

namespace recob { class Hit; }

namespace cluster{
   
  //--------------------------------------------------------------- 
  class DBScanService {
  public:
    
    
    DBScanService(edm::ParameterSet const& pset, edm::ActivityRegistry& reg);
    virtual ~DBScanService();
    
    void InitScan(edm::PtrVector<recob::Hit>& allhits);
    double getSimilarity(const std::vector<double> v1, const std::vector<double> v2); 
    std::vector<unsigned int> findNeighbors( unsigned int pid, double threshold, double threshold2);
    void computeSimilarity();
    void run_cluster();     
    double getSimilarity2(const std::vector<double> v1, const std::vector<double> v2); 
    void computeSimilarity2();
    double getWidthFactor(const std::vector<double> v1, const std::vector<double> v2); 
    void computeWidthFactor();
      

    std::vector<std::vector<unsigned int> > fclusters;               ///> collection of something
    std::vector<std::vector<double> >       fps;                     ///> the collection of points we are working on     
    std::vector<unsigned int>               fpointId_to_clusterId;   ///> mapping point_id -> clusterId     
    std::vector<std::vector<double> >       fsim;                    ///>
    std::vector<std::vector<double> >       fsim2;            	     ///>
    std::vector<std::vector<double> >       fsim3;            	     ///>
     
/*     friend std::ostream& operator <<(std::ostream& o, const DBScanService& c); */
/*     friend std::ostream& operator <<(std::ostream& o,const std::vector<double> & p); */
     
  private:
      
    // eps radius
    // Two points are neighbors if the distance 
    // between them does not exceed threshold value.
    double fEps;
    double fEps2;
     
    //minimum number of points
    unsigned int fMinPts;
      
    // noise vector
    std::vector<bool> fnoise;
    std::vector<bool> fvisited;
  };
} // namespace
