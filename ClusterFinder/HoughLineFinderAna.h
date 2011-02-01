#ifndef HOUGHLINEFINDERANA_H
#define HOUGHLINEFINDERANA_H

#include "TMath.h"
#include <vector>
#include <string>

#include "art/Framework/Core/EDAnalyzer.h"

class TH1F;
class TTree;
namespace cluster {
   
  class HoughLineFinderAna : public art::EDAnalyzer {
    
  public:
    
    explicit HoughLineFinderAna(fhicl::ParameterSet const& pset); 
    ~HoughLineFinderAna();
         
    void analyze(const art::Event&);
    void beginJob();
     
  private:

    std::string fHoughModuleLabel;
    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fDBScanModuleLabel;       
    TTree* ftree;
    int fm_run;          // Run number
    unsigned long int fm_run_timestamp;          // Run number
    int fm_event;        // Event number
    int fm_plane;        // Plane number
    int fm_dbsize;
    int fm_clusterid;    // Cluster ID
    int fm_wirespan;    // Wire spanned by track
    int fm_sizeClusterZ;  //Number of clusters
    int fm_sizeHitZ;      //Number of Hits
    double fm_clusterslope;
    double fm_clusterintercept;
    int *fm_wireZ;
    int *fm_hitidZ;
    double *fm_mipZ;
    double *fm_drifttimeZ;
    double *fm_widthZ;
    double *fm_upadcZ;
      
  };
  
  
} // end namespace cluster



#endif // HoughLineFinderAna_H
