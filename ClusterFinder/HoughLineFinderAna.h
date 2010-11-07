#ifndef HOUGHLINEFINDERANA_H
#define HOUGHLINEFINDERANA_H

#include "TMath.h"
#include <vector>
#include <string>

#include "FWCore/Framework/interface/EDAnalyzer.h"

class TH1F;
class TTree;
namespace cluster {
   
  class HoughLineFinderAna : public edm::EDAnalyzer {
    
  public:
    
    explicit HoughLineFinderAna(edm::ParameterSet const& pset); 
    ~HoughLineFinderAna();
         
    void analyze(const edm::Event& , edm::EventSetup const&);
    void beginJob(edm::EventSetup const&);
     
  private:

    std::string fHoughModuleLabel;
    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fDBScanModuleLabel;       
    TTree* ftree;
    Int_t fm_run;          // Run number
    ULong64_t fm_run_timestamp;          // Run number
    Int_t fm_event;        // Event number
    Int_t fm_plane;        // Plane number
    Int_t fm_dbsize;
    Int_t fm_clusterid;    // Cluster ID
    Int_t fm_wirespan;    // Wire spanned by track
    Int_t fm_sizeClusterZ;  //Number of clusters
    Int_t fm_sizeHitZ;      //Number of Hits
    Float_t fm_clusterslope;
    Float_t fm_clusterintercept;
    Int_t *fm_wireZ;
    Int_t *fm_hitidZ;
    Float_t *fm_mipZ;
    Float_t *fm_drifttimeZ;
    Float_t *fm_widthZ;
    Float_t *fm_upadcZ;
      
  };
  
  
} // end namespace cluster



#endif // HoughLineFinderAna_H
