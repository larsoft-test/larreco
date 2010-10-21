#ifndef HOUGHLINEFINDERANA_H
#define HOUGHLINEFINDERANA_H

#include "TMath.h"
#include "TObject.h"
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
    std::string fHitModuleLabel;
      TTree* ftree;
      Int_t fm_run;          // Run number
      Double_t fm_run_timestamp;          // Run number
      Int_t fm_event;        // Event number
      Int_t fm_plane;        // Plane number
      Int_t fm_clusterid;    // Cluster ID
      Int_t fm_wirespan;    // Wire spanned by track
      Int_t fm_sizeClusterZ;  //Number of clusters
      Int_t fm_sizeHitZ;      //Number of Hits
      Int_t *fm_hitidZ;
      Double_t *fm_mipZ;
      Double_t *fm_drifttimeZ;
      Double_t *fm_widthZ;
      Double_t *fm_upadcZ;
      
  };
  
  
} // end namespace cluster



#endif // HoughLineFinderAna_H
