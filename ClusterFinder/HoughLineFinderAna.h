#ifndef HOUGHLINEFINDERANA_H
#define HOUGHLINEFINDERANA_H

#include "ClusterFinder.h"
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
         
    void analyze(edm::Event& , edm::EventSetup const&);
    void beginJob(edm::EventSetup const&);
     
  private:

    TH1F* fGenerated[6];  ///< Spectra as generated. Unused, it'd seem.

    std::string fHoughModuleLabel;
    std::string fDigitModuleLabel;
    std::string fHitModuleLabel;
      TTree* tree;
      Int_t m_run;          // Run number
      Double_t m_run_timestamp;          // Run number
      Int_t m_event;        // Event number
      Int_t m_plane;        // Plane number
      Int_t m_clusterid;    // Cluster ID
      Int_t m_wirespan;    // Wire spanned by track
      Int_t m_sizeClusterZ;  //Number of clusters
      Int_t m_sizeHitZ;      //Number of Hits
      Int_t m_hitidZ;
      Int_t m_mipZ;
      Int_t m_drifttimeZ;
      Int_t m_widthZ;
      Int_t m_upadcZ;
      
  };
  
  
} // end namespace cluster



#endif // HoughLineFinderAna_H
