////////////////////////////////////////////////////////////////////////
//
// ClusterFinderAna base class:
// Supply the basic methods for cluster finding algorithms.
// Users must define cluster finding specifics in inheriting classes.
//
// kinga.partyka@yale.edu
//
// joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////////
#ifndef CLUSTERFINDERANA_H
#define CLUSTERFINDERANA_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RawData/DAQHeader.h"

#include <vector>
#include <string>


class TH1F;
class TH2F;
///Cluster finding and building 
namespace cluster {

   
  class ClusterFinderAna : public edm::EDAnalyzer {

  public:
          
    explicit ClusterFinderAna(edm::ParameterSet const& pset); 
    virtual ~ClusterFinderAna();
 
    /// read access to event
      void analyze(const edm::Event& evt,  edm::EventSetup const&);
    // analyze() Must go into new module. EC, 5-Oct-2010.
    //    void analyze(const edm::Event& evt);
      void beginJob(edm::EventSetup const&);

  private:
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
    
    std::string fDigitModuleLabel;
    std::string fHitModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fClusterFinderModuleLabel;
    std::string fDetSimModuleLabel;
    std::string fGenieGenModuleLabel;

    
      	 
  }; // class ClusterFinderAna

}

#endif // CLUSTERFINDERANA_H
