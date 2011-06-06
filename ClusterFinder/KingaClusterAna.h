////////////////////////////////////////////////////////////////////////
//
// 
// \author kinga.partyka@yale.edu
//
// 
////////////////////////////////////////////////////////////////////////
#ifndef KingaClusterAna_H
#define KingaClusterAna_H

#include "art/Framework/Core/EDAnalyzer.h"

#include <vector>
#include <string>


class TH1F;
class TH2F;
namespace recob { class Hit; }
///Cluster finding and building 
namespace cluster {

   
  class KingaClusterAna : public art::EDAnalyzer {

  public:
          
    explicit KingaClusterAna(fhicl::ParameterSet const& pset); 
    virtual ~KingaClusterAna();
 
    /// read access to event
    void analyze(const art::Event& evt);
    void beginJob();

  private:
  art::PtrVector<recob::Hit> allhits;
  std::vector<unsigned int> fwire_vertex;
  std::vector<unsigned int> ftime_vertex;
  double ftimetick;
  double fdriftvelocity; 
  std::string fKingaModuleLabel;
  std::string fGenieGenModuleLabel;
  std::string fLArG4ModuleLabel;
  std::string fHitsModuleLabel;
  std::string fClusterFinderModuleLabel;
   
    
    
    TH1F* Mu_theta;
    TH1F* Mu_phi;
    TH1F* Mu_phi_oneside;
    
    TH1F* pion_theta;
    TH1F* pion_phi;
    TH1F* pion_phi_oneside;
    
    TH1F* Energy_in_Sphere;
    TH1F* M_Delta_plus_plus;
    TH1F* M_Delta_plus_plus2;
   TH1F* M_Delta_plus_plus_Mother;
    TH1F* Number_protons;

TH1F* Ind_eng_rectangle;
TH1F* Coll_eng_rectangle;
TH1F* Ind_eng_rectangle2;
TH1F* Coll_eng_rectangle2;
TH1F* Ind_eng_rectangle3;
TH1F* Coll_eng_rectangle3;

TH1F* Vertex_x;
TH1F* Vertex_y;
TH1F* Vertex_z;





  }; // class KingaClusterAna

}

#endif 
