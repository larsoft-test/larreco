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
   std::string fKingaModuleLabel;
    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fClusterFinderModuleLabel;
    std::string fCalDataModuleLabel;
    std::string fGenieGenModuleLabel;

  }; // class KingaClusterAna

}

#endif 
