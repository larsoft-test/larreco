////////////////////////////////////////////////////////////////////////
/// \file  PrimaryVertexFinder.h
/// \brief Module to find vertices based on 2-d clusters
///
/// \tjyang@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef VertexFinder2D_H
#define VertexFinder2D_H

#include "art/Framework/Core/EDProducer.h"
#include "TMath.h"
#include <vector>
#include <string>
#include "RecoBase/recobase.h"

class TH1F;
class TH2F;

///vertex reconstruction
namespace vertex {
   
 class VertexFinder2D :  public art::EDProducer {
    
  public:
    
    explicit VertexFinder2D(fhicl::ParameterSet const& pset); 
    virtual ~VertexFinder2D();        
    void beginJob();
    void reconfigure(fhicl::ParameterSet p);

    
    void produce(art::Event& evt);

  private:

  
    std::string fClusterModuleLabel;

  };
    
}



#endif // VertexFinder2D_H
