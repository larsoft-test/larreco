////////////////////////////////////////////////////////////////////////
/// \file  HarrisVertexFinder.h
/// \brief Module to find vertices
///
/// \author  joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////////

#ifndef VertexModule_H
#define VertexModule_H

#include "art/Framework/Core/EDProducer.h"
#include "TMath.h"
#include <vector>
#include <string>

///vertex reconstruction
namespace vertex {
   
 class VertexModule :  public art::EDProducer {
    
  public:
    
    explicit VertexModule(fhicl::ParameterSet const& pset); 
    virtual ~VertexModule();        

    void produce(art::Event& evt);
    
  private:

    std::string fDBScanModuleLabel;
   
  };
    
}



#endif // VertexModule_H
