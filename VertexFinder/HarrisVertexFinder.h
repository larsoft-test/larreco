////////////////////////////////////////////////////////////////////////
/// \file  DriftElectrons.h
/// \brief Module to find vertices
///
/// \author  joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////////

#ifndef HarrisVertexFinder_H
#define HarrisVertexFinder_H

#include "FWCore/Framework/interface/EDProducer.h"
#include "TMath.h"
#include <vector>
#include <string>

namespace vertex {
   
 class HarrisVertexFinder :  public edm::EDProducer {
    
  public:
    
    explicit HarrisVertexFinder(edm::ParameterSet const& pset); 
    ~HarrisVertexFinder();        

    void produce(edm::Event& evt, edm::EventSetup const&);
    
  private:

    double Gaussian(int x, int y, double sigma);
    double GaussianDerivativeX(int x, int y);
    double GaussianDerivativeY(int x, int y);
  
    std::string fDBScanModuleLabel;
    int         fTimeBins;
    int         fMaxCorners;
    double      fGsigma;
    int         fWindow;
    double      fThreshold;
    int         fSaveVertexMap;
  };
    
}



#endif // HarrisVertexFinder_H
