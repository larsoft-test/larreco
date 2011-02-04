////////////////////////////////////////////////////////////////////////
/// \file  HarrisVertexFinder.h
/// \brief Module to find vertices
///
/// \author  joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////////

#ifndef VertexService_H
#define VertexService_H

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "TMath.h"
#include <vector>
#include <string>

///vertex reconstruction
namespace recob { class Cluster; }

namespace recob { class Vertex; }

namespace vertex {
   
 class VertexService {
    
  public:
    
    explicit VertexService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg); 
    virtual ~VertexService();        

    size_t Vertex(art::PtrVector<recob::Cluster>& clusIn, std::vector<recob::Vertex>& vtxcol);
    
  private:

    double Gaussian(int x, int y, double sigma);
    double GaussianDerivativeX(int x, int y);
    double GaussianDerivativeY(int x, int y);
  void VSSaveBMPFile(const char *fileName, unsigned char *pix, int dx, int dy);

    
    int         fTimeBins;
    int         fMaxCorners;
    double      fGsigma;
    int         fWindow;
    double      fThreshold;
    int         fSaveVertexMap;
  };
    
}



#endif // VertexService_H
