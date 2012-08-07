////////////////////////////////////////////////////////////////////////
/// \file  EndPointClusAlg.h
/// \brief algorithm to find 2D endpoints
///
/// \author  joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////////

#ifndef ENDPOINTALG_H
#define ENDPOINTALG_H

#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "TMath.h"
#include <vector>
#include <string>

namespace recob { 
  class Cluster;
  class EndPoint2D; 
  class Hit;
}

namespace cluster {
   
  ///Algorithm to find 2D end points
 class EndPointClusAlg {
    
  public:
    
    explicit EndPointClusAlg(fhicl::ParameterSet const& pset); 
    virtual ~EndPointClusAlg();        

    void   reconfigure(fhicl::ParameterSet const& pset);

    size_t EndPoint(art::PtrVector<recob::Cluster>                 & clusIn, 
		    std::vector<recob::EndPoint2D>                 & vtxcol,
		    std::vector< art::PtrVector<recob::Hit> >      & vtxHitsOut,
		    art::Event                                const& evt,
		    std::string                               const& label);
    
    size_t EndPoint(art::PtrVector<recob::Hit>& hits, std::vector<unsigned int> *corners ); 
    
  private:

    double Gaussian(int x, int y, double sigma);
    double GaussianDerivativeX(int x, int y);
    double GaussianDerivativeY(int x, int y);
    void VSSaveBMPFile(const char *fileName, unsigned char *pix, int dx, int dy);

    
    int          fTimeBins;
    int          fMaxCorners;
    double       fGsigma;
    int          fWindow;
    double       fThreshold;
    int          fSaveVertexMap;
  };
    
}



#endif // ENDPOINTALG_H
