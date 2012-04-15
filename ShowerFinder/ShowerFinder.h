////////////////////////////////////////////////////////////////////////
/// \file  ShowerFinder.h
/// \brief
///
///
/// \version $Id: SingleGen.cxx,v 1.4 2010/03/29 09:54:01 brebel Exp $
/// \author  roxanne
////////////////////////////////////////////////////////////////////////
#ifndef SHOWERFINDER_H
#define SHOWERFINDER_H

#include <vector>
#include <string>

#include "art/Framework/Core/EDProducer.h" 

///shower finding
namespace shwf {
   
  class ShowerFinder : public art::EDProducer {
    
  public:
    
    explicit ShowerFinder(fhicl::ParameterSet const&); 
    virtual ~ShowerFinder();
         
    void reconfigure(fhicl::ParameterSet const& p);
    void produce(art::Event& evt); 

  private:

    std::string fVertexModuleLabel;         ///< label of module finding 2D endpoint
    std::string fClusterModuleLabel;        ///< label of module finding clusters
    std::string fHoughLineModuleLabel;      ///< label of module finding hough line
    std::string fVertexStrengthModuleLabel; ///< label of module finding 2D endpoint 
    double      fRcone;                     ///< radious of cone for method
    double      fLcone;                     ///< length of the cone

  protected: 
    
    
  }; // class ShowerFinder
  
  
}

#endif // SHOWERFINDER_H
