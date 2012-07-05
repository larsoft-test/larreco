////////////////////////////////////////////////////////////////////////
/// \file  AxisRecoAlg.h
/// \brief algorithm to reconstruct 3D angle from 2d angles
///
/// \author  andrzej.szelc@yale.edu
////////////////////////////////////////////////////////////////////////

#ifndef AXISRECOALG_H
#define AXISRECOALG_H

#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "TMath.h"
#include <vector>
#include <string>

#include "Geometry/geo.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"

namespace shwf {
   
  ///Algorithm to find 2D end points
 class AxisRecoAlg {
    
  public:
    
    explicit AxisRecoAlg(); 
    virtual ~AxisRecoAlg();        

    void   reconfigure(fhicl::ParameterSet const& pset);

    int Get3DaxisN(int iplane0,int iplane1,double omega0, double omega1,double &phi,double &theta);
    
    double Get2Dangle(double wire,double time);
    double Get2Dangle(double wireend,double wirestart,double timeend,double timestart);
    
  private:

    art::ServiceHandle<geo::Geometry> geom; 
    art::ServiceHandle<util::DetectorProperties> detp; 
    art::ServiceHandle<util::LArProperties> larp; 
    
  };
    
}



#endif // AXISRECOALG_H
