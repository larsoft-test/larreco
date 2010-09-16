////////////////////////////////////////////////////////////////////////
//
// FFTHitFinder class
//
// pagebri3@msu.edu
//
//  This algorithm is designed to find hits on wires after deconvolution
//  with an average shape used as the input response.
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h" 

#include "HitFinder/FFTHitFinder.h"

namespace hit{

  DEFINE_FWK_MODULE(FFTHitFinder);

} // end of hit namespace
