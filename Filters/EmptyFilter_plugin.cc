////////////////////////////////////////////////////////////////////////
//
// EmptyFilter class:
// Algorith to produce event files with the
// blank events removed using only hit information.
//
// pagebri3@msu.edu
//
////////////////////////////////////////////////////////////////////////

/// Framework includes
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h" 
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h" 
 
// LArSoft includes
#include "EmptyFilter.h"

namespace filt {

  DEFINE_FWK_MODULE(EmptyFilter);

} //namespace filt
