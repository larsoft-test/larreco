////////////////////////////////////////////////////////////////////////
//
// ScanFilter class:
// Tells the downstream rreconstruction to not process events that
// do not pass certain hand scan criteria
//
// joshua.spitz@yale.edu
//
////////////////////////////////////////////////////////////////////////

/// Framework includes
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h" 
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h" 
 
// LArSoft includes
#include "ScanFilter.h"

namespace filt {

  DEFINE_FWK_MODULE(ScanFilter);

} //namespace filt
