////////////////////////////////////////////////////////////////////////
/// \file  ShowerFinder.cxx
/// \brief Reconstruct Showers 
///
/// \version $Id: ShowerFinder.cxx,v 0.1 2010/10/12 06:08:30 brossi strauss $
/// \author brossi@lhep.unibe.ch
/// \author strauss@lhep.unibe.ch
////////////////////////////////////////////////////////////////////////
// This class solves the following problem:
//
// Autoitc shower reconstruction and distinguish between electron and gamma

// Framework includes
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"

#include "ShowerFinder/ShowerReco.h"


namespace shwf {

  DEFINE_FWK_MODULE(Shower);

}
