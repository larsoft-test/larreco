////////////////////////////////////////////////////////////////////////
//
// ChannelFilter class:
//
// This class provides methods for returning the condition of
// a wire, as to whether it is bad and is to be ignored or perhaps
// if it has some known problem.  This allows the removal of detector
// specific code from a few places in LArSoft.  Right now is is only 
// implemented for Argoneut.
//  
//
// pagebri3@msu.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef CHANNELFILTER_H
#define CHANNELFILTER_H

// Framework includes
//////////////////////////////////////////////////////
// Of course these are all different
//////////////////////////////////////////////////////
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/ParameterSet/interface/ParameterSet.h" 
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h" 
#include "DataFormats/Common/interface/Handle.h" 
#include "DataFormats/Common/interface/Ptr.h" 
#include "DataFormats/Common/interface/PtrVector.h" 
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "FWCore/Services/interface/TFileService.h" 
#include "FWCore/Framework/interface/TFileDirectory.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h" 
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h" 
#include "Geometry/inc/Geometry.h"
#include "Geometry/inc/PlaneGeo.h"

namespace filter {

  class ChannelFilter {

  public:

    ChannelFilter();
    ~ChannelFilter();

    bool BadChannel(unsigned int channel);
    bool NoisyChannel(unsigned int channel);

  private:

  }; //class ChannelFilter
}
#endif // CHANNELFILTER_H

    

