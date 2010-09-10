///////////////////////////////////////////////////////
//
// ChannelFilter Class
//
//
//  pagebri3@msu.edu
//
///////////////////////////////////////////////////////
#include "Filters/inc/ChannelFilter.h"
#include "Geometry/inc/Geometry.h"

filter::ChannelFilter::ChannelFilter()
{
}

filter::ChannelFilter::~ChannelFilter()
{
}

bool filter::ChannelFilter::BadChannel(unsigned int channel) 
{
  //geo::Geometry * geom = geo::Geometry::Instance();  
  edm::Service<geo::Geometry> geom;
  if(geom->IsDetector("argoneut")) {
    
    switch(channel) {
    case 22:
    case 65:
    case 171:
    case 237:
    case 307:
    case 308:
    case 309:
    case 310:
    case 311:
    case 410:
      return true;
      break;
    default:
      return false;
      break;
    }
  }
  return false;
   
}

bool filter::ChannelFilter::NoisyChannel(unsigned int channel) 
{
  //geo::Geometry * geom = geo::Geometry::Instance();  
  edm::Service<geo::Geometry> geom;
  if(geom->IsDetector("argoneut")) {
     switch(channel) {
    case 31:
    case 41:
    case 108:
    case 120:
    case 121:
    case 124:
    case 392:
    case 399:
      return true;
      break;
    default:
      return false;
      break;
    }
  }
  return false;
}


