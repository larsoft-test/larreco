////////////////////////////////////////////////////////////////////////
//
// GausHitFinder class
//
// jaasaadi@syr.edu
//
//  Hit Finder algorithm based on FFTHitFinder to find
//  hits on wires after deconvolution with an average shape 
//  used as the input response. Using "seed" hits to better calculate
//  errors associated with the Hit and calculate the mulitplicity 
//  of the the found Hit
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 

#include "HitFinder/GausHitFinder.h"

namespace hit{

  DEFINE_ART_MODULE(GausHitFinder);

} // end of hit namespace
