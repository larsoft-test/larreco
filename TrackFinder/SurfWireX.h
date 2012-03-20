////////////////////////////////////////////////////////////////////////
///
/// \file   SurfWireX.h
///
/// \brief  Planar surface defined by readout wire and x-axis.
///
/// \author H. Greenlee 
///
/// This class derives from SurfYZPlane.  This class does not add any
/// new members, but has a constructor that allows construction from
/// a readout channel number.
///
////////////////////////////////////////////////////////////////////////

#ifndef SURFWIREX_H
#define SURFWIREX_H

#include "TrackFinder/SurfYZPlane.h"

namespace trkf {

  class SurfWireX : public SurfYZPlane
  {
  public:

    /// Constructor.
    SurfWireX(unsigned int channel);

    /// Destructor.
    virtual ~SurfWireX();
  };
}

#endif
