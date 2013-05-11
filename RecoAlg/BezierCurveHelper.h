////////////////////////////////////////////////////////////////////////
///
/// \file   BezierCurveHelper
///
/// \brief  Service for finding 3D track seeds
///
/// \author B J P Jones
///
/// Helper object for interpolating betweeen bezier
/// curve points, specified by seed objects
///


#ifndef BEZIERCURVEHELPER_H
#define BEZIERCURVEHELPER_H

#include <vector>

#include "RecoBase/Seed.h"

#include "fhiclcpp/ParameterSet.h"

#include "TVector3.h"

namespace trkf {



  class BezierCurveHelper {
  public:

    // Constructor.
    BezierCurveHelper();
    explicit BezierCurveHelper(int fCurveRes);

    // Destructor.
    ~BezierCurveHelper();

    // Update configuration parameters.
    void reconfigure(const fhicl::ParameterSet& pset);

    std::vector<TVector3> GetBezierPoints(recob::Seed * s1, recob::Seed * s2, int N=100);
    double   GetSegmentLength(recob::Seed * s1, recob::Seed * s2);
    void     GetBezierPointXYZ(recob::Seed * s1, recob::Seed * s2, float t, double * xyz);
    TVector3 GetBezierPoint(recob::Seed * s1, recob::Seed * s2, float t);
    void     GetDirectionScales(double * Pt1, double * Pt2, double * Dir1, double * Dir2, double *Scales);


    void SetCurveResolution(int CurveRes) {fCurveResolution=CurveRes;}
    int  GetCurveResolution()             {return fCurveResolution;}

  private:
    int fCurveResolution;


  };



}

#endif