#include "RecoBase/BezierTrackBase.h"

#ifndef BEZIERTRACK_h
#define BEZIERTRACK_h

class TVector3;

namespace recob
{
  class Seed;
  class Hit;
  class SpacePoint;
}

namespace trkf {
  
  class BezierTrack: public recob::BezierTrackBase
  {
  public:
    BezierTrack();
    BezierTrack(recob::BezierTrackBase btb);
    BezierTrack(std::vector<recob::Seed*> );
     
    ~BezierTrack();

    double GetLength()  const { return fTrackLength;}
    void   GetTrackPoint(     double s, double* xyz )           const;

    void CalculateSegments();    
    
    void   GetProjectedPointUVWX( double s, double* uvw, double * x,  int c, int t ) const;  
    void   GetProjectedPointUVWT( double s, double* uvw, double * ticks, int c, int t ) const;  

    double GetClosestApproach( recob::Hit* hit,       double *s,  double* Distance);
    double GetClosestApproach( recob::SpacePoint* sp, double *s,  double* Distance);
    double GetClosestApproach( TVector3 vec,          double *s,  double* Distance);
    
    //  void   GetTrackGradient(  double s, double* xyz )           const;
    //  void   GetTrackCurvature( double s, double* xyz )           const;
    //  void   GetTicks(          double s, double* t)              const;
 
     
   
  private:
    void FillSeedVector();
    
    std::vector<double>       fSegmentLength;
    std::vector<double>       fCumulativeLength;

    std::vector<recob::Seed* > fSeedCollection;
    
    double fTrackLength;

    
    
  };
}

#endif
