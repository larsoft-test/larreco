#include "Utilities/AssociationUtil.h"
#include "RecoBase/Track.h"

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
  
  class BezierTrack: public recob::Track
  {
  public:
    BezierTrack();
    BezierTrack(recob::Track btb);
    BezierTrack(std::vector<recob::Seed*> );
     
    ~BezierTrack();

    double GetLength()                                  const;
    double GetRMSCurvature()                            const;

    TVector3 GetTrackPointV     (  double s )           const;
    TVector3 GetTrackDirectionV (  double s )           const;
    double   GetCurvature(double s)                     const;
 
    int NSegments() { return NumberTrajectoryPoints();} 

    void CalculateSegments();    
    
    void   GetProjectedPointUVWX( double s, double* uvw, double * x,  int c, int t ) const;  
    void   GetProjectedPointUVWT( double s, double* uvw, double * ticks, int c, int t ) const;  

    void GetClosestApproach( recob::Hit* hit,       double &s,  double& Distance) const;
    void GetClosestApproach( art::Ptr<recob::Hit> hit,       double &s,  double& Distance) const;
    void GetClosestApproach( recob::SpacePoint* sp, double &s,  double& Distance) const;
    void GetClosestApproach( TVector3 vec,          double &s,  double& Distance) const;
    

    void   GetTrackPoint    (  double s, double* xyz )           const;
    void   GetTrackDirection(  double s, double* xyz )           const;


  private:
    void FillSeedVector();
    
    std::vector<double>       fSegmentLength;
    std::vector<double>       fCumulativeLength;

    std::vector<recob::Seed* > fSeedCollection;
    
    double  fTrackLength;
    int     fBezierResolution;
    
  };
}

#endif
