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
    BezierTrack(std::vector<TVector3> Pos, 
		std::vector<TVector3> Dir, 
		std::vector<std::vector<double> > dQdx);
    BezierTrack(std::vector<recob::Seed*> );
     
    ~BezierTrack();

    int NSegments()                                     const;

    double GetLength()                                  const;
    double GetRMSCurvature()                            const;

    double GetTotalCharge( unsigned int View )          const;
    double GetViewdQdx(    unsigned int View )          const; 

    TVector3 GetTrackPointV     (  double s )           const;
    TVector3 GetTrackDirectionV (  double s )           const;
    void   GetTrackPoint    (  double s, double* xyz )  const;
    void   GetTrackDirection(  double s, double* xyz )  const;  
 
    double GetCurvature(double s)                       const;
    double GetdQdx(double s, unsigned int View)         const ;
    
    void   GetProjectedPointUVWX( double s, double* uvw, double * x,  int c, int t )    const;  
    void   GetProjectedPointUVWT( double s, double* uvw, double * ticks, int c, int t ) const;  

    void   GetClosestApproach( recob::Hit* hit,            double &s,  double& Distance) const;
    void   GetClosestApproach( art::Ptr<recob::Hit> hit,   double &s,  double& Distance) const;
    void   GetClosestApproach( recob::SpacePoint* sp,      double &s,  double& Distance) const;
    void   GetClosestApproach( TVector3 vec,               double &s,  double& Distance) const;

    
    void   CalculatedQdx(art::PtrVector<recob::Hit>);   
    void   FillMySpacePoints(int N);
  
    recob::Track GetBaseTrack();

    
  private:
 
    int WhichSegment(double S) const;
    void CalculateSegments();    
   
  
    void FillSeedVector();
    
    std::vector<double>       fSegmentLength;
    std::vector<double>       fCumulativeLength;

    std::vector<recob::Seed* > fSeedCollection;
    
    double  fTrackLength;
    int     fBezierResolution;
    
  };
}

#endif
