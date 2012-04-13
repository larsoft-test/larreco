#include "TrackFinder/BezierTrack.h"
#include "TrackFinder/BezierCurveHelper.h"
#include "Utilities/DetectorProperties.h"
#include "Geometry/geo.h"
#include "TVector3.h"
#include "cetlib/exception.h"

namespace trkf {

  //----------------------------------------------------------------------
  // Constructor from Base object (for analysis)
  //
  BezierTrack::BezierTrack(recob::BezierTrackBase btb):
    recob::BezierTrackBase(btb)
  {
    fBezierResolution=1000;
    
    CalculateSegments();
  }


  //----------------------------------------------------------------------
  // Constructor from seed vector (for production)
  //
  BezierTrack::BezierTrack(std::vector<recob::Seed*> SeedCol):
    recob::BezierTrackBase()
  {
    fSeedCollection = SeedCol;
    CalculateSegments();
    fBezierResolution=1000;
  }



  //----------------------------------------------------------------------
  // Default constructor
  //
  BezierTrack::~BezierTrack()
  {
  }


  //----------------------------------------------------------------------
  // Given track points, fill the seed vector.  The seeds are then used 
  // for geomtry calculations / bezier fitting
  //
  void BezierTrack::FillSeedVector()
  {
    int NSeg=NSegments();
    
    std::cout<<"Filling seed vector of bezier track with n= "<<NSeg<<std::endl;
  
    double Pt[3], Dir[3];
    for(int i=0; i!=NSeg; i++)
      {


	Pt[0]=fPtX.at(i);
	Pt[1]=fPtY.at(i);
	Pt[2]=fPtZ.at(i);

	Dir[0]=fDirX.at(i);
	Dir[1]=fDirY.at(i);
	Dir[2]=fDirZ.at(i);
	
	fSeedCollection.push_back(new recob::Seed(Pt,Dir));
      }
    
  }
    


  //----------------------------------------------------------------------
  // Calculate the lengths of each bezier segment and fill useful 
  // calculational data members
  //

  void BezierTrack::CalculateSegments()
  {
    if(fSeedCollection.size()==0)
      { 
	if(NSegments()!=0) FillSeedVector();
	else 
	  {
	    throw cet::exception("no points in track")
	      <<"CalculateSegments method of Bezier track called with no"
	      <<" track information loaded.  You must fill track with "
	      <<"poisition and direction data before calling this method."
	      <<std::endl;
	    return;
	  }
      }
    std::cout<<"Loading bezier track with ID "<<  ID()<<std::endl;
    BezierCurveHelper bhlp(100);
    
    int Segments = NSegments();
    
    fTrackLength=0;
    bool FirstSeg=true;
    for(int i=0; i!=Segments; i++)
      {
	if(!FirstSeg)
	  {
	    float SegmentLength=bhlp.GetSegmentLength(fSeedCollection.at(i-1),fSeedCollection.at(i));
	  
	    fTrackLength+=SegmentLength;
	    fSegmentLength.push_back(SegmentLength);
	  }
	FirstSeg=false;
	fCumulativeLength.push_back(fTrackLength);
      }      
    std::cout<<"Total track length : " <<fTrackLength<<std::endl;
  }


  //----------------------------------------------------------------------
  // Find the point which is fraction s along the track and return
  // as a double[3]
  //
  void BezierTrack::GetTrackPoint(double s, double * xyz) const
  {
    if((s>1)||(s<0))
      {
	throw cet::exception("track point out of range")<<" s = "<<s <<" out of range \n";
      }
    else
      {
	BezierCurveHelper bhlp;
        for(unsigned int i=1; i!=fCumulativeLength.size(); i++)
          {
            if(  (   (fCumulativeLength.at(i-1) / fTrackLength) <=  s)
                 &&( (fCumulativeLength.at(i)   / fTrackLength) > s))
              {
		
                double locals = (s * fTrackLength - fCumulativeLength[i-1])/fSegmentLength[i-1];
		bhlp.GetBezierPointXYZ(fSeedCollection.at(i-1),fSeedCollection.at(i),locals, xyz);
              }

          }
      }
  }
  




  //----------------------------------------------------------------------
  //  Find the point which is fraction s along the track and get its 
  //   projected point in the wire view, in system uvwx
  //

  void BezierTrack::GetProjectedPointUVWX(double s, double *uvw, double*x, int t=0, int c=0) const
  {

    art::ServiceHandle<geo::Geometry>   geo;

    double  xyz[3];

    // Get point in 3D space
    GetTrackPoint(s , xyz);

    unsigned int c1,t1,p1,wirepoint;
    int NPlanes=geo->Cryostat(c).TPC(t).Nplanes();    

    for(int p=0; p!=NPlanes; p++)
      {
	unsigned int channelpoint = geo->NearestChannel(xyz,p,t,c);
	geo->ChannelToWire(channelpoint,c1,t1,p1,wirepoint);
	uvw[p]=wirepoint;
      }
    x[0]=xyz[0];
  }






  //----------------------------------------------------------------------
  //  Find the point which is fraction s along the track and get its 
  //   projected point in the wire view, in system uvwt
  //

  void BezierTrack::GetProjectedPointUVWT(double s, double *uvw, double*ticks, int t=0, int c=0) const
  {
    art::ServiceHandle<util::DetectorProperties> det;
    art::ServiceHandle<geo::Geometry>            geo;

    double  xyz[3];

    // Get point in 3D space
    GetTrackPoint(s , xyz);

    unsigned int c1,t1,p1,wirepoint;
    int NPlanes=geo->Cryostat(c).TPC(t).Nplanes();

    for(int p=0; p!=NPlanes; p++)
      {
	unsigned int channelpoint = geo->NearestChannel(xyz,p,t,c);
	geo->ChannelToWire(channelpoint,c1,t1,p1,wirepoint);
	uvw[p]=wirepoint;
	ticks[p]=det->ConvertXToTicks(xyz[0],p,t,c);
      }
  }



  //----------------------------------------------------------------------
  //  Calculate the closest approach of this track to a given hit
  //   and also the point where this occurs
  //

  void BezierTrack::GetClosestApproach( recob::Hit* hit,       double& s,  double& Distance) const 
  {
    art::ServiceHandle<util::DetectorProperties> det;
    art::ServiceHandle<geo::Geometry>            geo;

    unsigned int c1, t1, p1, w1;
    
    double xyzend1[3], xyzend2[3];
    
    double channel = hit->Channel();
    geo->ChannelToWire(channel,c1,t1,p1,w1);
    geo->WireEndPoints(c1,t1,p1,w1,xyzend1,xyzend2);
    
    xyzend1[0] = xyzend2[0] = det->ConvertTicksToX(hit->PeakTime(),p1,t1,c1);
    
    double iS, xyz[3], MinDistanceToPoint=10000, MinS=0;

    for(int i=0; i!=fBezierResolution; ++i)
      {
	iS=float(i)/fBezierResolution;
	GetTrackPoint(iS, xyz);
	// calculate line to point distance in 3D
	TVector3 end1(xyzend1[0],xyzend1[1],xyzend1[2]);
	TVector3 end2(xyzend2[0],xyzend2[1],xyzend2[2]);
	TVector3 trackpt(xyz[0],xyz[1],xyz[2]);
	
	float d = ((trackpt-end1).Cross(trackpt-end2)).Mag()/(end2-end1).Mag();
	
	if(d<MinDistanceToPoint)
	  {
	    MinDistanceToPoint=d;
	    MinS=0;
	  }
      }
   
    s = MinS;
    Distance = MinDistanceToPoint;
        
  }



  //----------------------------------------------------------------------
  //  Calculate the closest approach of this track to a given spacepoint
  //   and also the point where this occurs

  void BezierTrack::GetClosestApproach( recob::SpacePoint* sp, double& s,  double& Distance) const
  {
    const double* xyz = sp->XYZ();
    TVector3 Vec(xyz[0],xyz[1],xyz[2]);
    GetClosestApproach(Vec, s, Distance);
  }



  //----------------------------------------------------------------------
  //  Calculate the closest approach of this track to a given position
  //   and also the point where this occurs

  void BezierTrack::GetClosestApproach( TVector3 vec,          double& s,  double& Distance) const
  {
    art::ServiceHandle<util::DetectorProperties> det;
    art::ServiceHandle<geo::Geometry>            geo;

     double iS, xyz[3], MinDistanceToPoint=10000, MinS=0;

    for(int i=0; i!=fBezierResolution; ++i)
      {
	iS=float(i)/fBezierResolution;
	GetTrackPoint(iS, xyz);
	TVector3 trackpt(xyz[0],xyz[1],xyz[2]);

	float d = (vec-trackpt).Mag();
	  
	if(d<MinDistanceToPoint)
	  {
	    MinDistanceToPoint=d;
	    MinS=0;
	  }
      }
   
    s = MinS;
    Distance = MinDistanceToPoint;


  }



  //----------------------------------------------------------------------
  // Calculate the direction of the track by finding the difference
  //  in position between two points. Dir is normalized to 1
  //

  void BezierTrack::GetTrackDirection(double s, double * xyz) const
  {

    if((s<0.5/fBezierResolution)||(s>(1.-0.5/fBezierResolution)))
      {
	throw cet::exception("BezierTrack error: s out of range")<<
	  " cannot query gradient within "<< 0.5/fBezierResolution<<
	  " of track end.  You asked for s = "<<s <<
	  ", which is out of range \n";
      }
    double xyz1[3], xyz2[3];
    GetTrackPoint(s - 0.5/fBezierResolution, xyz1);
    GetTrackPoint(s + 0.5/fBezierResolution, xyz2);
    
    double dx = pow(pow(xyz1[0]-xyz2[0],2)+
		    pow(xyz1[1]-xyz2[1],2)+
		    pow(xyz1[2]-xyz2[2],2),0.5);
    
    for(int i=0; i!=3; ++i)
      {
	xyz[i] = (xyz1[0]-xyz2[0])/dx;
      }
      
  }


  //----------------------------------------------------------------------
  // Friendly versions returning TVector3's 
  // 
  
     
  TVector3 BezierTrack::GetTrackPointV(double s) const
  {
    double xyz[3];
    GetTrackPoint(s,xyz);
    TVector3 ReturnVec(xyz[0],xyz[1],xyz[2]);
    return ReturnVec;
  }

  TVector3 BezierTrack::GetTrackDirectionV(double s) const
  {
    double xyz[3];
    GetTrackDirection(s,xyz);
    TVector3 ReturnVec(xyz[0],xyz[1],xyz[2]);
    return ReturnVec;
  }



  //----------------------------------------------------------------------
  // Method for finding rate at which the track direction changes
  //  (output is local dtheta/dx)

  double BezierTrack::GetCurvature(double s) const
  {
    if((s<1./fBezierResolution)||(s>(1.-1./fBezierResolution)))
      {
	throw cet::exception("BezierTrack error: s out of range")<<
	  " cannot query curvature within "<< 1./fBezierResolution<<
	  " of track end.  You asked for s = "<<s <<
	  ", which is out of range \n";
      }
     
    TVector3 Pos1 = GetTrackPointV(s - 0.5/fBezierResolution);
    TVector3 Pos3 = GetTrackPointV(s + 0.5/fBezierResolution);

    TVector3 Grad1 = GetTrackDirectionV(s - 0.5/fBezierResolution);
    TVector3 Grad2 = GetTrackDirectionV(s );
    TVector3 Grad3 = GetTrackDirectionV(s + 0.5/fBezierResolution);
   
    double dx     = (Pos3-Pos1).Mag();
    double dtheta = (Grad3-Grad2).Angle(Grad2-Grad1);
    
    return (dtheta/dx);
  }


  //----------------------------------------------------------------------
  // Find the RMS curvature along the entire track
  //  (possible measure of multiple scattering)

  double BezierTrack::GetRMSCurvature() const
  {
    double RMS;
    for(int i=1; i!=(fBezierResolution-1); ++i)
      {
	RMS += pow(GetCurvature( float(i)/fBezierResolution),2);
      }
    return (pow(RMS/(fBezierResolution-2),0.5));
  }



  //----------------------------------------------------------------------
  //  return the track length (already calculated, so easy!)
  //  
  
  double BezierTrack::GetLength() const
  {
    return fTrackLength;
  }
  

    

    

}
  




