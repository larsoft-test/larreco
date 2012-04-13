#include "TrackFinder/BezierTrack.h"
#include "TrackFinder/BezierCurveHelper.h"
#include "Utilities/DetectorProperties.h"
#include "Geometry/geo.h"



namespace trkf {

  //----------------------------------------------------------------------
  // Constructor from Base object (for analysis)
  //
  BezierTrack::BezierTrack(recob::BezierTrackBase btb):
    recob::BezierTrackBase(btb)
  {
    CalculateSegments();
    fBezierResolution=1000;
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
	//	std::cout<<"Adding in seed with coord " << fPtX.at(i)<<" " <<fPtY.at(i)<<" " <<fPtZ.at(i)<<std::endl; 
	//	std::cout<<"                          " << fDirX.at(i)<<" " <<fDirY.at(i)<<" " <<fDirZ.at(i)<<std::endl; 

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
	    std::cout<<"No points in bezier track to make segments!"<<std::endl;
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
	    //  std::cout<<"Pushing back segment with length" << SegmentLength<<std::endl;
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
	std::cout<<"Error: Bezier track expects coordinate between 0 and 1"<<std::endl;
        for(int i=0; i!=3; xyz[i]=0.);
      }
    else
      {
	BezierCurveHelper bhlp;
        for(unsigned int i=1; i!=fCumulativeLength.size(); i++)
          {
	    //  std::cout<<"Seeking track point: " << fCumulativeLength.at(i-1)/fTrackLength << " " << fCumulativeLength.at(i)/fTrackLength<<std::endl;
            if(  (   (fCumulativeLength.at(i-1) / fTrackLength) <=  s)
                 &&( (fCumulativeLength.at(i)   / fTrackLength) > s))
              {
		
                double locals = (s * fTrackLength - fCumulativeLength[i-1])/fSegmentLength[i-1];
		std::cout<<"s, fCum, seglength, locals " <<s<<" " << fCumulativeLength[i-1]<< " " <<fSegmentLength[i]<<" " << locals<<std::endl;
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

  void BezierTrack::GetClosestApproach( recob::Hit* hit,       double& s,  double& Distance)
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

  void BezierTrack::GetClosestApproach( recob::SpacePoint* sp, double& s,  double& Distance)
  {
    const double* xyz = sp->XYZ();
    TVector3 Vec(xyz[0],xyz[1],xyz[2]);
    GetClosestApproach(Vec, s, Distance);
  }



  //----------------------------------------------------------------------
  //  Calculate the closest approach of this track to a given position
  //   and also the point where this occurs

  void BezierTrack::GetClosestApproach( TVector3 vec,          double& s,  double& Distance)
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

   


    

    

}
  




