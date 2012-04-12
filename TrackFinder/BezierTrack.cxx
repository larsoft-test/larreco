#include "TrackFinder/BezierTrack.h"
#include "TrackFinder/BezierCurveHelper.h"
#include "Utilities/DetectorProperties.h"
#include "Geometry/geo.h"



namespace trkf {


  BezierTrack::BezierTrack(recob::BezierTrackBase btb):
    recob::BezierTrackBase(btb)
  {
    CalculateSegments();
  }

  BezierTrack::BezierTrack(std::vector<recob::Seed*> SeedCol):
    recob::BezierTrackBase()
  {
    fSeedCollection = SeedCol;
    CalculateSegments();
  }

  

  void BezierTrack::FillSeedVector()
  {
    int NSeg=NSegments();
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
	    FirstSeg=true;
	    float SegmentLength=bhlp.GetSegmentLength(fSeedCollection.at(i-1),fSeedCollection.at(i));
	    
	    fTrackLength+=SegmentLength;
	    fSegmentLength.push_back(SegmentLength);
	  }
	fCumulativeLength.push_back(fTrackLength);
      }      
  }


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
            if(  (   (fCumulativeLength.at(i-1) / fTrackLength) <  s)
                 &&( (fCumulativeLength.at(i)   / fTrackLength) >= s))
              {

                double locals = (s * fTrackLength - fCumulativeLength[i-1])/fSegmentLength[i];
                bhlp.GetBezierPointXYZ(fSeedCollection.at(i-1),fSeedCollection.at(i),locals, xyz);
              }

          }
      }
  }
  

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


   


    

    

}
  




