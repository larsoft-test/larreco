#include "TrackFinder/BezierTrack.h"

#include "Calorimetry/CalorimetryAlg.h"
#include "Geometry/geo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/SpacePoint.h"
#include "SimpleTypesAndConstants/PhysicalConstants.h"
#include "TrackFinder/BezierCurveHelper.h"
#include "Utilities/DetectorProperties.h"

#include "TVector3.h"

#include "cetlib/exception.h"

namespace trkf {

  //----------------------------------------------------------------------
  // Constructor from Base object (for analysis)
  //
  BezierTrack::BezierTrack(recob::Track btb)
    : recob::Track(btb)
  {
    fID = btb.ID();
    fBezierResolution=1000;
    CalculateSegments();
  }


  //----------------------------------------------------------------------
  // Constructor from seed vector (for production)
  //
  BezierTrack::BezierTrack(std::vector<recob::Seed*> const SeedCol )
    : recob::Track()
  {
    if(SeedCol.size()>0)
      {
  double Pt[3], Dir[3], PtErr[3], DirErr[3];
  double FirstPt[3];
  double FirstDir[3];
  SeedCol.at(0)->GetPoint(     Pt,  PtErr  );
  SeedCol.at(0)->GetDirection( Dir, DirErr );
  for(int i=0; i!=3; ++i)
    {
      double BigN=10000;
      FirstPt[i]=Pt[i]-((1.+1/BigN) * Dir[i]);
      FirstDir[i]=Dir[i]/BigN;
    }
  fSeedCollection.push_back(new recob::Seed(FirstPt, FirstDir, PtErr, DirErr));

      }
    for(size_t i=0; i!=SeedCol.size(); ++i)
      {
  fSeedCollection.push_back(new recob::Seed(*SeedCol.at(i)));
      }
    if(SeedCol.size()>0)
      {
  double Pt[3], Dir[3], PtErr[3], DirErr[3];
  double LastPt[3], LastDir[3];
  SeedCol.at(SeedCol.size()-1)->GetPoint(     Pt,  PtErr  );
  SeedCol.at(SeedCol.size()-1)->GetDirection( Dir, DirErr );
  for(int i=0; i!=3; ++i)
    {
      double BigN=10000;
      LastPt[i]=Pt[i]+((1.+1/BigN) * Dir[i]);
      LastDir[i]=Dir[i]/BigN;
    }
  fSeedCollection.push_back(new recob::Seed(LastPt, LastDir, PtErr, DirErr));
      }
    std::cout<<"Constructing BTrack from " << fSeedCollection.size()<<" seeds."<<std::endl;
    CalculateSegments();
    fBezierResolution=1000;
  }


  //----------------------------------------------------------------------
  // Constructor from track coordinates
  //
  BezierTrack::BezierTrack(std::vector<TVector3> Pos,
         std::vector<TVector3> Dir,
         std::vector<std::vector<double> > dQdx,
         int const& id)
    : recob::Track(Pos, Dir, dQdx, std::vector<double>(2, util::kBogusD), id)
  {
    fBezierResolution=1000;
    CalculateSegments();
  }


  //----------------------------------------------------------------------
  // Default constructor
  //
  BezierTrack::~BezierTrack()
  {
  }



  //----------------------------------------------------------------------
  // Get the track pitch at some s for a particular view
  //

  double BezierTrack::GetTrackPitch(geo::View_t view, double s, double WirePitch, unsigned int c, unsigned int t)
  {
    static std::map<geo::View_t, bool>          DoneCalc;
    static std::map<geo::View_t, TVector3>      PitchVecs;

    if(!(DoneCalc[view]))
      {
  unsigned int p;
  unsigned int pTop;

  art::ServiceHandle<geo::Geometry> geom;

  pTop=geom->Cryostat(c).TPC(t).Nplanes();
  for(p=0; p!=pTop; ++p)
    if(geom->Cryostat(c).TPC(t).Plane(p).View() == view) break;

  double WireEnd1[3], WireEnd2[3];

  geom->WireEndPoints(c, t, p, 0, WireEnd1, WireEnd2);

  double WirePitch = geom->WirePitch(view);

  TVector3 WireVec;
  for(size_t i=0; i!=3; ++i)
    WireVec[i]= WireEnd2[i] -WireEnd1[i];
  PitchVecs[view] = (1.0 / WirePitch) * WireVec.Cross(TVector3(1,0,0)).Unit();

  DoneCalc[view]=true;
      }

    TVector3 TrackDir = GetTrackDirectionV(s).Unit();

    return 1. / (PitchVecs[view].Dot(TrackDir));
  }



  //------------------------------------------------------
  // Build an anab::Calorimetry object for this track

  anab::Calorimetry BezierTrack::GetCalorimetryObject(std::vector<art::Ptr<recob::Hit> > Hits, geo::SigType_t sigtype, calo::CalorimetryAlg const& calalg)
  {

    art::ServiceHandle<geo::Geometry> geom;

    std::vector<double>  dEdx;
    std::vector<double>  dQdx;
    std::vector<double>  resRange;
    std::vector<double>  deadwire;
    double Range = GetLength();
    std::vector<double>  TrkPitch;



    double s, d;
    TVector3 PlanePitchDirection;


    for(size_t i=0; i!=Hits.size(); ++i)
      {
  if(Hits.at(i)->SignalType()==sigtype)
    {
      geo::View_t view = Hits.at(i)->View();
      double WirePitch = geom->WirePitch(view);
      GetClosestApproach(Hits.at(i),s,d);
      double ThisPitch = GetTrackPitch( view, s, WirePitch );
      TrkPitch.push_back(ThisPitch);
      resRange.push_back(Range - s*Range);
      dQdx.push_back(Hits.at(i)->Charge()/ThisPitch);
      dEdx.push_back( calalg.dEdx_AMP(Hits.at(i), ThisPitch, false) );
      std::cout<<Hits.at(i)->Charge()<< " " << ThisPitch<< " " << dEdx.at(dEdx.size()-1)<<std::endl;
    }


      }


    double KineticEnergy;

    for(size_t i=0; i!=dEdx.size(); ++i)
      KineticEnergy+=dEdx.at(i);

    return anab::Calorimetry(KineticEnergy, dEdx, dQdx, resRange, deadwire, Range, TrkPitch);


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


  Pt[0]=(fXYZ.at(i))[0];
  Pt[1]=(fXYZ.at(i))[1];
  Pt[2]=(fXYZ.at(i))[2];

  Dir[0]=(fDir.at(i))[0];
  Dir[1]=(fDir.at(i))[1];
  Dir[2]=(fDir.at(i))[2];

  fSeedCollection.push_back(new recob::Seed(Pt,Dir));
      }

  }




  //----------------------------------------------------------------------
  void BezierTrack::FillTrajectoryVectors()
  {
    art::ServiceHandle<geo::Geometry> geom;
    int NPlanes = geom->Nplanes();
    fdQdx.resize(NPlanes) ;
    double Pt[3], Dir[3], ErrPt[3], ErrDir[3];
    for(std::vector<recob::Seed*>::const_iterator it=fSeedCollection.begin();
  it!=fSeedCollection.end(); ++it)
      {

  (*it)->GetPoint(Pt,ErrPt);
  (*it)->GetDirection(Dir, ErrDir);
  TVector3 Point(Pt[0],Pt[1],Pt[2]);
  TVector3 Direction(Dir[0],Dir[1],Dir[2]);

  for(int i=0; i!=NPlanes; ++i)
    fdQdx.at(i).push_back(0);
  fXYZ.push_back(Point);
  fDir.push_back(Direction);
      }
  }


  //----------------------------------------------------------------------
  // Calculate the lengths of each bezier segment and fill useful
  // calculational data members
  //

  void BezierTrack::CalculateSegments()
  {
    std::cout<<"Entering calculatesegments method"<<std::endl;
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
    if(NSegments()==0)
      {
  std::cout<<"Filling traj vectors             "<<std::endl;
  if(fSeedCollection.size()!=0)  FillTrajectoryVectors();
  std::cout<<"Filed traj vectors               "<<std::endl;
      }

    std::cout<<"Loading bezier track with ID "<<  ID()<<std::endl;
    BezierCurveHelper bhlp(100);

    int Segments = NSegments();

    fTrackLength=0;
    bool FirstSeg=true;
    for(int i=0; i!=Segments; ++i)
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
    std::cout<<"Total track length : " <<fTrackLength<< ", N = " << NSegments()<<std::endl;
  }


  //----------------------------------------------------------------------
  // Find the point which is fraction s along the track and return
  // as a double[3]
  //
  void BezierTrack::GetTrackPoint(double s, double * xyz) const
  {
    if((s>=1)||(s<=0))
      {
  // catch these easy floating point errors
  if((s>0.9999)&&(s<1.00001))
    {
      TVector3 End1 = fXYZ.at(fXYZ.size()-1);
      xyz[0]=End1[0];
      xyz[1]=End1[1];
      xyz[2]=End1[2];

      return;
    }
  else if((s<0.0001)&&(s>-0.00001))
    {
      TVector3 End0 = fXYZ.at(0);
      xyz[0]=End0[0];
      xyz[1]=End0[1];
      xyz[2]=End0[2];

      return;
    }

  // otherwise complain about the mistake
  else
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
  //   - this fast version is for ennumerating which members of a hit
  //     collection are within d, how far each is and where the closest
  //     approach occurs.  This version is optimized for speed for this
  //     application
  void BezierTrack::GetClosestApproaches( art::PtrVector<recob::Hit> hits,     std::vector<double>& s,  std::vector<double>& Distances) const
  {
    std::vector< art::Ptr<recob::Hit> > hitv;
    for(size_t h = 0; h < hits.size(); ++h) hitv.push_back(hits[h]);
    this->GetClosestApproaches(hitv, s, Distances);
  }

  void BezierTrack::GetClosestApproaches( std::vector< art::Ptr<recob::Hit> > hits,     std::vector<double>& s,  std::vector<double>& Distances) const
  {
    art::ServiceHandle<util::DetectorProperties> det;
    art::ServiceHandle<geo::Geometry>            geo;

    s.clear();
    Distances.clear();


    // Pull all the relevant information out of hits and into simple vectors

    std::vector<TVector3> HitEnd1s, HitEnd2s;
    std::vector<double> WireLengths;

    HitEnd1s.resize(hits.size());
    HitEnd2s.resize(hits.size());
    WireLengths.resize(hits.size());

    unsigned int c1, t1, p1, w1;
    double End1[3], End2[3];

    size_t NHits = hits.size();

    for(size_t i=0; i!=NHits; i++)
      {
  Distances.push_back(10000);
  s.push_back(-1);

  geo->ChannelToWire(hits.at(i)->Channel(),c1,t1,p1,w1);

  geo->WireEndPoints(c1,t1,p1,w1,End1,End2);

  HitEnd1s.at(i)[0]= HitEnd2s.at(i)[0]= det->ConvertTicksToX(hits.at(i)->PeakTime(),p1,t1,c1);
  HitEnd1s.at(i)[1]= End1[1];
  HitEnd2s.at(i)[1]= End2[1];
  HitEnd1s.at(i)[2]= End1[2];
  HitEnd2s.at(i)[2]= End2[2];

  WireLengths.at(i)=((HitEnd1s.at(i)-HitEnd2s.at(i)).Mag());
      }


    double iS;

    for(int ipt=0; ipt!=fBezierResolution; ++ipt)
      {
  iS=float(ipt)/fBezierResolution;
  TVector3 trackpt = GetTrackPointV(iS);

  for(size_t ihit=0; ihit!=NHits; ++ihit)
    {
      float d = ((trackpt-HitEnd1s.at(ihit)).Cross(trackpt-HitEnd2s.at(ihit))).Mag()/WireLengths.at(ihit);

      if(d<Distances.at(ihit))
        {
    Distances.at(ihit)=d;
    s.at(ihit)=iS;
        }
    }
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
      MinS=iS;
    }
      }

    s = MinS;
    Distance = MinDistanceToPoint;

  }




  //----------------------------------------------------------------------
  //  Calculate the closest approach of this track to a given hit
  //   and also the point where this occurs
  //

  void BezierTrack::GetClosestApproach( art::Ptr<recob::Hit> hit,       double& s,  double& Distance) const
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
      MinS=iS;
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
      MinS=iS;
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
  if     (( s < 0.5 / fBezierResolution ) && ( s > -0.00001 ) )
    s = 0.5 / fBezierResolution;
  else if((s>(1.-0.5/fBezierResolution)) && ( s < 1.00001)) s = 1.-0.5/fBezierResolution;
  else
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
  xyz[i] = (xyz1[i]-xyz2[i])/dx;
      }

  }


  //----------------------------------------------------------------------
  // Friendly versions returning TVector3's
  //


  //----------------------------------------------------------------------
  TVector3 BezierTrack::GetTrackPointV(double s) const
  {
    double xyz[3];
    GetTrackPoint(s,xyz);
    TVector3 ReturnVec(xyz[0],xyz[1],xyz[2]);
    return ReturnVec;
  }

  //----------------------------------------------------------------------
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
    TVector3 Pos2 = GetTrackPointV(s );
    TVector3 Pos3 = GetTrackPointV(s + 0.5/fBezierResolution);

    double dx     = (1./fBezierResolution)*fTrackLength;
    double dtheta = (Pos3-Pos2).Angle(Pos2-Pos1);

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



  //----------------------------------------------------------------------
  //  Fill the dQdx vector for the track, based on a set of hits
  //    provided

  void BezierTrack::CalculatedQdx(art::PtrVector<recob::Hit> Hits)
  {
    fdQdx.clear();

    std::map<int, std::map<int, double> > hitmap;
    //        ^              ^       ^
    //       view            seg    charge

    for(size_t i=0; i!=Hits.size(); ++i)
      {
  double Distance, S;
  GetClosestApproach(Hits.at(i), S, Distance);
  //  std::cout<<"hit " << i <<"  " <<  S << " " << WhichSegment(S) << " " << Hits.at(i)->View()<<std::endl;

  (hitmap[Hits.at(i)->View()])[WhichSegment(S)] += Hits.at(i)->Charge();
      }

    int NSeg = NSegments();

    for(std::map<int,std::map<int,double> >::const_iterator itview=hitmap.begin();
  itview!=hitmap.end(); ++itview)
      {
  std::vector<double> ThisViewdQdx;
  ThisViewdQdx.resize(NSeg);
  for(std::map<int,double>::const_iterator itseg=itview->second.begin();
      itseg!=itview->second.end(); ++itseg)
    {
      //   std::cout << itview->first<<" " <<itseg->first << " " <<itseg->second<<std::endl;
      int seg = itseg->first;

      // need to nudge hits which fell outside the track

      if((seg>-1) && (seg<NSeg))
        ThisViewdQdx[seg] = (itseg->second / fSegmentLength[seg]);
    }
  fdQdx.push_back(ThisViewdQdx);
      }
    std::cout<<"Size of dQdx structure :" <<fdQdx.size()<<" : ";
    for(size_t i=0; i!=fdQdx.size(); i++)
      std::cout<<fdQdx.at(i).size()<<", ";
    std::cout<<std::endl;

  }



  //----------------------------------------------------------------------
  //  Fill the dQdx vector for the track, based on a set of hits
  //    provided - optimized version

  void BezierTrack::CalculatedQdx(art::PtrVector<recob::Hit> Hits, std::vector<double> SValues)
  {
    fdQdx.clear();

    std::map<int, std::map<int, double> > hitmap;
    //        ^              ^       ^
    //       view            seg    charge

    for(size_t i=0; i!=Hits.size(); ++i)
      {

  (hitmap[Hits.at(i)->View()])[WhichSegment(SValues.at(i))] += Hits.at(i)->Charge();
      }

    int NSeg = NSegments();

    for(std::map<int,std::map<int,double> >::const_iterator itview=hitmap.begin();
  itview!=hitmap.end(); ++itview)
      {
  std::vector<double> ThisViewdQdx;
  ThisViewdQdx.resize(NSeg);
  for(std::map<int,double>::const_iterator itseg=itview->second.begin();
      itseg!=itview->second.end(); ++itseg)
    {
      //   std::cout << itview->first<<" " <<itseg->first << " " <<itseg->second<<std::endl;
      int seg = itseg->first;

      // need to nudge hits which fell outside the track

      if((seg>-1) && (seg<NSeg))
        ThisViewdQdx[seg] = (itseg->second / fSegmentLength[seg]);
    }
  fdQdx.push_back(ThisViewdQdx);
      }
    std::cout<<"Size of dQdx structure :" <<fdQdx.size()<<" : ";
    for(size_t i=0; i!=fdQdx.size(); i++)
      std::cout<<fdQdx.at(i).size()<<", ";
    std::cout<<std::endl;

  }


  //----------------------------------------------------------------------
  //  Get dQdx for a particular S value
  //

  double BezierTrack::GetdQdx(double s, unsigned int view) const
  {
    view++;
    if((s<0)||(s>1))
      {
  throw cet::exception("Bezier dQdx: S out of range")
    <<"Bezier track S value of " << s <<" is not in the range 0<S<1"
    <<std::endl;
      }
    if((view<0)||view>(fdQdx.size()-1))
      {
  throw cet::exception("Bezier dQdx: view out of range")
    <<"Bezier track view value of " << view <<" is not in the range "
    <<"of stored views, 0 < view < " << fdQdx.size()
    <<std::endl;
      }
    int Segment = WhichSegment(s);
    //    if((Segment>=fdQdx[view].size())||(Segment<0))
    //      {
    //  throw cet::exception("Bezier dQdx: segment of range")
    //    <<"Bezier track segment value of " << Segment <<" is not in the range "
    //    <<"of stored segments, 0 < seg < " << fdQdx[view].size()
    //    <<std::endl;
    //     }
    if( ( Segment < (int)fdQdx[view].size() ) && (Segment>1))
      return fdQdx[view][Segment];
    else
      {
  std::cout<<"GetdQdx : Bad s  " <<s<<", " <<Segment<<std::endl;
  return 0;
      }
  }


  //----------------------------------------------------------------------
  //  Get RMS dQdx for a particular view
  //

  double BezierTrack::GetViewdQdx(unsigned int view) const
  {
    view++;
    if((view<0)||view>fdQdx.size()-1)
      {
  throw cet::exception("view out of range")
    <<"Bezier track view value of " << view <<" is not in the range "
    <<"of stored views, 0 < view < " << fdQdx.size()
    <<std::endl;
      }

    size_t NSeg = NSegments();
    double TotaldQdx;
    for(size_t i=0; i!=NSeg; ++i)
      {
  TotaldQdx+=fdQdx.at(view).at(i)*fSegmentLength[i];
      }
    return TotaldQdx / fTrackLength;
  }



  //----------------------------------------------------------------------
  //  Get total charge for a particular view
  //

  double BezierTrack::GetTotalCharge(unsigned int View) const
  {
    View++;
    return GetViewdQdx(View) * fTrackLength;
  }



  //----------------------------------------------------------------------
  //  Find which track segment a particular S lies in
  //
  int BezierTrack::WhichSegment(double s) const
  {
    int ReturnVal=-1;
    for(size_t i=0; i!=fCumulativeLength.size()-1; i++)
      {
  if( (fCumulativeLength.at(i)/fTrackLength <= s)
      &&(fCumulativeLength.at(i+1)/fTrackLength >= s))
    ReturnVal=i;
      }
    if( ( s<(fCumulativeLength.at(0)/fTrackLength)) && (s>-0.0001)) ReturnVal=0;
    return ReturnVal;
  }


  //----------------------------------------------------------------------
  // Return the number of track segments
  //

  int BezierTrack::NSegments() const
  {
    return NumberTrajectoryPoints();
  }


  //----------------------------------------------------------------------
  // Return a fresh copy of the RecoBase track object
  //

  recob::Track BezierTrack::GetBaseTrack()
  {
    std::vector<double> mom(2, util::kBogusD);
    recob::Track TheTrack(fXYZ, fDir, fdQdx, mom, ID());
    return TheTrack;
  }


  std::vector<recob::SpacePoint> BezierTrack::FillMySpacePoints(int N)
  {
    std::vector<recob::SpacePoint> spts(N);
    for(int i=0; i!=N; ++i)
      {
  double xyz[3];
  double ErrXYZ[3]={0.1,0.1,0.1};
  GetTrackPoint(float(i)/N,xyz);
  recob::SpacePoint TheSP(xyz, ErrXYZ, util::kBogusD, i);

  spts[i] = TheSP;

      }

    return spts;
  }


  recob::Track BezierTrack::GetJoinedBaseTrack(BezierTrack * Other)
  {

    std::vector<TVector3> XYZ = fXYZ;
    std::vector<TVector3> Dir = fDir;
    std::vector<std::vector<double> > dQdx = fdQdx;

    if((dQdx.size()!=Other->fdQdx.size()))
      {
  throw cet::exception("BezierTrack combination error ")<<
    "Two vectors do not share the same dQdx vector size" <<
    dQdx.size()<<", "<< Other->fdQdx.size()<<std::endl;
      }

    XYZ.insert(XYZ.end(), Other->fXYZ.begin(), Other->fXYZ.end());
    Dir.insert(Dir.end(), Other->fDir.begin(), Other->fDir.end());
    for(size_t i=0; i!=dQdx.size(); i++)
      dQdx.at(i).insert(dQdx.at(i).end(), Other->fdQdx.at(i).begin(), Other->fdQdx.at(i).end());

    return recob::Track(XYZ,Dir,dQdx);

  }


  recob::Track BezierTrack::GetJoinedPartBaseTrack(BezierTrack * Other, int mybegin, int myend, int theirbegin, int theirend)
  {

    std::vector<TVector3> XYZ ;
    std::vector<TVector3> Dir ;
    std::vector<std::vector<double> > dQdx ;

    dQdx.resize(fdQdx.size());

    if((dQdx.size()!=Other->fdQdx.size()))
      {
  throw cet::exception("BezierTrack combination error ")<<
    "Two vectors do not share the same dQdx vector size" <<
    dQdx.size()<<", "<< Other->fdQdx.size()<<std::endl;
      }


    XYZ.insert(XYZ.end(), fXYZ.begin()+mybegin, fXYZ.begin()+myend);
    Dir.insert(Dir.end(), fDir.begin()+mybegin, fDir.begin()+myend);
    for(size_t i=0; i!=dQdx.size(); i++)
      dQdx.at(i).insert(dQdx.at(i).end(), fdQdx.at(i).begin()+mybegin, fdQdx.at(i).begin()+myend);

    XYZ.insert(XYZ.end(), Other->fXYZ.begin()+theirbegin, Other->fXYZ.begin()+theirend);
    Dir.insert(Dir.end(), Other->fDir.begin()+theirbegin, Other->fDir.begin()+theirend);
    for(size_t i=0; i!=dQdx.size(); i++)
      dQdx.at(i).insert(dQdx.at(i).end(), Other->fdQdx.at(i).begin()+theirbegin, Other->fdQdx.at(i).end()+theirend);

    return recob::Track(XYZ,Dir,dQdx);

  }


  recob::Track BezierTrack::GetReverseBaseTrack()
  {
    std::vector<TVector3> XYZ = fXYZ;
    std::vector<TVector3> Dir = fDir;
    std::vector<std::vector<double> > dQdx = fdQdx;

    for(size_t i=0; i!=dQdx.size(); i++)
      reverse(dQdx.at(i).begin(), dQdx.at(i).end());

    reverse(XYZ.begin(), XYZ.end());
    reverse(Dir.begin(), Dir.end());

    // flip trajectory points
    for(size_t i=0; i!=Dir.size(); i++) Dir[i]=-Dir[i];

    std::vector<double> mom(2, util::kBogusD);
    recob::Track TheTrack(XYZ, Dir, dQdx, mom, ID());
    return TheTrack;
  }
}





