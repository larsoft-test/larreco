//////////////////////////////////////////////////////////////////////
///
/// \file   BezierCurveHelper.cxx
///
/// \brief  Helper object for interpolating between track segments
///
/// \author B J P Jones
///
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <cmath>
#include <map>
#include <algorithm>
#include "cetlib/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TrackFinder/BezierCurveHelper.h"
#include "Geometry/geo.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "TVector.h"


namespace trkf{
//----------------------------------------------------------------------
// Constructor.
//
BezierCurveHelper::BezierCurveHelper()
{
  fCurveResolution=100;
}

BezierCurveHelper::BezierCurveHelper(int CurveRes)
{
  fCurveResolution=CurveRes;
}

//----------------------------------------------------------------------
// Destructor.
//
BezierCurveHelper::~BezierCurveHelper()
{
}

double BezierCurveHelper::GetSegmentLength(recob::Seed * s1, recob::Seed * s2)
{
  double Length=0;
  std::vector<TVector3> BezierPoints = GetBezierPoints(s1, s2, fCurveResolution);
  for(unsigned int p=0; p!=BezierPoints.size()-1; ++p)
    {
      /*
      std::cout<<"Adding segment " << BezierPoints.at(p+1)[2]<<" "
	       << BezierPoints.at(p+1)[0]<<" "
	       << BezierPoints.at(p+1)[1]<<", "
	       << BezierPoints.at(p)[2]<<" "
	       << BezierPoints.at(p)[0]<<" "
	       << BezierPoints.at(p)[1]<<" "
	       << (BezierPoints.at(p+1)-BezierPoints.at(p)).Mag()<<std::endl;
      */
      Length+=(BezierPoints.at(p+1)-BezierPoints.at(p)).Mag();
    }
  return Length;
}

void BezierCurveHelper::GetDirectionScales(double * Pt1, double * Pt2, double * Dir1, double * Dir2, double * Scales)
{

  double PtSepVec[3];
  for(int i=0; i!=3; ++i) PtSepVec[i]=Pt2[0]-Pt1[0];
  
  double PtSep = pow(
    pow(PtSepVec[0],2) +
    pow(PtSepVec[1],2) +
    pow(PtSepVec[2],2),0.5);

  double Dir1Length = pow(Dir1[0]*Dir1[0]+
			  Dir1[1]*Dir1[1]+
			  Dir1[2]*Dir1[2],0.5);

  double Dir2Length = pow(Dir2[0]*Dir2[0]+
			  Dir2[1]*Dir2[1]+
			  Dir2[2]*Dir2[2],0.5);
  
  if((PtSepVec[0]*Dir1[0]+PtSepVec[1]*Dir1[1]+PtSepVec[2]*Dir1[2])>0)
    Scales[0] =     PtSep / (4.*Dir1Length);
  else
    Scales[0] = 0 - PtSep / (4.*Dir1Length);

  if((PtSepVec[0]*Dir2[0]+PtSepVec[1]*Dir2[1]+PtSepVec[2]*Dir2[2])>0)
    Scales[1] = 0 - PtSep / (4.*Dir2Length);
  else
    Scales[1] =     PtSep / (4.*Dir2Length);

}


void BezierCurveHelper::GetBezierPointXYZ(recob::Seed * s1, recob::Seed * s2, float t, double * xyz)
{
  TVector3 BezierPoint = GetBezierPoint(s1, s2, t);
  xyz[0]=BezierPoint[0];
  xyz[1]=BezierPoint[1];
  xyz[2]=BezierPoint[2];

}

TVector3 BezierCurveHelper::GetBezierPoint(recob::Seed * s1, recob::Seed * s2, float t)
{
  TVector3 ReturnVec3;
  double Pt1[3], Pt2[3], Dir1[3], Dir2[3], Mid1[3], Mid2[3], Dummy[3];
  s1->GetDirection(Dir1, Dummy);
  s2->GetDirection(Dir2, Dummy);
  s1->GetPoint(Pt1,      Dummy);
  s2->GetPoint(Pt2,      Dummy);
  
  double DirScales[2];
  GetDirectionScales(Pt1,Pt2,Dir1,Dir2,DirScales);

  float nt=1.-t;
  for(int i=0; i!=3; i++)
    {
      Mid1[i]=Pt1[i]+Dir1[i]*DirScales[0];
      Mid1[i]=Pt1[i]+Dir1[i]*DirScales[0];
      Mid1[i]=Pt1[i]+Dir1[i]*DirScales[0];
      Mid2[i]=Pt2[i]+Dir2[i]*DirScales[1];
      Mid2[i]=Pt2[i]+Dir2[i]*DirScales[1];
      Mid2[i]=Pt2[i]+Dir2[i]*DirScales[1];
      ReturnVec3[i]=
	nt * nt * nt *        Pt1[i]
	+ 3.* nt * nt * t *   Mid1[i]
	+ 3.* nt * t * t *    Mid2[i]
	+ t * t * t *         Pt2[i];     
    }
  return ReturnVec3;

  //  std::cout<<"BCurve calc: " << Pt1[0]<<" "<<Pt1[1]<<" "<<Pt1[2]<<" ";
  //  std::cout<<"             " << Pt2[0]<<" "<<Pt2[1]<<" "<<Pt2[2]<<" ";
  //  std::cout<<"             " << Dir1[0]<<" "<<Dir1[1]<<" "<<Dir1[2]<<" ";
  // std::cout<<"             " << Dir2[0]<<" "<<Dir2[1]<<" "<<Dir2[2]<<" ";
  //  std::cout<<"             " << t<<" " <<ReturnVec3[0]<<" "<<ReturnVec3[1]<<" "<<ReturnVec3[2]<<" "<<std::endl;

}

std::vector<TVector3> BezierCurveHelper::GetBezierPoints(recob::Seed * s1, recob::Seed * s2, int N)
{
  
  std::vector<TVector3> ReturnVec;
  ReturnVec.resize(N);
  
  double Pt1[3], Pt2[3], Dir1[3], Dir2[3], Mid1[3], Mid2[3], Dummy[3];
  s1->GetDirection(Dir1, Dummy);
  s2->GetDirection(Dir2, Dummy);
  s1->GetPoint(Pt1,      Dummy);
  s2->GetPoint(Pt2,      Dummy);

  double DirScales[2];
  GetDirectionScales(Pt1,Pt2,Dir1,Dir2,DirScales);
  
  //  std::cout<<"BCurve calc: " << Pt1[0]<<" "<<Pt1[1]<<" "<<Pt1[2]<<" ";
  //  std::cout<<"             " << Pt2[0]<<" "<<Pt2[1]<<" "<<Pt2[2]<<" ";
  //  std::cout<<"             " << Dir2[0]<<" "<<Dir2[1]<<" "<<Dir2[2]<<" ";
  //  std::cout<<"             " << Dir2[0]<<" "<<Dir2[1]<<" "<<Dir2[2]<<std::endl;
  for(int i=0; i!=3; i++)
    {
      Mid1[i]=Pt1[i]+Dir1[i]*DirScales[0];
      Mid1[i]=Pt1[i]+Dir1[i]*DirScales[0];
      Mid1[i]=Pt1[i]+Dir1[i]*DirScales[0];
      Mid2[i]=Pt2[i]+Dir2[i]*DirScales[1];
      Mid2[i]=Pt2[i]+Dir2[i]*DirScales[1];
      Mid2[i]=Pt2[i]+Dir2[i]*DirScales[1];
      for(int p=0; p!=N; p++)
	{
	  float t  = float(p) / N;
	  float nt = 1.-t;
	  ReturnVec.at(p)[i] =
	    nt*nt*nt        * Pt1[i]
	    + 3.*nt*nt*t    * Mid1[i]
	    + 3.*nt*t*t     * Mid2[i]
	    + t*t*t         * Pt2[i];
	}
    }
  
  return ReturnVec;
}

  
}
