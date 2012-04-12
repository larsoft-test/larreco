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
  double Length;
  std::vector<TVector3> BezierPoints = GetBezierPoints(s1, s2, fCurveResolution);
  for(unsigned int p=0; p!=BezierPoints.size()-1; p++)
    {
      Length+=(BezierPoints.at(p+1)-BezierPoints.at(p)).Mag();
    }
  return Length;
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
  s1->GetPoint(Pt2,      Dummy);
  

  for(int i=0; i!=3; i++)
    {
      Mid1[i]=Pt1[i]+Dir1[i]/2.0;
      Mid1[i]=Pt1[i]+Dir1[i]/2.0;
      Mid1[i]=Pt1[i]+Dir1[i]/2.0;
      Mid2[i]=Pt2[i]-Dir2[i]/2.0;
      Mid2[i]=Pt2[i]-Dir2[i]/2.0;
      Mid2[i]=Pt2[i]-Dir2[i]/2.0;
      ReturnVec3[i]=
	pow(1.-t,3) * Pt1[i]
	+ 3.*pow(1.-t, 2)*t   * Mid1[i]
	+ 3.*(1.-t)*pow(t, 2) * Mid2[i]
	+ pow(t,3) * Pt2[2];     
    }
  return ReturnVec3;
}

std::vector<TVector3> BezierCurveHelper::GetBezierPoints(recob::Seed * s1, recob::Seed * s2, int N)
{
  
  std::vector<TVector3> ReturnVec;
  ReturnVec.resize(N);
  
  double Pt1[3], Pt2[3], Dir1[3], Dir2[3], Mid1[3], Mid2[3], Dummy[3];
  s1->GetDirection(Dir1, Dummy);
  s2->GetDirection(Dir2, Dummy);
  s1->GetPoint(Pt1,      Dummy);
  s1->GetPoint(Pt2,      Dummy);
  

  for(int i=0; i!=3; i++)
    {
      Mid1[i]=Pt1[i]+Dir1[i]/2.0;
      Mid1[i]=Pt1[i]+Dir1[i]/2.0;
      Mid1[i]=Pt1[i]+Dir1[i]/2.0;
      Mid2[i]=Pt2[i]-Dir2[i]/2.0;
      Mid2[i]=Pt2[i]-Dir2[i]/2.0;
      Mid2[i]=Pt2[i]-Dir2[i]/2.0;
      for(int p=0; p!=N; p++)
	{
	  float t = float(p) / N;
	  ReturnVec.at(p)[i]=
	    pow(1.-t,3) * Pt1[i]
	    + 3.*pow(1.-t, 2)*t   * Mid1[i]
	    + 3.*(1.-t)*pow(t, 2) * Mid2[i]
	    + pow(t,3) * Pt2[2];
	}
    }
  
  return ReturnVec;
}

  
}
