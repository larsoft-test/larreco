//////////////////////////////////////////////////////////////////////
///
/// \file   SeedFinderService.cxx
///
/// \brief  Service for 3D Seeds from Spacepoints
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
#include "TrackFinder/SeedFinderService.h"
#include "Geometry/geo.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "TH1F.h"
#include "TH2D.h"

namespace trkf{
//----------------------------------------------------------------------
// Constructor.
//
SeedFinderService::SeedFinderService(const fhicl::ParameterSet& pset,
					   art::ActivityRegistry& reg)
{
  reconfigure(pset);
}

//----------------------------------------------------------------------
// Destructor.
//
SeedFinderService::~SeedFinderService()
{
}

//----------------------------------------------------------------------
// Update configuration parameters.
//
void trkf::SeedFinderService::reconfigure(const fhicl::ParameterSet& p)
{
  // Get configuration parameters.

  fSeedLength            = p.get<double>("SeedLength");
  fMinPointsInCluster    = p.get<unsigned int>("MinPointsInCluster");
  fMinPointsInSeed       = p.get<unsigned int>("MinPointsInSeed");
  fAngularDev            = p.get<double>("AngularDev");
}


//
// Set the event ID - only important for organizing
// histogram output.
//
void SeedFinderService::SetEventID(int EventID)
{
  fEventID=EventID;
}



//
// Take a collection of spacepoints sorted in Z and find
// exactly one straight seed, as high in Z as possible.
//
std::vector<recob::Seed*> SeedFinderService::FindSeedExhaustively(std::vector<recob::SpacePoint> Points)
{
  recob::Seed* TrackSeedThisCombo;
  bool KeepChopping=true;
  bool FoundValidSeed=false;
  if(Points.size()<=fMinPointsInCluster) KeepChopping=false;
  while(KeepChopping)
    {
      TrackSeedThisCombo = FindSeedAtEnd(Points);
      if(TrackSeedThisCombo->IsValid())
	{
	  KeepChopping=false;
	  FoundValidSeed=true;
	}
      else
        {
          // if we found an invalid or weak seed, chop off some points and try again
	  if(Points.size() > fMinPointsInSeed)
	    Points.resize(Points.size()-fMinPointsInSeed);
	  else  // or if there is nothing left to chop, leave the loop
	    KeepChopping=false;
        }
    }
  if(!FoundValidSeed)
    {
      TrackSeedThisCombo = new recob::Seed();
    }
  else
    {
      double SeedPoint[3], SeedErr[3];
      TrackSeedThisCombo->GetPoint(SeedPoint,SeedErr);  
      std::cout<<"FindSeedEx : " << SeedPoint[0]<< " " <<
	SeedPoint[1]<< " " <<
	SeedPoint[2]<< " " << std::endl;
    }
  std::vector<recob::Seed*> ReturnVector;
  ReturnVector.push_back(TrackSeedThisCombo);
  return ReturnVector;

}

//
// Take a collection of spacepoints and return a 
// vector of as many straight line seeds as possible
//
std::vector<recob::Seed *> SeedFinderService::FindAsManySeedsAsPossible(std::vector<recob::SpacePoint> Points)
{
  std::vector<recob::Seed*> ReturnVector;
   
  recob::Seed* TrackSeedThisCombo;
  bool KeepChopping=true;
  if(Points.size()<fMinPointsInCluster) KeepChopping=false;
  while(KeepChopping)
    {
      TrackSeedThisCombo = FindSeedAtEnd(Points);
      if(TrackSeedThisCombo->IsValid())
        {
	  ReturnVector.push_back(TrackSeedThisCombo);
        }
      // if enough left, chop off some points and try again
      if(Points.size() > fMinPointsInSeed)
	Points.resize(Points.size()-fMinPointsInSeed);
      else  // or if there is nothing left to chop, leave the loop
	KeepChopping=false;
      
    }
  
  return ReturnVector;
 

}

//
// Start with a vector of hits ordered in Z.
// Put together a seed using the last N hits, find its
// centre, direction and strength. Return this seed for
// further scrutiny.
//
// This method is private since it produces only a possible
// seed, and not a finished product which has passed quality
// checks
//
recob::Seed * SeedFinderService::FindSeedAtEnd(std::vector<recob::SpacePoint> Points)
{
  
  recob::Seed * ReturnSeed;

  std::vector<int> PointsInRange;
  PointsInRange.clear();

  TVector3 HighestZPoint( Points.at(Points.size()-1).XYZ()[0],
			  Points.at(Points.size()-1).XYZ()[1],
			  Points.at(Points.size()-1).XYZ()[2] );

  TVector3 CentreOfPoints(0.,0.,0.);

  for(int index=Points.size()-1; index!=-1; index--)
    {
      // first check z, then check total distance
      //  (much faster, since most will be out of range in z anyway)
 
     
      if( ( HighestZPoint[2] - Points.at(index).XYZ()[2] ) < fSeedLength)
        {
          if( ( ( pow(HighestZPoint[0]-Points.at(index).XYZ()[0],2) +
                  pow(HighestZPoint[1]-Points.at(index).XYZ()[1],2) +
                  pow(HighestZPoint[2]-Points.at(index).XYZ()[2],2) ) < pow(fSeedLength,2)))
	    {
	      PointsInRange.push_back(index);
	      CentreOfPoints[0]+=Points.at(index).XYZ()[0];
	      CentreOfPoints[1]+=Points.at(index).XYZ()[1];
	      CentreOfPoints[2]+=Points.at(index).XYZ()[2];
	    }
        }
    }

  unsigned int NPoints = PointsInRange.size();
  if(NPoints==0) return new recob::Seed();
  CentreOfPoints = (1./float(NPoints)) * CentreOfPoints;

  if(NPoints>fMinPointsInSeed)
    {
      float costheta=0, phi=0, costheta2=0, phi2=0;


      for(unsigned int i=0; i!=PointsInRange.size(); i++)
        {
          TVector3 ThisPoint(Points.at(PointsInRange.at(i)).XYZ()[0],
                             Points.at(PointsInRange.at(i)).XYZ()[1],
                             Points.at(PointsInRange.at(i)).XYZ()[2]);
          TVector3 Step = CentreOfPoints-ThisPoint;
	  
	  if(Step.Z()<0) Step=-Step;
 

          phi       += Step.Phi();
	  phi2      += pow(Step.Phi(),2);

	  // Need to check range of costheta to avoid getting a nan.
	  //  (if CosTheta=1, floating point errors mess us up)
	  if(Step.CosTheta()<0.9999)
	    {
	      costheta     += Step.CosTheta();
	      costheta2    += pow(Step.CosTheta(),2);
	    } 
	  else
	    {
	      costheta+=1;
	      costheta2+=1;
	    }
	  //	  std::cout<<phi<< " " <<phi2<< " " <<costheta << " " <<costheta2<<std::endl;
	}

      
      float sigtheta    = pow((costheta2 / NPoints - pow( costheta/NPoints, 2 )),0.5);
      float sigphi      = pow((phi2   / NPoints - pow( phi/NPoints,   2 )),0.5);
      float meantheta   = acos(costheta      / NPoints);
      float meanphi     = phi        / NPoints;
      
      float AngularDev =  pow(pow(sigtheta,2)+pow(sin(sigtheta)*sigphi,2),0.5); 
      
      std::cout<<"SeedFinderService debug: " << AngularDev<<std::endl;
      
      if(AngularDev<fAngularDev)
	{
     
	  double PtArray[3], DirArray[3];
          
	  PtArray[0]=CentreOfPoints.X();
	  PtArray[1]=CentreOfPoints.Y();
	  PtArray[2]=CentreOfPoints.Z();
	  
	  // calculate direction from average over hits
	  
	  TVector3 SeedDirection;
	  
	  SeedDirection.SetMagThetaPhi(fSeedLength, meantheta, meanphi);
	  
	  DirArray[0]=SeedDirection.X();
	  DirArray[1]=SeedDirection.Y();
	  DirArray[2]=SeedDirection.Z();
	  
	  ReturnSeed = new recob::Seed(PtArray,DirArray);
	}
	else ReturnSeed = new recob::Seed();
    }
  else ReturnSeed = new recob::Seed();
  
  return ReturnSeed;
}


//
// This method will probably become deprecated, which is why it has 
// not been tidied up properly.  Produce output histograms of the
// spacepoints used by the SeedFinderService.
//
void SeedFinderService::ProduceSpacePointPlots(std::vector<std::vector<recob::SpacePoint> > Points)
{
  std::stringstream HistoName("");
  HistoName.flush();

  int xdivs=500;
  int ydivs=500;
  int zdivs=500;
  float xlow=-250;
  float ylow=-250;
  float zlow=0;
  float xhigh=150;
  float yhigh=250;
  float zhigh=1200;

  art::ServiceHandle<art::TFileService> tfs;

  for(unsigned int i=0; i!=Points.size(); i++)
    {
      if(Points.at(i).size()>fMinPointsInCluster)
	{
	  HistoName<<"sp"<<fEventID<<"c"<<i<<"xy";
	  TH2D * hxy = (TH2D*)tfs->make<TH2D>(HistoName.str().c_str(),HistoName.str().c_str(),xdivs,xlow,xhigh,ydivs,ylow,yhigh);
	  HistoName.str("");
	  HistoName.flush();
	  
	  HistoName<<"sp"<<fEventID<<"c"<<i<<"xz";
	  TH2D * hxz = (TH2D*)tfs->make<TH2D>(HistoName.str().c_str(),HistoName.str().c_str(),xdivs,xlow,xhigh,zdivs,zlow,zhigh);
	  HistoName.str("");
	  HistoName.flush();
	  
	  HistoName<<"sp"<<fEventID<<"c"<<i<<"yz";
	  TH2D * hyz = (TH2D*)tfs->make<TH2D>(HistoName.str().c_str(),HistoName.str().c_str(),ydivs,ylow,yhigh,zdivs,zlow,zhigh);
	  HistoName.str("");
	  HistoName.flush();
	  
	  for(unsigned int j=0; j!=Points.at(i).size(); j++)
	    {
	      
	      std::cout <<  Points.at(i).at(j).XYZ()[0]<< ", "
			<<  Points.at(i).at(j).XYZ()[1]<< ", "
			<<  Points.at(i).at(j).XYZ()[2]<< " "
			<<  std::endl;
	      hxy->Fill(Points.at(i).at(j).XYZ()[0],Points.at(i).at(j).XYZ()[1]);
	      hxz->Fill(Points.at(i).at(j).XYZ()[0],Points.at(i).at(j).XYZ()[2]);
		      hyz->Fill(Points.at(i).at(j).XYZ()[1],Points.at(i).at(j).XYZ()[2]);
	    }
	}
    }
}


// Find the nearby hits to this seed

std::vector<std::vector<recob::Hit*> > SeedFinderService::CollectSeedHits(std::vector<recob::Seed *> TheSeed, std::vector<recob::Hit*> TheHits)
{
  
  std::vector<std::vector<recob::Hit*> > ReturnVec;
  for(std::vector<recob::Hit*>::const_iterator it=TheHits.begin();
      it!=TheHits.end();
      it++)
    {
      

    }

  return ReturnVec;
}


// Refit seeds using nearby hits
std::vector<recob::Seed*> SeedFinderService::RefitSeeds(std::vector<recob::Seed *> TheSeeds, std::vector<recob::Hit*> TheHits)
{
  std::vector<recob::Seed*> ReturnVec;
  std::vector<std::vector<recob::Hit*> > NearbyHits = CollectSeedHits(TheSeeds,TheHits);

  return ReturnVec;
}  

  
  
}
