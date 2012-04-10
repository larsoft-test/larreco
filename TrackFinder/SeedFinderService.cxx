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
  fMinDirectionStrength  = p.get<float>("MinDirectionStrength");

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
	  if(Points.size() >= fMinPointsInSeed)
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
          KeepChopping=false;
	  ReturnVector.push_back(TrackSeedThisCombo);
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

  float HighestZ = Points.at(Points.size()-1).XYZ()[2];
  float HighestZX = Points.at(Points.size()-1).XYZ()[0];
  float HighestZY = Points.at(Points.size()-1).XYZ()[1];

  double CentreOfPointsX=0,CentreOfPointsY=0,CentreOfPointsZ=0;

  for(int index=Points.size()-1; index!=-1; index--)
    {
      // first check z, then check total distance
      //  (much faster, since most will be out of range in z anyway)
 
     
      if( ( HighestZ - Points.at(index).XYZ()[2] ) < fSeedLength)
        {
          if( ( ( pow(HighestZ-Points.at(index).XYZ()[2],2) +
                  pow(HighestZX-Points.at(index).XYZ()[0],2) +
                  pow(HighestZY-Points.at(index).XYZ()[1],2) ) < pow(fSeedLength,2)))
	    {
	      PointsInRange.push_back(index);
	      CentreOfPointsX+=Points.at(index).XYZ()[0];
	      CentreOfPointsY+=Points.at(index).XYZ()[1];
	      CentreOfPointsZ+=Points.at(index).XYZ()[2];
	    }
        }
    }

  int NPoints = PointsInRange.size();
  if(NPoints==0) return new recob::Seed();
  CentreOfPointsX /= NPoints;
  CentreOfPointsY /= NPoints;
  CentreOfPointsZ /= NPoints;

  art::ServiceHandle<art::TFileService> tfs;
  std::stringstream HistoName("");

   
  
  std::vector<std::map<int,std::map<int, int> > > HitMapVec;
  HitMapVec.resize(4);
  for(int i=0; i!=4; i++) HitMapVec.at(i).clear();
 
  int CosThetaBins=20;
  int PhiBins=20;

  std::cout<<"Points in range " << PointsInRange.size()<<std::endl;
  TVector3 HighestZPoint(HighestZX,HighestZY,HighestZ);
  TVector3 CentreOfPoints(CentreOfPointsX,CentreOfPointsY,CentreOfPointsZ);

  if(PointsInRange.size()>fMinPointsInSeed)
    {

      for(unsigned int i=0; i!=PointsInRange.size(); i++)
        {
          TVector3 ThisPoint(Points.at(PointsInRange.at(i)).XYZ()[0],
                             Points.at(PointsInRange.at(i)).XYZ()[1],
                             Points.at(PointsInRange.at(i)).XYZ()[2]);
          TVector3 Step = CentreOfPoints-ThisPoint;
	  
	  if(Step.Z()<0) Step=-Step;
 
          float costheta = Step.CosTheta();
          float phi      = Step.Phi();
	  
	  int costhetabin, phibin;	  

          costhetabin = (int) ((costheta+1.)/2. * CosThetaBins - 0.5);
	  phibin = (int) ((phi/6.284)*PhiBins-0.5);
 	  HitMapVec.at(0)[costhetabin][phibin] = HitMapVec.at(0)[costhetabin][phibin] + 1;

          costhetabin = (int) ((costheta+1.)/2. * CosThetaBins );
	  phibin = (int) ((phi/6.284)*PhiBins-0.5);
 	  HitMapVec.at(1)[costhetabin][phibin] = HitMapVec.at(1)[costhetabin][phibin] + 1;

          costhetabin = (int) ((costheta+1.)/2. * CosThetaBins - 0.5);
	  phibin = (int) ((phi/6.284)*PhiBins );
 	  HitMapVec.at(2)[costhetabin][phibin] = HitMapVec.at(2)[costhetabin][phibin] + 1;

          costhetabin = (int) ((costheta+1.)/2. * CosThetaBins);
	  phibin = (int) ((phi/6.284)*PhiBins);
 	  HitMapVec.at(3)[costhetabin][phibin] = HitMapVec.at(3)[costhetabin][phibin] + 1;
	  
        }
      
      double MaxValues[4]       = {0.};
      double MaxCosThetaBins[4] = {0.};
      double MaxPhiBins[4]   = {0.};
      for(int mapid=0; mapid!=4; mapid++)
	{
	  
	  // find most populated bin in hitmap
	  for(std::map<int,std::map<int, int> >::const_iterator itcostheta=HitMapVec.at(mapid).begin();
	      itcostheta!=HitMapVec.at(mapid).end();
	      itcostheta++)
	    for(std::map<int, int>::const_iterator itphi=itcostheta->second.begin();
		itphi!=itcostheta->second.end();
		itphi++)
	      {
		if(itphi->second > MaxValues[mapid])
		  {
		    MaxValues[mapid]       = itphi->second;
		    MaxCosThetaBins[mapid]  = itcostheta->first;
		    MaxPhiBins[mapid]       = itphi->first;
		  }
	      }   
	}
            
      //      std::cout<<"Winning seed strengths : " <<float(MaxValues[0])/PointsInRange.size()<<" " << float(MaxValues[1])/PointsInRange.size()
      //	       <<float(MaxValues[2])/PointsInRange.size()<<" " << float(MaxValues[3])/PointsInRange.size()<<std::endl;

      double theta, phi;
      if(MaxValues[0] > int(fMinDirectionStrength * PointsInRange.size()))
	{
	  theta = acos(2.0*float(MaxCosThetaBins[0])/CosThetaBins-1.0);
	  phi = float(MaxPhiBins[0])/PhiBins*6.284;
	}
      else if(MaxValues[1] > int(fMinDirectionStrength * PointsInRange.size()))
	{
	  theta = acos(2.0*float(0.5+MaxCosThetaBins[1])/CosThetaBins-1.0);
	  phi = float(MaxPhiBins[1])/PhiBins*6.284;
	}
      else if(MaxValues[2] > int(fMinDirectionStrength * PointsInRange.size()))
	{
	  theta = acos(2.0*float(0.5+MaxCosThetaBins[2])/CosThetaBins-1.0);
	  phi = float(0.5+MaxPhiBins[2])/PhiBins*6.284;
	}
      else if(MaxValues[3] > int(fMinDirectionStrength * PointsInRange.size()))
	{
	  theta = acos(2.0*float(0.5+MaxCosThetaBins[3])/CosThetaBins-1.0);
	  phi = float(0.5+MaxPhiBins[3])/PhiBins*6.284;
	}
      else
	{
	  return new recob::Seed();
	}
      
      TVector3 SeedDirection;
      SeedDirection.SetMagThetaPhi(fSeedLength, theta, phi);
      
      double PtArray[3], DirArray[3];
           
      // std::cout<<"SeedFinderService debug : Point " << 
      // 	HighestZPoint.X() << " " <<
      //	HighestZPoint.Y() << " " <<
      //	HighestZPoint.Z() << std::endl;
      
      PtArray[0]=CentreOfPoints.X();
      PtArray[1]=CentreOfPoints.Y();
      PtArray[2]=CentreOfPoints.Z();
      
      DirArray[0]=SeedDirection.X();
      DirArray[1]=SeedDirection.Y();
      DirArray[2]=SeedDirection.Z();
      
      ReturnSeed = new recob::Seed(PtArray,DirArray);
    }
  else ReturnSeed = new recob::Seed();


  //  std::cout<<"This is the seed finding routine" <<std::endl;
  //  std::cout<<"Called on a vector of space points of size " <<Points.size()<<std::endl;

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
  
  
}
