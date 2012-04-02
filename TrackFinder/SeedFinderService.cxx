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

void SeedFinderService::SetEventID(int EventID)
{
  fEventID=EventID;
}


recob::Seed SeedFinderService::FindSeedExhaustively(std::vector<recob::SpacePoint> Points)
{
  recob::Seed TrackSeedThisCombo;
  bool KeepChopping=true;
  while(KeepChopping)
    {
      TrackSeedThisCombo = FindSeedAtEnd(Points);
      if(TrackSeedThisCombo.IsValid()) KeepChopping=false;
      else
        {
          // if we found an invalid or weak seed, chop off some points and try again
  if(Points.size() >= fMinPointsInSeed)
    Points.resize(Points.size()-fMinPointsInSeed);
  else  // or if there is nothing left to chop, leave the loop
    KeepChopping=false;
        }
    }
  return TrackSeedThisCombo;

}


recob::Seed SeedFinderService::FindSeedAtEnd(std::vector<recob::SpacePoint> Points)

{
  static int OldEventID = -1;
  static int CallsSoFar =  0;
  if(OldEventID!=fEventID)
    {
      CallsSoFar=0;
      OldEventID=fEventID;
    }
  else
    CallsSoFar++;

  recob::Seed ReturnSeed;

  std::vector<int> PointsInRange;

  float HighestZ = Points.at(Points.size()-1).XYZ()[2];
  float HighestZX = Points.at(Points.size()-1).XYZ()[0];
  float HighestZY = Points.at(Points.size()-1).XYZ()[1];

  for(int index=Points.size()-1; index!=-1; index--)
    {
      // first check z, then check total distance
      //  (much faster, since most will be out of range in z anyway)
      if( ( HighestZ - Points.at(index).XYZ()[2] ) < fSeedLength)
        {
          if( ( ( pow(HighestZ-Points.at(index).XYZ()[2],2) +
                  pow(HighestZX-Points.at(index).XYZ()[0],2) +
                  pow(HighestZY-Points.at(index).XYZ()[1],2) ) < pow(fSeedLength,2)))
            PointsInRange.push_back(index);
        }
    }


  art::ServiceHandle<art::TFileService> tfs;
  std::stringstream HistoName("");

  std::map<int,std::map<int, int> > HitMap;
  HitMap.clear();
  int CosThetaBins=20;
  int PhiBins=20;

  HistoName<<"SeedThetaPhi"<<fEventID<<"c"<<CallsSoFar;
  TH2D * hThetaPhi = (TH2D*)tfs->make<TH2D>(HistoName.str().c_str(),HistoName.str().c_str(),CosThetaBins,-1,1,PhiBins,0,6.28);
  HistoName.flush();



  TVector3 HighestZPoint(HighestZX,HighestZY,HighestZ);

  if(PointsInRange.size()>fMinPointsInSeed)
    {

      for(unsigned int i=0; i!=PointsInRange.size(); i++)
        {
          TVector3 ThisPoint(Points.at(PointsInRange.at(i)).XYZ()[0],
                             Points.at(PointsInRange.at(i)).XYZ()[1],
                             Points.at(PointsInRange.at(i)).XYZ()[2]);
          TVector3 Step = HighestZPoint-ThisPoint;

          float costheta = Step.CosTheta();
          float phi = Step.Phi();

          int costhetabin = (int) ((costheta + 1.)*CosThetaBins/2.-0.5);
          int phibin = (int) ((phi/6.284)*PhiBins-0.5);

          HitMap[costhetabin][phibin] = HitMap[costhetabin][phibin] + 1;

          hThetaPhi->Fill(costheta,phi);
        }

      int MaxValue=0;
      int MaxCosThetaBin=0;
      int MaxPhiBin=0;

      // find most populated bin in hitmap
      for(std::map<int,std::map<int, int> >::const_iterator itcostheta=HitMap.begin();
          itcostheta!=HitMap.end();
          itcostheta++)
	for(std::map<int, int>::const_iterator itphi=itcostheta->second.begin();
            itphi!=itcostheta->second.end();
            itphi++)
          {
            if(itphi->second > MaxValue)
              {
                MaxValue       = itphi->second;
                MaxCosThetaBin = itcostheta->first;
                MaxPhiBin      = itphi->first;
              }
          }

      if(MaxValue > int(fMinDirectionStrength * PointsInRange.size()))
	{
          TVector3 SeedDirection;
          SeedDirection.SetMagThetaPhi(fSeedLength, acos( 2.0 * float(MaxCosThetaBin) / CosThetaBins - 1.0), float(MaxPhiBin)/PhiBins * 6.284);

	  double PtArray[3], DirArray[3];
	 
	  HighestZPoint.GetXYZ(    PtArray );
	  SeedDirection.GetXYZ(    DirArray );

          ReturnSeed.SetPoint(     PtArray );
          ReturnSeed.SetDirection( DirArray );
        }
      else ReturnSeed.SetValidity(false);
    }
  else
    ReturnSeed.SetValidity(false);



  std::cout<<"This is the seed finding routine" <<std::endl;
  std::cout<<"Called on a vector of space points of size " <<Points.size()<<std::endl;

  return ReturnSeed;
}

// This method will probably become deprecated, which is why it has 
// not been tidied up properly.

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
