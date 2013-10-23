//
// Name: SeedFinderAlgorithm.cxx
//
// Purpose: Implementation file for module SeedFinderAlgorithm.
//
// Ben Jones, MIT
//

#include <vector>
#include <stdint.h>

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "RecoAlg/SeedFinderAlgorithm.h"
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TPrincipal.h"
#include "TTree.h"
#include "art/Framework/Services/Optional/TFileService.h" 


namespace trkf {

  //----------------------------------------------------------------------------
  SeedFinderAlgorithm::SeedFinderAlgorithm(const fhicl::ParameterSet& pset)
  {
    

    reconfigure(pset);
  
    CalculateGeometricalElements();
    
    art::ServiceHandle<art::TFileService> tfs;
    ftMonitoringTree = (TTree*)tfs->make<TTree>("SeedTree","SeedTree");
 
    ftMonitoringTree->Branch("ThetaXZ",   &ftThetaXZ,   "ThetaXZ/F");
    ftMonitoringTree->Branch("ThetaYZ",   &ftThetaYZ,   "ThetaYZ/F");
    ftMonitoringTree->Branch("Theta",     &ftTheta,     "Theta/F");
     ftMonitoringTree->Branch("NSpts",     &ftNSpts,     "NSpts/I");
    ftMonitoringTree->Branch("NUHits",    &ftNUHits,    "NUHits/I");
    ftMonitoringTree->Branch("NVHits",    &ftNVHits,    "NVHits/I");
    ftMonitoringTree->Branch("NWHits",    &ftNWHits,    "NWHits/I");
    ftMonitoringTree->Branch("URMSb",     &ftURMSb,     "URMSb/F");
    ftMonitoringTree->Branch("VRMSb",     &ftVRMSb,     "VRMSb/F");
    ftMonitoringTree->Branch("WRMSb",     &ftWRMSb,     "WRMSb/F");
    ftMonitoringTree->Branch("URMS",      &ftURMS,      "URMS/F");
    ftMonitoringTree->Branch("VRMS",      &ftVRMS,      "VRMS/F");
    ftMonitoringTree->Branch("WRMS",      &ftWRMS,      "WRMS/F");
    ftMonitoringTree->Branch("Keep",      &ftKeep,      "Keep/B");

    

  

  }

  //----------------------------------------------------------------------------
  SeedFinderAlgorithm::~SeedFinderAlgorithm()
  {
  }

  //----------------------------------------------------------------------------
  void SeedFinderAlgorithm::reconfigure(fhicl::ParameterSet const& pset)
  {
 
    fSptalg                = new trkf::SpacePointAlg(pset.get<fhicl::ParameterSet>("SpacePointAlg"));
    
    fInitSeedLength        = pset.get<double>("InitSeedLength");
    fMinPointsInSeed       = pset.get<int>("MinPointsInSeed");
  
    fRefits                = pset.get<double>("Refits");

    fHitResolution         = pset.get<double>("HitResolution");

    fOccupancyCut          = pset.get<double>("OccupancyCut");
    fExtendSeeds           = pset.get<bool>("ExtendSeeds");
	

    fMaxViewRMS.resize(3);
    fMaxViewRMS            = pset.get<std::vector<double> >("MaxViewRMS"); 
 
    
    //    fPCAThreshold          = pset.get<double>("PCAThreshold");
    // Deprecated (but left in, in case we want to bring it back)

  }



  
  

  //------------------------------------------------------------
  // Given a set of spacepoints, find seeds, and catalogue
  //  spacepoints by the seeds they formed
  //
  std::vector<recob::Seed> SeedFinderAlgorithm::FindSeeds( art::PtrVector<recob::Hit> const& HitsFlat, std::vector<art::PtrVector<recob::Hit> >& CataloguedHits)
  {


    // Vector of seeds found to return
    std::vector<recob::Seed>       ReturnVector;
    
    // First check if there is any overlap
    std::vector<recob::SpacePoint> spts;
    
    fSptalg->clearHitMap();
    fSptalg->makeSpacePoints(HitsFlat, spts);

    if(int(spts.size()) < fMinPointsInSeed)
      {
	return std::vector<recob::Seed>();
      }


    // If we got here, we have enough spacepoints to potentially make seeds from these hits.

    // This is table will let us quickly look up which hits are in a given view / channel.
    //  structure is OrgHits[View][Channel] = {index1,index2...}
    // where the indices are of HitsFlat[index].
 
    std::map<geo::View_t, std::map<uint32_t, std::vector<int> > > OrgHits;
  
    // These two tables contain the hit to spacepoint mappings
    
    std::vector<std::vector<int> > SpacePointsPerHit(HitsFlat.size(), std::vector<int>());
    std::vector<std::vector<int> > HitsPerSpacePoint(spts.size(), std::vector<int>());

    // This vector records the status of each hit.
    // The possible values are:
    //  0. hit unused
    //  1. hit used in a seed
    // Initially none are used.

    std::vector<char> HitStatus(HitsFlat.size(),0);

    
    // Fill the organizing map

    for(size_t i=0; i!=HitsFlat.size(); ++i)
      {
	OrgHits[HitsFlat.at(i)->View()][HitsFlat.at(i)->Channel()].push_back(i);
      }


    // Fill the spacepoint / hit lookups

    for(size_t iSP =0; iSP!=spts.size(); ++iSP)
      {
	art::PtrVector<recob::Hit> HitsThisSP = fSptalg->getAssociatedHits(spts.at(iSP));

	for(size_t iH = 0; iH != HitsThisSP.size(); ++iH)
	  {
	    geo::View_t ThisView     = HitsThisSP.at(iH)->View();
	    uint32_t    ThisChannel  = HitsThisSP.at(iH)->Channel();
	    float       ThisTime     = HitsThisSP.at(iH)->PeakTime();
	    float       eta          = 0.001;


	    for(size_t iOrg = 0; iOrg!= OrgHits[ThisView][ThisChannel].size(); ++iOrg)
	      {
		if(fabs(ThisTime - HitsFlat.at( OrgHits[ThisView][ThisChannel][iOrg] )->PeakTime() ) < eta)
		  {
		    SpacePointsPerHit.at(OrgHits[ThisView][ThisChannel][iOrg]).push_back(iSP);
		    HitsPerSpacePoint.at(iSP).push_back(OrgHits[ThisView][ThisChannel][iOrg]);
		  }
	      }
	  }
      }


    // The general idea will be:
    //  A. Make spacepoints from remaining hits
    //  B. Look for 1 seed in that vector
    //  C. Discard used hits
    //  D. Repeat




    // This vector keeps track of the status of each space point.  
    // The key is the position in the AllSpacePoints vector.
    // The value is 
    //    0: point unused
    //    1: point used in seed (no longer implemented)
    //    2: point thrown but unused 
    //    3: flag to terminate seed finding
    
    std::vector<char> PointStatus(spts.size(),0);
    
    std::vector<std::map<geo::View_t, std::vector<int> > > WhichHitsPerSeed;
    

    bool KeepChopping = true;
    
    while(KeepChopping)
      {
	
	
	// This vector keeps a list of the points used in this seed
	std::vector<int>                PointsUsed;
	    
	// Find exactly one seed, starting at high Z
	recob::Seed TheSeed = FindSeedAtEnd(spts, PointStatus, PointsUsed);
	
	// If it was a good seed, collect up the relevant spacepoints
	// and add the seed to the return vector 
	    
	
	if(TheSeed.IsValid())
	  {
	    
	    for(size_t iP=0; iP!=PointsUsed.size(); ++iP)
	      {
		for(size_t iH=0; iH!=HitsPerSpacePoint.at(PointsUsed.at(iP)).size(); ++iH)
		  {
		    int UsedHitID = HitsPerSpacePoint.at(PointsUsed.at(iP)).at(iH);
		    HitStatus[UsedHitID] = 2;
		  }
	      }
	    PointStatus[PointsUsed.at(0)] = 1;
	    ConsolidateSeed(TheSeed, HitsFlat, HitStatus, OrgHits, fExtendSeeds);
	    
	  }
	
	if(TheSeed.IsValid())
	  {
	    if(fRefits>0)
	      {
		std::vector<char> HitStatusGood;
		recob::Seed      SeedGood;
		for(size_t r=0; r!=uint(fRefits); ++r)
		  {
		    SeedGood =      TheSeed;
		    HitStatusGood = HitStatus;

		    std::vector<int> PresentHitList;
		    for(size_t iH=0; iH!=HitStatus.size(); ++iH)
		      {
			if(HitStatus[iH]==2)
			  {
			    PresentHitList.push_back(iH);
			  }
		      }
		    double pt[3], dir[3], err[3];
		    
		    TheSeed.GetPoint(pt,err);
		    TheSeed.GetDirection(dir,err);
		    
		    
		    TVector3 Center, Direction;
		    GetCenterAndDirection(HitsFlat, PresentHitList, Center, Direction);
		    
		    
		    Direction = Direction.Unit() * TheSeed.GetLength();
		    
		    
		    
		    for(size_t n=0; n!=3; ++n)
		      {
			pt[n] = Center[n];
			dir[n] = Direction[n];
			
			
			TheSeed.SetPoint(pt, err);
			TheSeed.SetDirection(dir, err);
		      }
		    
		    ConsolidateSeed(TheSeed, HitsFlat, HitStatus, OrgHits, fExtendSeeds);	  
		    
		    // if we accidentally invalidated the seed, go back to the old one and escape
		    if(!TheSeed.IsValid()) 
		      {
			HitStatus = HitStatusGood; 
			TheSeed   = SeedGood;
			break;
		      }
		  }
	      } 
	  }
	
      
	
	if(TheSeed.IsValid())
	  {
	    WhichHitsPerSeed.push_back(std::map<geo::View_t, std::vector<int> >());
	    
	    art::PtrVector<recob::Hit> HitsWithThisSeed;       
	    for(size_t iH=0; iH!=HitStatus.size(); ++iH)
	      {
		if(HitStatus.at(iH)==2)
		  {
		    WhichHitsPerSeed.at(WhichHitsPerSeed.size()-1)[HitsFlat[iH]->View()].push_back(iH);
		    HitsWithThisSeed.push_back(HitsFlat.at(iH));
		    HitStatus.at(iH)=1;
		    
		    for(size_t iSP=0; iSP!=SpacePointsPerHit.at(iH).size(); ++iSP)
		      {
			PointStatus[SpacePointsPerHit.at(iH).at(iSP)] = 1;
		      }
		  }
	      }
	    
	    
	    // Record that we used this set of hits with this seed in the return catalogue
	    ReturnVector.push_back(TheSeed);
	    CataloguedHits.push_back(HitsWithThisSeed);
	    
	    
	    // Tidy up
	    HitsWithThisSeed.clear();
	  }
	else
	  {
	    
	    // If it was not a good seed, throw out the top SP and try again
	    PointStatus.at(PointsUsed.at(0))=2;
	  }
	
	int TotalSPsUsed=0;
	for(size_t i=0; i!=PointStatus.size(); ++i)
	  {
	    if(PointStatus[i]!=0) TotalSPsUsed++;
	  }
	    
	if((int(spts.size()) - TotalSPsUsed) < fMinPointsInSeed)
	  KeepChopping=false;
	
	if((PointStatus[0]==3)||(PointStatus.size()==0)) KeepChopping=false;
	
	PointsUsed.clear();

      } // end spt-level loop
    

    // If we didn't find any seeds, we can make one last ditch attempt. See
    //  if the whole collection is colinear enough to make one megaseed
    //  (good for patching up high angle tracks)


    if(ReturnVector.size()==0)
      {
	std::vector<int> ListAllPoints;
	for(size_t i=0; i!=spts.size(); ++i)
	  ListAllPoints.push_back(i);
	
	TVector3 SeedCenter(    0, 0, 0 );
	TVector3 SeedDirection( 0, 0, 0 );
	double   SeedStrength = 0;

	GetCenterAndDirection(spts, ListAllPoints, SeedCenter, SeedDirection);


	if(SeedStrength > fPCAThreshold)
	  {
	    double PtArray[3], DirArray[3];
	    for(size_t n=0; n!=3; ++n)
	      {
		PtArray[n] = SeedCenter[n];
		DirArray[n] = SeedDirection[n];
	      }
	    recob::Seed TheSeed(PtArray,DirArray);
	
	    //	    if(fRefits>0) RefitSeed(TheSeed, spts);	
	    
	    std::vector<double> RMS = GetHitRMS(TheSeed, spts);
	    
	    bool ThrowOutSeed = false;
	    if(fMaxViewRMS.at(0)>0)
	      {
		for(size_t j=0; j!=fMaxViewRMS.size(); j++)
		  {
		    if(fMaxViewRMS.at(j)<RMS.at(j))
		      {
			ThrowOutSeed=true;
		      }
		  }
	      }
	    if(!ThrowOutSeed)
	      {
		ReturnVector.push_back(TheSeed);
		CataloguedHits.push_back(HitsFlat);
		
	      }
	  
	  }
      }
    
    // Tidy up
    SpacePointsPerHit.clear();
    HitsPerSpacePoint.clear();
    PointStatus.clear();	  
    OrgHits.clear();
    HitStatus.clear();
    
    

    return ReturnVector;


  }



  //------------------------------------------------------------
  // Latest extendseed method
  //

  void SeedFinderAlgorithm::ConsolidateSeed(recob::Seed& TheSeed, art::PtrVector<recob::Hit> const& HitsFlat, std::vector<char>& HitStatus,
					    std::map<geo::View_t, std::map<uint32_t, std::vector<int> > >& OrgHits, bool Extend)
  {
    
    
    
    bool ThrowOutSeed = false;
    
    // This will keep track of what hits are in this seed
    std::map<geo::View_t, std::map<uint32_t, std::vector<int> > > HitsInThisSeed;
    
    int NHitsThisSeed=0;

    double MinS = 1000, MaxS=-1000;
    for(size_t i=0; i!=HitStatus.size(); ++i)
      {
	if(HitStatus.at(i)==2)
	  {
	    double disp, s;
	    GetHitDistAndProj(TheSeed, HitsFlat.at(i),disp, s);
	    if(fabs(s)>1.2)
	      {
		// This hit is not rightfully part of this seed, toss it.
		HitStatus[i]=0;
	      }
	    else
	      {
		NHitsThisSeed++;
			
		if(s<MinS) MinS = s;
		if(s>MaxS) MaxS = s;
		HitsInThisSeed[HitsFlat.at(i)->View()][HitsFlat.at(i)->Channel()].push_back(i);
	      }
	  }
      }

    double LengthRescale = (MaxS - MinS)/2.;
    double PositionShift = (MaxS + MinS)/2.;

    double pt[3], dir[3], err[3];
    TheSeed.GetPoint(pt, err);
    TheSeed.GetDirection(dir, err);
    

    for(size_t n=0; n!=3; ++n)
      {
	pt[n]  += dir[n] * PositionShift;
	dir[n] *= LengthRescale;
      }

    TheSeed.SetPoint(pt, err);
    TheSeed.SetDirection(dir, err);

    
    // Run through checking if we missed any hits
    for(auto itP = HitsInThisSeed.begin(); itP!=HitsInThisSeed.end(); ++itP)
      {
	double dist, s;
	geo::View_t View = itP->first;
	uint32_t LowestChan = itP->second.begin()->first;
	uint32_t HighestChan = itP->second.rbegin()->first;
	for(uint32_t c=LowestChan; c!=HighestChan; ++c)
	  {
	    for(size_t h=0; h!=OrgHits[View][c].size(); ++h)
	      {
		if(HitStatus[OrgHits[View][c].at(h)]==0)
		  {
		    GetHitDistAndProj(TheSeed, HitsFlat[OrgHits[View][c].at(h)], dist, s);
		    if(dist < fHitResolution)
		      {
			NHitsThisSeed++;
			
			HitStatus[OrgHits[View][c].at(h)]=2;
			HitsInThisSeed[View][c].push_back(OrgHits[View][c].at(h));
		      }
		  }
	      }
	  }
      }

    if(NHitsThisSeed==0) ThrowOutSeed=true;

    
    // Check seed occupancy
    
    uint32_t LowestChanInSeed[3], HighestChanInSeed[3], LowestChanInView[3], HighestChanInView[3];
    double Occupancy[3];
    
    for(auto itP = HitsInThisSeed.begin(); itP!=HitsInThisSeed.end(); ++itP)
      {
	int ViewID;
	geo::View_t View = itP->first;
	if(View==geo::kU) ViewID=0;
	if(View==geo::kV) ViewID=1;
	if(View==geo::kW) ViewID=2;

        LowestChanInSeed[ViewID]  = itP->second.begin()->first;
	HighestChanInSeed[ViewID] = itP->second.rbegin()->first;
	LowestChanInView[ViewID]  =  OrgHits[View].begin()->first;
	HighestChanInView[ViewID] = OrgHits[View].rbegin()->first;
	
	int FilledChanCount=0;
	
	for(size_t c  = LowestChanInSeed[ViewID]; c!=HighestChanInSeed[ViewID]; ++c)
	  {
	    if(itP->second[c].size()>0) ++FilledChanCount;
	  }
	
	Occupancy[ViewID] = float(FilledChanCount) / float(HighestChanInSeed[ViewID]-LowestChanInSeed[ViewID]);
      }	
    
    
    for(size_t n=0; n!=3; ++n)
      if(Occupancy[n]<fOccupancyCut) ThrowOutSeed=true;
	 
    if((Extend)&&(!ThrowOutSeed))
      {
	std::vector<std::vector<double> > ToAddNegativeS(3,std::vector<double>());
	std::vector<std::vector<double> > ToAddPositiveS(3,std::vector<double>());
	std::vector<std::vector<int> >    ToAddNegativeH(3,std::vector<int>());
	std::vector<std::vector<int> >    ToAddPositiveH(3,std::vector<int>());
	
	for(auto itP = HitsInThisSeed.begin(); itP!=HitsInThisSeed.end(); ++itP)
	  {
	    double dist, s;
	    int ViewID=0;

	    geo::View_t View = itP->first;

	    if(View==geo::kU) ViewID=0;
	    if(View==geo::kV) ViewID=1;
	    if(View==geo::kW) ViewID=2;
	    
	    
	    if(LowestChanInSeed[ViewID]>LowestChanInView[ViewID])
	      {
		for(uint32_t c=LowestChanInSeed[ViewID]-1; c!=LowestChanInView[ViewID]; --c)
		  {
		    bool GotOneThisChannel=false;
		    for(size_t h=0; h!=OrgHits[View][c].size(); ++h)
		      {
			GetHitDistAndProj(TheSeed, HitsFlat[OrgHits[View][c].at(h)], dist, s);
			if(dist < fHitResolution)
			  {
			    GotOneThisChannel=true;
			    if(s<0) 
			      {
				ToAddNegativeS[ViewID].push_back(s);
				ToAddNegativeH[ViewID].push_back(OrgHits[View][c].at(h));
			      }
			    else
			      {
				ToAddPositiveS[ViewID].push_back(s);
				ToAddPositiveH[ViewID].push_back(OrgHits[View][c].at(h));
			      }
			  }
		      }
		    if(GotOneThisChannel==false) break;
		  }
	      }
	    if(HighestChanInSeed < HighestChanInView)
	      
	      for(uint32_t c=HighestChanInSeed[ViewID]+1; c!=HighestChanInView[ViewID]; ++c)
		{
		  bool GotOneThisChannel=false;
		  for(size_t h=0; h!=OrgHits[View][c].size(); ++h)
		    {
		      GetHitDistAndProj(TheSeed, HitsFlat[OrgHits[View][c].at(h)], dist, s);
		      if(dist < fHitResolution)
			{
			  GotOneThisChannel=true;
			  if(s<0) 
			    {
			      ToAddNegativeS[ViewID].push_back(s);
			      ToAddNegativeH[ViewID].push_back(OrgHits[View][c].at(h));
			    }
			  else
			    {
			      ToAddPositiveS[ViewID].push_back(s);
			      ToAddPositiveH[ViewID].push_back(OrgHits[View][c].at(h));
			    }
			}
		    }
		  if(GotOneThisChannel==false) break;
		}
	  }
	
	
	double ExtendPositiveS, ExtendNegativeS;

	if((ToAddPositiveS[0].size()>0)&&
	   (ToAddPositiveS[1].size()>0)&&
	   (ToAddPositiveS[2].size()>0))
	  {
	    for(size_t n=0; n!=3; ++n)
	      {
		int n1 = (n+1)%3;
		int n2 = (n+2)%3;
		
		if( (ToAddPositiveS[n].back() <= ToAddPositiveS[n1].back()) &&
		    (ToAddPositiveS[n].back() <= ToAddPositiveS[n2].back()) )
		  {
		    ExtendPositiveS = ToAddPositiveS[n].back();
		  }
	      }
	  }

	if((ToAddNegativeS[0].size()>0)&&
	   (ToAddNegativeS[1].size()>0)&&
	   (ToAddNegativeS[2].size()>0))
	  {
	    for(size_t n=0; n!=3; ++n)
	      {
		int n1 = (n+1)%3;
		int n2 = (n+2)%3;
		if( (ToAddNegativeS[n].back() >= ToAddNegativeS[n1].back()) &&
		    (ToAddNegativeS[n].back() >= ToAddNegativeS[n2].back()) )
		  {
		    ExtendNegativeS = ToAddNegativeS[n].back();
		  }
	      }
	  }
	
	if(fabs(ExtendNegativeS) < 1.) ExtendNegativeS=-1.;
	if(fabs(ExtendPositiveS) < 1.) ExtendPositiveS=1.;

	
	LengthRescale = (ExtendPositiveS - ExtendNegativeS)/2.;
	PositionShift = (ExtendPositiveS + ExtendNegativeS)/2.;
	
	for(size_t n=0; n!=3; ++n)
	  {
	    pt[n]  += dir[n] * PositionShift;
	    dir[n] *= LengthRescale;
	  }


	TheSeed.SetPoint(pt, err);
	TheSeed.SetDirection(dir, err);
	
      }


    
    if(ThrowOutSeed) TheSeed.SetValidity(false);
      
    
  }
  
  //------------------------------------------------------------


  void SeedFinderAlgorithm::GetHitDistAndProj( recob::Seed const& ASeed,  art::Ptr<recob::Hit> const& AHit, double& disp, double& s)
  {
    art::ServiceHandle<util::DetectorProperties> det;
    art::ServiceHandle<geo::Geometry> geom;
    
    double xyzStart[3], xyzEnd[3];
    
    geom->WireEndPoints(0,0, AHit->WireID().Plane, AHit->WireID().Wire, xyzStart, xyzEnd);
    
    double HitX = det->ConvertTicksToX(AHit->PeakTime(), AHit->WireID().Plane, 0, 0);
    
    double HitXHigh = det->ConvertTicksToX(AHit->EndTime(), AHit->WireID().Plane, 0, 0);    
    double HitXLow =  det->ConvertTicksToX(AHit->StartTime(), AHit->WireID().Plane, 0, 0);
    
    double HitWidth = HitXHigh - HitXLow;
    
    double pt[3], dir[3], err[3];
    
    ASeed.GetDirection(dir,err);
    ASeed.GetPoint(pt,err);
    


    TVector3 sPt  (pt[0],   pt[1],                 pt[2]);
    TVector3 sDir (dir[0],  dir[1],                dir[2]);
    TVector3 hPt   (HitX,    xyzStart[1],           xyzStart[2]);
    TVector3 hDir  (0,       xyzStart[1]-xyzEnd[1], xyzStart[2]-xyzEnd[2]);

    s = (sPt - hPt).Dot(hDir*(hDir.Dot(sDir))-sDir*(hDir.Dot(hDir))) / (hDir.Dot(hDir)*sDir.Dot(sDir)-pow(hDir.Dot(sDir),2));
    
    disp = fabs((sPt - hPt).Dot(sDir.Cross(hDir))/(sDir.Cross(hDir)).Mag()) / HitWidth;
  }

  //------------------------------------------------------------
  // Try to find one seed at the high Z end of a set of spacepoints 
  //

    recob::Seed  SeedFinderAlgorithm::FindSeedAtEnd(std::vector<recob::SpacePoint> const& Points, std::vector<char>& PointStatus, std::vector<int>& PointsInRange)
  {
    // This pointer will be returned later
    recob::Seed ReturnSeed;
    
    // Keep track of spacepoints we used, not just their IDs
    std::vector<recob::SpacePoint> PointsUsed;
   
    // Clear output vector
    PointsInRange.clear();
   
    // Loop through hits looking for highest Z seedable point
    TVector3 HighestZPoint;
    bool NoPointFound=true;
    int counter = Points.size()-1; 
    while((NoPointFound==true)&&(counter>=0))
      {
	if(PointStatus[counter]==0)
	  {
	    HighestZPoint = TVector3(Points.at(counter).XYZ()[0],
				     Points.at(counter).XYZ()[1],
				     Points.at(counter).XYZ()[2]);
	    NoPointFound=false;
	  }
	else
	  counter--;
      }
    if(NoPointFound)
      {
	// We didn't find a high point at all
	//  - let the algorithm know to give up.
	PointStatus[0]=3;
      }

    // Now we have the high Z point, loop through collecting
    // near enough hits.  We look 2 seed lengths away, since 
    // the seed is bidirectional from the central point 
       
    double TwiceLength = 2.0*fInitSeedLength;
    
    for(int index=Points.size()-1; index!=-1; --index)
      {
	if(PointStatus[index]==0)
	  {
	    // first check z, then check total distance
	    //  (much faster, since most will be out of range in z anyway)
	    if( ( HighestZPoint[2] - Points.at(index).XYZ()[2] ) < TwiceLength)
	      {
		double DistanceToHighZ =pow(
					     pow(HighestZPoint[1]-Points.at(index).XYZ()[1],2) +
					     pow(HighestZPoint[2]-Points.at(index).XYZ()[2],2),0.5 ); 
		if( DistanceToHighZ < TwiceLength)
		  {
		    PointsInRange.push_back(index);
		    PointsUsed.push_back(Points.at(index));
		  }
	      }
	    else break;
	  }
      }
    
    TVector3 SeedCenter(    0, 0, 0 );
    TVector3 SeedDirection( 0, 0, 0 );
    double   SeedStrength = 0;
    

    // Check we have enough points in here to form a seed,
    // otherwise return a dud
    int NPoints = PointsInRange.size();
    ftNSpts = NPoints;
   
    if(NPoints<fMinPointsInSeed) return  recob::Seed();
    
    GetCenterAndDirection(Points, PointsInRange, SeedCenter, SeedDirection);
    
    // See if seed points have some linearity

    bool ThrowOutSeed = false;
     
   
    double PtArray[3], DirArray[3];
    
    PtArray[0] = SeedCenter.X();
    PtArray[1] = SeedCenter.Y();
    PtArray[2] = SeedCenter.Z();
    
    
    double AngleFactor = pow(pow(SeedDirection.Y(),2)+pow(SeedDirection.Z(),2),0.5)/SeedDirection.Mag();
    
    DirArray[0] = SeedDirection.X() * fInitSeedLength / AngleFactor;
    DirArray[1] = SeedDirection.Y() * fInitSeedLength / AngleFactor;
    DirArray[2] = SeedDirection.Z() * fInitSeedLength / AngleFactor;
    
    
    ReturnSeed = recob::Seed(PtArray,DirArray);
   
    if(SeedStrength < fPCAThreshold)
      {
	ThrowOutSeed=true;
      }
	

    std::vector<double> RMSb = GetHitRMS(ReturnSeed, PointsUsed);
    
    ftURMSb = RMSb.at(0);
    ftVRMSb = RMSb.at(1);
    ftWRMSb = RMSb.at(2);

   
    //    if((!ThrowOutSeed) && (fRefits>0))
    //      RefitSeed(ReturnSeed, PointsUsed);
    
	
    std::vector<double> RMS = GetHitRMS(ReturnSeed, PointsUsed);
    
    ftURMS = RMS.at(0);
    ftVRMS = RMS.at(1);
    ftWRMS = RMS.at(2);
    
    if(fMaxViewRMS.at(0)>0)
      {

	for(size_t j=0; j!=fMaxViewRMS.size(); j++)
	  {
	    if(fMaxViewRMS.at(j)<RMS.at(j)) 
	      {
	    
		ThrowOutSeed=true;
	      }
	    //   mf::LogVerbatim("SeedFinderAlgorithm") << RMS.at(j);		
	  }
      }

    ftKeep = ThrowOutSeed;
    ftMonitoringTree->Fill();

    // If the seed is marked as bad, return a dud, otherwise
    //  return the ReturnSeed pointer
    if(!ThrowOutSeed)
      return ReturnSeed;
    else
      return recob::Seed();
  }








  //-----------------------------------------------------
  // Given a seed, find which spacepoints from a collection
  //  are within some distance
  //

  std::vector<int> SeedFinderAlgorithm::DetermineNearbySPs(recob::Seed const& TheSeed, std::vector<recob::SpacePoint> const& AllSpacePoints, std::map<int, int> PointStatus, double MaxDistance)
  {
    std::vector<int> ReturnVector;

    double SeedPt[3], Err[3];
    TheSeed.GetPoint( SeedPt,  Err );
    
   
    
    double SeedLength = TheSeed.GetLength();

    for(int i=AllSpacePoints.size()-1; i!=-1; --i)
      {
	if(PointStatus[i]!=1)
	  {
	    // check Z separation first - helps prevent doing unnecessary calculations

	    if(( SeedPt[2]-AllSpacePoints.at(i).XYZ()[2] ) > SeedLength) 
	      break;
	    
 	    else if( TheSeed.GetDistanceFrom(AllSpacePoints.at(i))
		     < MaxDistance )
	      ReturnVector.push_back(i);
	    
	  }
	
      }

    return ReturnVector;
  }



  //-----------------------------------------------------
  // Find the RMS distance of the hits in a set of spacepoints
  //  from the spine of a seed in each view.
  //

  std::vector<double> SeedFinderAlgorithm::GetHitRMS(recob::Seed const& TheSeed, std::vector<recob::SpacePoint> const& SpacePoints)
  {
    // Vector to return later, of one RMS value for each view
    std::vector<double>    ReturnVec;

    // Services we need
    art::ServiceHandle<geo::Geometry>                geom;
    art::ServiceHandle<util::DetectorProperties>     det;

    // At present, algorithm works only for 1 TPC, 1 cryostat
    unsigned int c=0, t=0;

    std::map<int,art::PtrVector<recob::Hit> > HitMap;
    std::map<uint32_t, bool>                  HitsClaimed;
    size_t Planes = geom->TPC(0).Nplanes();

    // for each spacepoint
    for(std::vector<recob::SpacePoint>::const_iterator itSP=SpacePoints.begin();
        itSP!=SpacePoints.end(); itSP++)
      {
	art::PtrVector<recob::Hit> HitsThisSP = fSptalg->getAssociatedHits(*itSP);
	for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisSP.begin();
	    itHit!=HitsThisSP.end(); itHit++)
	  {
	    if(!HitsClaimed[(*itHit)->Channel()])
	      {

		HitMap[(*itHit)->View()].push_back(*itHit);
		HitsClaimed[(*itHit)->Channel()]=true;
	      }
	    
	  }
	HitsThisSP.clear();
      }
   
    // Find the RMS in each plane
    for(size_t planeno=0; planeno!=Planes; planeno++)
      {
	int view = geom->TPC(0).Plane(planeno).View();
	
	art::PtrVector<recob::Hit> HitsThisPlane = HitMap[view];
	double SeedCentralWire=0, SeedCentralTime=0;

	double SeedPt[3], Err[3], SeedDir[3];
	TheSeed.GetPoint(      SeedPt,   Err);
	TheSeed.GetDirection(  SeedDir,  Err);

	TVector3 SeedPoint(     SeedPt[0],     SeedPt[1],   SeedPt[2]);
        TVector3 SeedDirection( SeedDir[0],    SeedDir[1],  SeedDir[2]);

	// Determine useful directions for this plane
	double Wire0End1[3],Wire0End2[3];
	double Wire1End1[3],Wire1End2[3];


	geom->WireEndPoints(c,t,planeno,0, Wire0End1, Wire0End2);
	geom->WireEndPoints(c,t,planeno,1, Wire1End1, Wire1End2);

	TVector3 Wire0End1Vec(Wire0End1[0],Wire0End1[1],Wire0End1[2]);
	TVector3 Wire0End2Vec(Wire0End2[0],Wire0End2[1],Wire0End2[2]);
	TVector3 Wire1End1Vec(Wire1End1[0],Wire1End1[1],Wire1End1[2]);
	TVector3 Wire1End2Vec(Wire1End2[0],Wire1End2[1],Wire1End2[2]);

	TVector3 WireVec    = (Wire0End1Vec-Wire0End2Vec).Unit();
	TVector3 XVec(1,0,0);
	TVector3 PlaneNormDirection=(WireVec.Cross(XVec)).Unit();
	if( (Wire1End1Vec-Wire0End1Vec).Dot(PlaneNormDirection)<0)
	  PlaneNormDirection = -PlaneNormDirection;

	double   WirePitch  = geom->WirePitch();
	double   WireOffset = PlaneNormDirection.Dot(Wire1End1);

	SeedCentralWire = (SeedPoint.Dot(PlaneNormDirection) - WireOffset) / WirePitch;

        double SeedDirInPlane = SeedDirection.Dot(PlaneNormDirection);
        double SeedDirInTime  = SeedDirection.Dot(XVec);

        SeedCentralTime = det->ConvertXToTicks(SeedPt[0],planeno,t,c);

        double OneTickInX = det->GetXTicksCoefficient();

        double LSTot=0;
	int    Count=0;

	double a = SeedDirInTime / SeedDirInPlane;

	
	for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisPlane.begin(); itHit!=HitsThisPlane.end(); itHit++)
          {
            double DDisp = (double((*itHit)->WireID().Wire)-SeedCentralWire)*WirePitch/OneTickInX;
	    double Time = (*itHit)->PeakTime();
	    double Width = 0.5*(*itHit)->EndTime()-(*itHit)->StartTime();
	    double XDisp = (Time-SeedCentralTime);
	    LSTot += pow(XDisp - a * DDisp, 2)/pow(Width,2);
	    Count++;
	    
          }
	ReturnVec.push_back(pow(LSTot/Count,0.5));
      }

    // tidy up
    for(std::map<int,art::PtrVector<recob::Hit> >::iterator it = HitMap.begin();
	it!=HitMap.end(); ++it)
      it->second.clear();

    return ReturnVec;
  
  }


  //-----------------------------------------------------
  // Given a seed and some spacepoints, refit the seed
  //  to the hits in the SPs
  //

  void SeedFinderAlgorithm::RefitSeed(recob::Seed& TheSeed, std::vector<recob::SpacePoint> const& SpacePoints)
  {

    // Get the services we need
    art::ServiceHandle<geo::Geometry>            geom;
    art::ServiceHandle<util::DetectorProperties> det;


    if((geom->NTPC()!=1)||geom->NTPC()!=1)
      {
        throw cet::exception("SeedFinder : Refit only works for 1 tpc, 1 cryostat detector")<<
          "seed refitting feature not yet developped for multi cryostat or "<<
          "multi TPC detector - see TrackFinder/SeedFinder.cxx"<<std::endl;
      }
    size_t Planes = geom->TPC(0).Nplanes();

    // Get this hits in each view that made this seed

    std::map<int,art::PtrVector<recob::Hit> > HitMap;
    std::map<uint32_t, bool>                  HitsClaimed;


    // for each spacepoint
    for(std::vector<recob::SpacePoint>::const_iterator itSP=SpacePoints.begin();
        itSP!=SpacePoints.end(); itSP++)
      {
        // get hits from each plane (ensuring each used once only)
	art::PtrVector<recob::Hit> HitsThisSP = fSptalg->getAssociatedHits(*itSP);
	for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisSP.begin();
	    itHit!=HitsThisSP.end(); itHit++)
	  {
	    if(!HitsClaimed[(*itHit)->Channel()])
	      {
		HitMap[(*itHit)->View()].push_back(*itHit);
		HitsClaimed[(*itHit)->Channel()]=true;
	      }
	  }
      }

    double                  WirePitch  = geom->WirePitch();
    TVector3                XVec       = TVector3(1,0,0);

    std::vector<double>     WireOffsets;            WireOffsets.resize(Planes);
    std::vector<TVector3>   PlaneNormDirections;    PlaneNormDirections.resize(Planes);
    std::vector<TVector3>   WireVecs;               WireVecs.resize(Planes);
    std::vector<double>     HitCentralTimes;        HitCentralTimes.resize(Planes);
    std::vector<double>     HitCentralWires;        HitCentralWires.resize(Planes);
    std::vector<int>        NHits;                  NHits.resize(Planes);

    double AverageXWeighted;
    double AverageXDemocratic;
    int    TotalHits;

    // Calculate geometrical stuff for each plane
    for(size_t planeno=0; planeno!=Planes; planeno++)
      {
	// Determine useful directions for this plane
	int view = geom->TPC(0).Plane(planeno).View();

	double Wire0End1[3],Wire0End2[3];
	double Wire1End1[3],Wire1End2[3];
	geom->WireEndPoints(0,0,planeno,0, Wire0End1, Wire0End2);
	geom->WireEndPoints(0,0,planeno,1, Wire1End1, Wire1End2);
	
	TVector3 Wire0End1Vec(Wire0End1[0],Wire0End1[1],Wire0End1[2]);
	TVector3 Wire0End2Vec(Wire0End2[0],Wire0End2[1],Wire0End2[2]);
	TVector3 Wire1End1Vec(Wire1End1[0],Wire1End1[1],Wire1End1[2]);
	TVector3 Wire1End2Vec(Wire1End2[0],Wire1End2[1],Wire1End2[2]);
	
	WireVecs[view]            = (Wire0End1Vec-Wire0End2Vec).Unit();
        PlaneNormDirections[view] = (WireVecs[view].Cross(XVec)).Unit();
	if( (Wire1End1Vec-Wire0End1Vec).Dot(PlaneNormDirections[view])<0)
	  PlaneNormDirections[view] = -PlaneNormDirections[view];
	
	WireOffsets[view]         = PlaneNormDirections.at(view).Dot(Wire0End1Vec);
	
	art::PtrVector<recob::Hit> HitsThisPlane = HitMap[view]; 
	for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisPlane.begin(); itHit!=HitsThisPlane.end(); itHit++)
	  {
	    HitCentralWires[view]  += (*itHit)->WireID().Wire;
	    HitCentralTimes[view]  += (*itHit)->PeakTime();
	  }
	
	HitCentralTimes[view]/=HitsThisPlane.size();
	HitCentralWires[view]/=HitsThisPlane.size();
	NHits[view]=HitsThisPlane.size();
	
	AverageXWeighted      += det->ConvertTicksToX(HitCentralTimes[view],planeno,0,0) * NHits[view];
	AverageXDemocratic    += det->ConvertTicksToX(HitCentralTimes[view],planeno,0,0);
	TotalHits             += NHits[view];

	HitsThisPlane.clear();
      }
    
    // Times are easy - no coordinate degeneracy
    AverageXDemocratic /= Planes;
    AverageXWeighted   /= TotalHits;
    

    // Begin iterative refit procedure
    // Loop the prescribed number of times
    for(int loop=0; loop!=fRefits; loop++)
      {
	for(size_t view=0; view!=Planes; view++)
          {
	    // Adjust central seed point
	    // --------------------------

            // Get seed central point in this view



            double SeedPt[3], Err[3], SeedDir[3];
            TheSeed.GetPoint(    SeedPt,   Err);
            TheSeed.GetDirection(SeedDir,  Err);

            TVector3 SeedPoint(SeedPt[0],      SeedPt[1],  SeedPt[2]);
            TVector3 SeedDirection(SeedDir[0], SeedDir[1], SeedDir[2]);


	    // Find what the seed central wire was before
            double SeedCentralWire = (PlaneNormDirections.at(view).Dot(SeedPoint) - WireOffsets.at(view)) / WirePitch;
            

	    // Move it half way to that predicted from hits
	  
	    double InPlaneShift = 0.5*(HitCentralWires[view] - SeedCentralWire) * WirePitch;

	    SeedCentralWire = 0.5*(HitCentralWires[view] + SeedCentralWire);          



            SeedPt[0] = AverageXDemocratic;
            SeedPt[1] = SeedPt[1] + InPlaneShift * PlaneNormDirections.at(view)[1] ;
            SeedPt[2] = SeedPt[2] + InPlaneShift * PlaneNormDirections.at(view)[2] ;

	    // Update the seed
            TheSeed.SetPoint(SeedPt);


            for(int i=0; i!=3; i++)
              SeedPoint[i]=SeedPt[i];



            // Adjust direction
            //-----------------------

            double d2=0, d=0, N=0, dx=0, x=0;

	    art::PtrVector<recob::Hit> HitsThisPlane = HitMap[view]; 
	    unsigned int c=0, t=0, p=0;
            for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisPlane.begin(); itHit!=HitsThisPlane.end(); itHit++)
              {
		double DDisp = (double((*itHit)->WireID().Wire)-SeedCentralWire)*WirePitch; 
	       double Time  = (*itHit)->PeakTime();
	       double XDisp = det->ConvertTicksToX(Time,p,t,c) - AverageXDemocratic;

	       x  += XDisp;
	       d  += DDisp;
	       dx += DDisp*XDisp;
	       N  += 1;
	       d2 += DDisp*DDisp;

              }
	    TMatrixD LeastSquaresLHS(2,2);
            LeastSquaresLHS[0][0]=d2;
	    LeastSquaresLHS[0][1]=d;
            LeastSquaresLHS[1][0]=d;
            LeastSquaresLHS[1][1]=N;

	    TVectorD LeastSquaresRHS(2);
	    LeastSquaresRHS[0]=dx;
            LeastSquaresRHS[1]=x;
	    try{
	      LeastSquaresLHS.Invert();
	    }
	    catch(...)
	      {
		return;
	      }
            TVectorD LeastSquaresResult(2);
	    LeastSquaresResult = LeastSquaresLHS*LeastSquaresRHS;
	    // parameters for t = ac + b
            double a = LeastSquaresResult[0];
            //   double b = LeastSquaresResult[1];


	    // Project seed direction vector into parts parallel and perp to pitch

            double SeedPlaneComp     = SeedDirection.Dot(PlaneNormDirections.at(view).Unit());
            double SeedTimeComp      = SeedDirection.Dot(XVec.Unit());
	    double SeedOutOfPlaneComp= SeedDirection.Dot(WireVecs[view].Unit());


            double Newa=0.5*(a+SeedTimeComp/SeedPlaneComp);

            double OriginalSeedLength = SeedDirection.Mag();

            // Set seed direction to this theta without changing length
            SeedPlaneComp  = SeedTimeComp / Newa;

            // build the new direction from these 3 orthogonal components
            SeedDirection =
              SeedTimeComp       * XVec +
              SeedPlaneComp      * PlaneNormDirections.at(view) +
              SeedOutOfPlaneComp * WireVecs.at(view);

            SeedDirection.SetMag(OriginalSeedLength);

            for(int i=0; i!=3; i++)
              SeedDir[i]=SeedDirection[i];

            TheSeed.SetDirection(SeedDir);
	    
	    HitsThisPlane.clear();
          } // next plane
      } // next iteration

  
  }





  //------------------------------------------------------------

  std::vector<recob::SpacePoint> SeedFinderAlgorithm::ExtractSpacePoints(std::vector<recob::SpacePoint> const& AllPoints, std::vector<int> IDsToExtract)
  {
    std::vector<recob::SpacePoint> ReturnVec;
    for(size_t i=0; i!=IDsToExtract.size(); ++i)
      ReturnVec.push_back(AllPoints.at(IDsToExtract.at(i)));

    return ReturnVec;
  }

  //-----------------------------------------------------------

  void  SeedFinderAlgorithm::GetCenterAndDirection(art::PtrVector<recob::Hit> const& HitsFlat, std::vector<int>&  HitsToUse, TVector3& Center, TVector3& Direction)
  {
    // Initialize the services we need
    art::ServiceHandle<util::DetectorProperties>     det;

  
 
    std::map<uint32_t, bool>   HitsClaimed;

    // We'll store hit coordinates in each view into this vector

    std::vector<std::vector<double> > HitTimes(3);
    std::vector<std::vector<double> > HitWires(3);
    std::vector<std::vector<double> > HitWidths(3);
    std::vector<double> MeanWireCoord(3,0);
    std::vector<double> MeanTimeCoord(3,0);

    // Run through the collection getting hit info for these spacepoints

    std::vector<double> x(3,0), y(3,0), xx(3,0), xy(3,0), yy(3,0),  sig(3,0), N(3,0);
    
    for(size_t i=0; i!=HitsToUse.size(); ++i)
      {
	auto itHit = HitsFlat.begin()+HitsToUse[i];
	
	size_t ViewIndex;
	
	if(     (*itHit)->View() == geo::kU) ViewIndex=0;
	else if((*itHit)->View() == geo::kV) ViewIndex=1;
	else if((*itHit)->View() == geo::kW) ViewIndex=2;
	double WireCoord = (*itHit)->WireID().Wire * fPitches.at(ViewIndex);
	double TimeCoord = det->ConvertTicksToX((*itHit)->PeakTime(),ViewIndex,0,0);
	double TimeUpper = det->ConvertTicksToX((*itHit)->EndTime(), ViewIndex,0,0);
	double TimeLower = det->ConvertTicksToX((*itHit)->StartTime(), ViewIndex,0,0);
	double Width = fabs(0.5*(TimeUpper-TimeLower));
	double Width2 = pow(Width,2);
	
	HitWires.at(ViewIndex).push_back(WireCoord);
	HitTimes.at(ViewIndex).push_back(TimeCoord);
	HitWidths.at(ViewIndex).push_back(fabs(0.5*(TimeUpper-TimeLower)));
	
	MeanWireCoord.at(ViewIndex) += WireCoord;
	MeanTimeCoord.at(ViewIndex) += TimeCoord;
	
	// Elements for LS
	x.at(ViewIndex)   += WireCoord / Width2;
	y.at(ViewIndex)   += TimeCoord / Width2;
	xy.at(ViewIndex)  += (TimeCoord * WireCoord) / Width2;
	xx.at(ViewIndex)  += (WireCoord * WireCoord) / Width2;
	yy.at(ViewIndex)  += (TimeCoord * TimeCoord) / Width2;
	sig.at(ViewIndex) += 1./Width2;
	N.at(ViewIndex)++;
      }

    std::vector<double>   ViewRMS(3);
    std::vector<double>   ViewGrad(3);
    std::vector<double>   ViewOffset(3);
    
    std::vector<double> CenterTime(3);
    std::vector<double> CenterWireAtT0(3);
    


    for(size_t n=0; n!=3; ++n)
      {
	MeanWireCoord[n] /= N[n];
	MeanTimeCoord[n] /= N[n];
	
	ViewGrad.at(n)   = (y[n]/sig[n] - xy[n]/x[n])/(x[n]/sig[n]-xx[n]/x[n]);
	ViewOffset.at(n) = (y[n] - ViewGrad[n]*x[n])/sig[n];
	ViewRMS.at(n)    = (yy[n] + pow(ViewGrad[n],2) * xx[n] + pow(ViewOffset[n],2) - 2*ViewGrad[n]*xy[n] - 2*ViewOffset[n]*y[n] + 2*ViewGrad[n]*ViewOffset[n]*x[n]);

      }


    std::vector<TVector3> ViewDir2D(3);

    std::vector<TVector3> Center2D(3);

    // Calculate centers and directions from pairs
    for(size_t n=0; n!=3; ++n)
      {
	int n1 = (n+1)%3;
	int n2 = (n+2)%3;


	ViewDir2D[n] =
	  (fXDir + fPitchDir[n1] *(1./ViewGrad[n1])
	   + fWireDir[n1]  * (((1./ViewGrad[n2]) - fPitchDir[n1].Dot(fPitchDir[n2]) * (1./ViewGrad[n1])) / fWireDir[n1].Dot(fPitchDir[n2])) ).Unit();
	
	Center2D[n] =
	  fXDir * 0.5 * (MeanTimeCoord[n1]+MeanTimeCoord[n2])
	  + fPitchDir[n1] * (MeanWireCoord[n1] + fWireZeroOffset[n1])
	  + fWireDir[n1] *  ( ((MeanWireCoord[n2] + fWireZeroOffset[n2]) - ( MeanWireCoord[n1] + fWireZeroOffset[n1] )*fPitchDir[n1].Dot(fPitchDir[n2]))/(fPitchDir[n2].Dot(fWireDir[n1])) );
      }
    
    
    
    for(size_t n=0; n!=3; ++n)
      {
	size_t n1 = (n+1)%3;
	size_t n2 = (n+2)%3;
	if( (N[n] <= N[n1]) &&
	    (N[n] <= N[n2]) )
	  {
	    
	    Center    = Center2D[n];
	    Direction = ViewDir2D[n];
	  
	    ftTheta = Direction.Theta();
	    TVector3 NX(1,0,0), NY(0,1,0);
	    ftThetaYZ = (Direction - Direction.Dot(NX)*NX).Theta();
	    ftThetaXZ = (Direction - Direction.Dot(NY)*NY).Theta();
	    
	  }
      }
  }
  




  //------------------------------------------------------------
  // Given a set of spacepoints, find preferred direction
  //  and center

  void  SeedFinderAlgorithm::GetCenterAndDirection(std::vector<recob::SpacePoint> const& Points, std::vector<int> const& PointsInRange, TVector3& Center, TVector3& Direction)
  {

    // Initialize the services we need
    
    art::ServiceHandle<geo::Geometry>                geom;
    art::ServiceHandle<util::DetectorProperties>     det;
    
    
    std::map<uint32_t, bool>   HitsClaimed;
    
    // We'll store hit coordinates in each view into this vector
    
    std::vector<std::vector<double> > HitTimes(3);
    std::vector<std::vector<double> > HitWires(3);
    std::vector<std::vector<double> > HitWidths(3);
    std::vector<double> MeanWireCoord(3,0);
    std::vector<double> MeanTimeCoord(3,0);
    
    // Run through the collection getting hit info for these spacepoints
    
    std::vector<double> x(3,0), y(3,0), xx(3,0), xy(3,0), yy(3,0),  sig(3,0), N(3,0);
    
    std::vector<art::PtrVector<recob::Hit> > HitsInThisCollection(3);
    for(unsigned int i=0; i!=PointsInRange.size(); i++)
      {
	art::PtrVector<recob::Hit> HitsThisSP = fSptalg->getAssociatedHits((Points.at(PointsInRange.at(i))));
	for(art::PtrVector<recob::Hit>::const_iterator itHit=HitsThisSP.begin();
	    itHit!=HitsThisSP.end(); ++itHit)
	  {
	    if(!HitsClaimed[(*itHit)->Channel()])
	      {
		HitsClaimed[(*itHit)->Channel()]=true;
		
		size_t ViewIndex;
		
		if(     (*itHit)->View() == geo::kU) ViewIndex=0;
		else if((*itHit)->View() == geo::kV) ViewIndex=1; 
		else if((*itHit)->View() == geo::kW) ViewIndex=2;
		double WireCoord = (*itHit)->WireID().Wire * fPitches.at(ViewIndex);
		double TimeCoord = det->ConvertTicksToX((*itHit)->PeakTime(),ViewIndex,0,0);
		double TimeUpper = det->ConvertTicksToX((*itHit)->EndTime(), ViewIndex,0,0);
		double TimeLower = det->ConvertTicksToX((*itHit)->StartTime(), ViewIndex,0,0);
		double Width = fabs(0.5*(TimeUpper-TimeLower));
		double Width2 = pow(Width,2);
		
		HitWires.at(ViewIndex).push_back(WireCoord);
		HitTimes.at(ViewIndex).push_back(TimeCoord);
		HitWidths.at(ViewIndex).push_back(fabs(0.5*(TimeUpper-TimeLower)));
		
		HitsInThisCollection.at(ViewIndex).push_back(*itHit);
		
		MeanWireCoord.at(ViewIndex) += WireCoord;
		MeanTimeCoord.at(ViewIndex) += TimeCoord;
		
		// Elements for PCA
		x.at(ViewIndex)   += WireCoord / Width2;
		y.at(ViewIndex)   += TimeCoord / Width2;
		xy.at(ViewIndex)  += (TimeCoord * WireCoord) / Width2;
		xx.at(ViewIndex)  += (WireCoord * WireCoord) / Width2;
		yy.at(ViewIndex)  += (TimeCoord * TimeCoord) / Width2;
		sig.at(ViewIndex) += 1./Width2;
		N.at(ViewIndex)++;
		
	      }
	  }
      }

    std::vector<double>   ViewRMS(3);
    std::vector<double>   ViewGrad(3);
    std::vector<double>   ViewOffset(3);
    
    std::vector<double> CenterTime(3);
    std::vector<double> CenterWireAtT0(3);
    
    
    
    for(size_t n=0; n!=3; ++n)
      {
	MeanWireCoord[n] /= N[n];
	MeanTimeCoord[n] /= N[n];
	
	ViewGrad.at(n)   = (y[n]/sig[n] - xy[n]/x[n])/(x[n]/sig[n]-xx[n]/x[n]);
	ViewOffset.at(n) = (y[n] - ViewGrad[n]*x[n])/sig[n];
	ViewRMS.at(n)    = (yy[n] + pow(ViewGrad[n],2) * xx[n] + pow(ViewOffset[n],2) - 2*ViewGrad[n]*xy[n] - 2*ViewOffset[n]*y[n] + 2*ViewGrad[n]*ViewOffset[n]*x[n]);	    
	
      }
    
    
    std::vector<TVector3> ViewDir2D(3);
    
    std::vector<TVector3> Center2D(3);

    // Calculate centers and directions from pairs
    for(size_t n=0; n!=3; ++n)
      {
	int n1 = (n+1)%3;
	int n2 = (n+2)%3;
	
	
	ViewDir2D[n] =
	  (fXDir + fPitchDir[n1] *(1./ViewGrad[n1])
	   + fWireDir[n1]  * (((1./ViewGrad[n2]) - fPitchDir[n1].Dot(fPitchDir[n2]) * (1./ViewGrad[n1])) / fWireDir[n1].Dot(fPitchDir[n2])) ).Unit();
	
	Center2D[n] =
	  fXDir * 0.5 * (MeanTimeCoord[n1]+MeanTimeCoord[n2])
	  + fPitchDir[n1] * (MeanWireCoord[n1] + fWireZeroOffset[n1])
	  + fWireDir[n1] *  ( ((MeanWireCoord[n2] + fWireZeroOffset[n2]) - ( MeanWireCoord[n1] + fWireZeroOffset[n1] )*fPitchDir[n1].Dot(fPitchDir[n2]))/(fPitchDir[n2].Dot(fWireDir[n1])) );
      }
    
    
    
    for(size_t n=0; n!=3; ++n)
      {
	size_t n1 = (n+1)%3;
	size_t n2 = (n+2)%3;
	if( (N[n] <= N[n1]) &&
	    (N[n] <= N[n2]) )
	  {
	    
	    Center    = Center2D[n];
	    Direction = ViewDir2D[n];
	  	    
	    ftTheta = Direction.Theta();
	    TVector3 NX(1,0,0), NY(0,1,0);
	    ftThetaYZ = (Direction - Direction.Dot(NX)*NX).Theta();
	    ftThetaXZ = (Direction - Direction.Dot(NY)*NY).Theta();
	    
	  }
      }
    
    
  }



  //-----------------------------------------------
  void SeedFinderAlgorithm::CalculateGeometricalElements()
  {
    art::ServiceHandle<geo::Geometry> geom;
    
    // Find pitch of each wireplane
    fPitches.resize(3);
    fPitches.at(0) = fabs(geom->WirePitch(geo::kU));
    fPitches.at(1) = fabs(geom->WirePitch(geo::kV));
    fPitches.at(2) = fabs(geom->WirePitch(geo::kW));

    // Setup basis vectors
    fXDir = TVector3(1,0,0);
    fYDir = TVector3(0,1,0);
    fZDir = TVector3(0,0,1);

    fWireDir.resize(3);
    fPitchDir.resize(3);
    fWireZeroOffset.resize(3);

    double xyzStart1[3], xyzStart2[3];
    double xyzEnd1[3], xyzEnd2[3];

    // Calculate wire coordinate systems
    for(size_t n=0; n!=3; ++n)
      {
	geom->WireEndPoints(0,0,n,0,xyzStart1,xyzEnd1);
	geom->WireEndPoints(0,0,n,1,xyzStart2,xyzEnd2);
	fWireDir[n] = TVector3(xyzEnd1[0] - xyzStart1[0],
			      xyzEnd1[1] - xyzStart1[1],
			      xyzEnd1[2] - xyzStart1[2]).Unit();
	fPitchDir[n] = fWireDir[n].Cross(fXDir).Unit();
	if(fPitchDir[n].Dot(TVector3(xyzEnd2[0] - xyzEnd1[0],
				     xyzEnd2[1] - xyzEnd1[1],
				     xyzEnd2[2] - xyzEnd1[2]))<0) fPitchDir[n] = -fPitchDir[n];
	
	fWireZeroOffset[n] =
	  xyzEnd1[0]*fPitchDir[n][0] +
	  xyzEnd1[1]*fPitchDir[n][1] +
	  xyzEnd1[2]*fPitchDir[n][2];
	
      }


  }


  //-----------------------------------------------

  std::vector<recob::Seed>    SeedFinderAlgorithm::GetSeedsFromUnSortedHits(art::PtrVector<recob::Hit> const & Hits, std::vector<art::PtrVector<recob::Hit> >& HitCatalogue)
  {  
    std::vector<recob::Seed> ReturnVec;
  
    ReturnVec = FindSeeds( Hits, HitCatalogue);          
    
    return ReturnVec;
  }


 

  //---------------------------------------------

  std::vector<std::vector<recob::Seed> > SeedFinderAlgorithm::GetSeedsFromSortedHits(std::map<geo::View_t, std::vector<art::PtrVector<recob::Hit> > >  const& SortedHits, std::vector<std::vector<art::PtrVector<recob::Hit> > >& HitsPerSeed)
  {


    std::vector<std::vector<recob::Seed> > ReturnVec;

    // This piece of code looks detector specific, but its not -
    //   it also works for 2 planes, but one vector is empty.

    if(!(fSptalg->enableU()&&fSptalg->enableV()&&fSptalg->enableW()))
	mf::LogWarning("BezierTrackerModule")<<"Warning: SpacePointAlg is does not have three views enabled. This may cause unexpected behaviour in the bezier tracker.";
    
    
      try
        {
          for(std::vector<art::PtrVector<recob::Hit> >::const_iterator itU = SortedHits.at(geo::kU).begin();
              itU !=SortedHits.at(geo::kU).end(); ++itU)
            for(std::vector<art::PtrVector<recob::Hit> >::const_iterator itV = SortedHits.at(geo::kV).begin();
                itV !=SortedHits.at(geo::kV).end(); ++itV)
	      for(std::vector<art::PtrVector<recob::Hit> >::const_iterator itW = SortedHits.at(geo::kW).begin();
                  itW !=SortedHits.at(geo::kW).end(); ++itW)
                {
		  art::PtrVector<recob::Hit>    HitsThisComboFlat;

                  if(fSptalg->enableU())
                    for(size_t i=0; i!=itU->size(); ++i)
		      HitsThisComboFlat.push_back(itU->at(i));
		    
                  if(fSptalg->enableV())
                    for(size_t i=0; i!=itV->size(); ++i)
		      HitsThisComboFlat.push_back(itV->at(i));
		    
                  if(fSptalg->enableW())
                    for(size_t i=0; i!=itW->size(); ++i)
		      HitsThisComboFlat.push_back(itW->at(i));
		    
		  std::vector<art::PtrVector<recob::Hit> > CataloguedHits;
		      
		  
		  std::vector<recob::Seed> Seeds
		    = FindSeeds(HitsThisComboFlat, CataloguedHits);
		      
		  // Add this harvest to return vectors
		  HitsPerSeed.push_back(CataloguedHits);		      
		  ReturnVec.push_back(Seeds);

		  // Tidy up
		  CataloguedHits.clear();
                  Seeds.clear();
		}
	}
      catch(...)
        {
	  mf::LogWarning("SeedFinderTrackerModule")<<" bailed during hit map lookup - have you enabled all 3 planes?";
          ReturnVec.push_back(std::vector<recob::Seed>());
        }

      return ReturnVec;

  }
  

}
