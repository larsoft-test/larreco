//
// Name: BezierTrackerAlgorithm.cxx
//
// Purpose: Implementation file for module BezierTrackerAlgorithm.
//
// Ben Jones, MIT
//

#include <vector>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "RecoAlg/BezierTrackerAlgorithm.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "RecoObjects/BezierTrack.h"
#include "RecoObjects/BezierCurveHelper.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"


#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 

#include "TVector3.h"

namespace trkf {

  bool BTrack_SeedCountComparator(std::vector<recob::Seed> const& s1, std::vector<recob::Seed> const& s2);


  BezierTrackerAlgorithm::BezierTrackerAlgorithm(const fhicl::ParameterSet& pset)
  {
    
    reconfigure(pset);
  }

  //------------------------------------------------

  BezierTrackerAlgorithm::~BezierTrackerAlgorithm()
  {
  }

  //------------------------------------------------


  void BezierTrackerAlgorithm::reconfigure(fhicl::ParameterSet const& pset)
  {
    
    CalculateGeometricalElements();
    
    fHitDistance       = pset.get<double>("HitDistance");
    fTrackJoinAngle    = pset.get<double>("TrackJoinAngle");
    fDirectJoinDistance= pset.get<double>("DirectJoinDistance");
    
    fVertexImpactThreshold = pset.get<double>("VertexImpactThreshold");
    fVertexExtrapDistance = pset.get<double>("VertexExtrapDistance");
    
    fOccupancyThresh.resize(3);
    fOccupancyThresh = pset.get<std::vector<double> >("OccupancyThresh");
    fTrackResolution = pset.get<double>("TrackResolution");
    fOverlapCut      = pset.get<double>("OverlapCut");



    fTheSeedFinder = new SeedFinderAlgorithm(pset.get<fhicl::ParameterSet>("SeedFinder"));
    
  }


  //-----------------------------------------------
  std::vector<trkf::BezierTrack> BezierTrackerAlgorithm::MakeTracksNew(std::map<geo::View_t, std::vector<art::PtrVector<recob::Hit> > >& SortedHits, std::vector<art::PtrVector<recob::Hit> >& HitAssocs)
  {
    
    std::vector<trkf::BezierTrack> ReturnVector;

    size_t UEntries = SortedHits[geo::kU].size();
    size_t VEntries = SortedHits[geo::kV].size();
    size_t WEntries = SortedHits[geo::kW].size();

    
    // Get seeds from the hit collection    
    std::vector<std::map<size_t, bool> > MapNeeded(3);
 

    std::vector<std::vector<art::PtrVector<recob::Hit> > > HitsPerSeed;
    std::vector<std::vector<recob::Seed> > Seeds = GetSeedFinderAlgorithm()->GetSeedsFromSortedHits(SortedHits, HitsPerSeed);


    // Find the longest seed
    for(size_t nU =0; nU != UEntries; ++nU)
      for(size_t nV =0; nV != VEntries; ++nV)
	for(size_t nW =0; nW != WEntries; ++nW)
	  
	  {
	    size_t iCombo = 
	      nU * VEntries * WEntries + 
	      nV * WEntries +
	      nW;
	    
	    
	    if(Seeds.at(iCombo).size()>0)
	      {
		MapNeeded[0][nU]=true;
		MapNeeded[1][nV]=true;
		MapNeeded[2][nW]=true;
	      }
	  }
    

    std::vector< std::map<uint32_t, std::vector<int> > > OrgHitsU(UEntries);
    std::vector< std::map<uint32_t, std::vector<int> > > OrgHitsV(VEntries);
    std::vector< std::map<uint32_t, std::vector<int> > > OrgHitsW(WEntries);

    for(size_t nU=0; nU!=UEntries; ++nU)
      {
	if(MapNeeded[0][nU])
	  for(size_t iH=0; iH!=SortedHits[geo::kU][nU].size(); ++iH)
	    {
	      OrgHitsU[nU][SortedHits[geo::kU][nU][iH]->Channel()].push_back(iH);
	    }
      }
    for(size_t nV=0; nV!=VEntries; ++nV)
      {
	if(MapNeeded[1][nV])
	  for(size_t iH=0; iH!=SortedHits[geo::kV][nV].size(); ++iH)
	    {
	      OrgHitsV[nV][SortedHits[geo::kV][nV][iH]->Channel()].push_back(iH);
	    }
      }
    for(size_t nW=0; nW!=WEntries; ++nW)
      {
	if(MapNeeded[2][nW])
	  for(size_t iH=0; iH!=SortedHits[geo::kW][nW].size(); ++iH)
	    {
	      OrgHitsW[nW][SortedHits[geo::kW][nW][iH]->Channel()].push_back(iH);
	    }
      }
    
    for(size_t nU =0; nU != UEntries; ++nU)
      for(size_t nV =0; nV != VEntries; ++nV)
	for(size_t nW =0; nW != WEntries; ++nW)
	  
	  {
	    size_t iCombo = 
	      nU * VEntries * WEntries + 
	      nV * WEntries +
	      nW;
	    
	    
	    if(Seeds.at(iCombo).size()>0)
	      {
		
		std::vector<art::PtrVector<recob::Hit>* > HitStruct(3);
		HitStruct[0] = & SortedHits[geo::kU].at(nU);
		HitStruct[1] = & SortedHits[geo::kV].at(nV);
		HitStruct[2] = & SortedHits[geo::kW].at(nW);
		
		std::vector<std::map<uint32_t, std::vector<int> >* > OrgHits(3);
		OrgHits[0] = & OrgHitsU[nU];
		OrgHits[1] = & OrgHitsV[nV];
		OrgHits[2] = & OrgHitsW[nW];
		
		std::vector<std::vector<recob::Seed> > SeedsForTracks = OrganizeSeedsIntoTracksNew(Seeds.at(iCombo), HitStruct, HitsPerSeed.at(iCombo), OrgHits);
		//	mf::LogVerbatim("BezierTrackerAlgorithm")
		//	std::cout<<"This combo had " << Seeds.at(iCombo).size()<<" seeds and got " << SeedsForTracks.size()<<" tracks"<<std::endl;
		
		for(size_t i=0; i!=SeedsForTracks.size(); ++i)
		  {
		    //		    mf::LogVerbatim("BezierTrackerAlgorithm")<<" " 
		    //  std::cout<< i <<" with " <<SeedsForTracks.at(i).size()<<" seeds"<<std::endl;
		    ReturnVector.push_back(trkf::BezierTrack(SeedsForTracks.at(i)));
		    // This line temporary - replace with real assoc list
		    HitAssocs.push_back(art::PtrVector<recob::Hit>());
		  }
	      }
	    
	  } // End of loop through combos
    return ReturnVector;
    
    
  }




  //------------------------------------------------
  
  void BezierTrackerAlgorithm::MakeDirectJoins(std::vector<trkf::BezierTrack>& BTracks, std::vector<art::PtrVector<recob::Hit> > & HitVecs)
  {
    mf::LogInfo("BezierTrackerAlgorithm")<<"Making Direct Track Joins";
    
    std::vector<TVector3> End1Directions, End1Points, End2Directions, End2Points;
    for(size_t i=0; i!=BTracks.size(); ++i)
      {
	End1Points.push_back(BTracks.at(i).GetTrackPointV(0));
	End2Points.push_back(BTracks.at(i).GetTrackPointV(1));
	End1Directions.push_back(BTracks.at(i).GetTrackDirectionV(0));
	End2Directions.push_back(BTracks.at(i).GetTrackDirectionV(1));
      }
    

    std::map<int,bool> ToErase;

    for(size_t t1=0; t1!=BTracks.size(); ++t1)
      {
	if(!ToErase[t1]) 
	  for(size_t t2=0; t2!=BTracks.size(); ++t2)
	    {
	      if((!ToErase[t2])&&(t1!=t2))
		{
		  if( ( (fabs(End2Directions.at(t1).Angle(End1Directions.at(t2)))) < fTrackJoinAngle )
		      && ( (End2Points.at(t1)-End1Points.at(t2) ).Mag()<fDirectJoinDistance)
		      && ( (End2Directions.at(t1).Dot(End1Directions.at(t2)))>0)
		      && ( ( End2Points.at(t1)-End2Points.at(t2)).Mag()>fDirectJoinDistance)  
		      && ( ( End1Points.at(t1)-End1Points.at(t2)).Mag()>fDirectJoinDistance) )
		    
		    {
		      mf::LogVerbatim("BezierTrackerAlgorithm")<<" Making track join " << t1<<", " << t2;
		      
		      // Get combined seed collection
		      std::vector<recob::Seed> SeedCol1 = BTracks.at(t1).GetSeedVector();
		      std::vector<recob::Seed> SeedCol2 = BTracks.at(t2).GetSeedVector();
		      SeedCol1.pop_back();
		      SeedCol2.erase(SeedCol2.begin());
		      SeedCol1.insert(SeedCol1.end(), SeedCol2.begin(), SeedCol2.end());
		      // Replace track1 and remove track2
		      
		      ToErase[t2]=true;
		      mf::LogVerbatim("BezierTrackerAlgorithm")
			<< "Marking " << t2 << " for erase"<<std::endl;
		      
		      BTracks.at(t1) = trkf::BezierTrack(SeedCol1);
		      
		      // Add the pointervectors into 1 and remove 2
		      AddPtrVectors(HitVecs.at(t1), HitVecs.at(t2));
		      HitVecs.at(t2).clear();
				      
		      // Reset t1 and excape t2 loop to go around again
		      t1=-1;
		      break;
		    }
		}
	    }
      } 
    
    // Remove the tracks we marked for deletion
    mf::LogVerbatim("BezierTrackerAlgorithm")<<"Removing deletable elements";
    for(int i=BTracks.size()-1; i >= 0; --i)
      {
	if(ToErase[i])
	  {
	    BTracks.erase(BTracks.begin()+i);
	    HitVecs.erase(HitVecs.begin()+i);
	  }
      }
    
  }
  
  //-----------------------------------------------
  void BezierTrackerAlgorithm::AddPtrVectors(art::PtrVector<recob::Hit>& Receiver, art::PtrVector<recob::Hit> const & ToAdd)
  {
    for(size_t i=0; i!=ToAdd.size(); ++i)
      {
	Receiver.push_back(ToAdd.at(i));
      }
  } 
					     


  

  //----------------------------------------------------------------------
  //
  // From a PtrVector of hits, determine which are nearby
  //

  std::vector<int> BezierTrackerAlgorithm::DetermineNearbyHits(art::PtrVector<recob::Hit> const& Hits, BezierTrack const& BTrack, std::vector<double>& SValues)
  {
    std::vector<int> ReturnVector;
    std::vector<double> s, distance;
    BTrack.GetClosestApproaches(Hits,SValues, distance);

    for(size_t i=0; i!=Hits.size(); ++i)
      {

	if((distance.at(i)<fHitDistance)&&(SValues.at(i)<1)&&(SValues.at(i)>0)) ReturnVector.push_back(i); 

      }
    return ReturnVector;
  }


  


  //-----------------------------------------------------

  std::vector<std::vector< recob::Seed > > BezierTrackerAlgorithm::OrganizeSeedsIntoTracksNew(std::vector<recob::Seed > & AllSeeds, std::vector<art::PtrVector<recob::Hit>*>& AllHits, std::vector<art::PtrVector<recob::Hit> >& WhichHitsPerSeed, std::vector<std::map<uint32_t, std::vector<int> >* >& OrgHits )
  {
    mf::LogVerbatim("BezierTrackerAlgorithm")<<"Organizing seed collection into track collections ";

    std::vector<std::vector<uint32_t> > MaxSeedChanPerView(3);
    std::vector<std::vector<uint32_t> > MinSeedChanPerView(3);
    
    // This vector will keep track of which hits we are still allowed
    //  to use.  Key : 
    //   0 - hit available
    //   1 - hit taken by an extrapolated bezier curve
    //   2 - hit taken by a seed on the track

    std::vector<std::vector<int> > HitStatus(3);
    for(size_t n=0; n!=3; ++n) HitStatus.at(n).resize(AllHits[n]->size());


    // This map keeps track of the indices of the hits (in the AllHits vector)
    //  which go with each seed, in each view
    //   SeedToHitMap[per_view][per_seed][1_to_n_with_this_seed] = index 

    std::vector<std::vector<std::vector<int> > >  SeedToHitMap(3);
    
    // Fill Seed <--> Hit lookups
    for(size_t iSeed=0; iSeed!=AllSeeds.size(); ++iSeed)
      {	
	// first make sure we have at least one hit in each view.
	
	bool CompleteSet = false;
	std::vector<bool> AtLeastOneHitThisView(3,false);
	for(size_t iHit=0; iHit!=WhichHitsPerSeed.at(iSeed).size(); ++iHit)
	  {
	    if(!AtLeastOneHitThisView[WhichHitsPerSeed.at(iSeed).at(iHit)->View()])
	      {
		AtLeastOneHitThisView[WhichHitsPerSeed.at(iSeed).at(iHit)->View()] = true;
		if(AtLeastOneHitThisView[0]&&AtLeastOneHitThisView[1]&&AtLeastOneHitThisView[2])
		  {
		    CompleteSet = true;
		    break;
		  }
	      }
	  }
	if(!CompleteSet)
	  {
	    // If we don't have at least one hit in each view from this 
	    //  seed, we're not going to be able to use it.
	    //  ditch it (and its corresponding hits) from the game.
	    AllSeeds.erase(AllSeeds.begin() + iSeed);
	    WhichHitsPerSeed.erase(WhichHitsPerSeed.begin() + iSeed);
	    --iSeed;
	    continue;
	  }
	else
	  {
	    // if we didn't throw out that seed, push back a new entry in 
	    //  the seed <-> hit map vector, and 	 	    
	    // Find the min and max channels for each seed
	    for(size_t n=0; n!=3; ++n)
	      {
		SeedToHitMap[n].push_back(std::vector<int>());
		MinSeedChanPerView[n].push_back(10000);
		MaxSeedChanPerView[n].push_back(0);
	      }
	    for(size_t iHit=0; iHit!=WhichHitsPerSeed.at(iSeed).size(); ++iHit)
	      {
		uint32_t Channel = WhichHitsPerSeed.at(iSeed).at(iHit)->Channel();
		geo::View_t View = WhichHitsPerSeed.at(iSeed).at(iHit)->View();
		
		int ViewID=0;
		if(View==geo::kU)      ViewID=0; 
		else if(View==geo::kV) ViewID=1; 
		else if(View==geo::kW) ViewID=2;	    
		
		double      eta     = 0.01;
		
		for(size_t iH=0; iH!=OrgHits[ViewID]->operator[](Channel).size();++iH)
		  {
		    double PeakTime1 = AllHits[ViewID]->at(OrgHits[ViewID]->operator[](Channel).at(iH))->PeakTime();
		    double PeakTime2 = WhichHitsPerSeed.at(iSeed).at(iHit)->PeakTime();
		    
		    if( fabs(PeakTime1-PeakTime2)<eta)
		      
		      SeedToHitMap[ViewID][iSeed].push_back(OrgHits[ViewID]->operator[](Channel).at(iH));
		  }
		
	    
		if(Channel>MaxSeedChanPerView[ViewID][iSeed]) MaxSeedChanPerView[ViewID][iSeed] = Channel;
		if(Channel<MinSeedChanPerView[ViewID][iSeed]) MinSeedChanPerView[ViewID][iSeed] = Channel;
		
	      }
	  }
      }
    
    // Now we've set up the book keeping, we get into the real
    //  work of track finding.

    // This will store our proposed track vectors
    std::vector<std::vector<recob::Seed > > OrganizedByTrack;
     
    TVector3 AverageDir;

    std::vector<TVector3> VSeedDirs;
    std::vector<TVector3> VSeedPoss;
    std::vector<std::vector<int> >       SeedSDirs;
    std::vector<std::vector<uint32_t> >  SeedHighSChans;
    std::vector<std::vector<uint32_t> >  SeedLowSChans;
    
      

    for(size_t i=0; i!=AllSeeds.size(); ++i)
      {
	double SeedDir[3], SeedPos[3], Err[3];

	AllSeeds.at(i).GetDirection(SeedDir,Err);
	AllSeeds.at(i).GetPoint(SeedPos,Err);
	TVector3 VSeedDir = TVector3(SeedDir[0],SeedDir[1],SeedDir[2]);
	VSeedDirs.push_back(VSeedDir);
	TVector3 VSeedPos = TVector3(SeedPos[0],SeedPos[1],SeedPos[2]);
	VSeedPoss.push_back(VSeedPos);
	
	// This is a two directional degeneracy
	//  This is one way to resolve it (maybe not the best)

	if(VSeedDir.Dot(fZDir)<0)
	  AverageDir -= VSeedDir;
	else
	  AverageDir += VSeedDir;
	
	// Work out which direction the seed points out in each plane
	std::vector<double> TimeDir, WireDir;
	
	GetSeedDirProjected(AllSeeds.at(i), WireDir, TimeDir);
	
	SeedSDirs.push_back(std::vector<int>(3));
	SeedHighSChans.push_back(std::vector<uint32_t>(3));
	SeedLowSChans.push_back(std::vector<uint32_t>(3));
	

	for(size_t n=0; n!=3; ++n)
	  {
	    if(WireDir[n]>0)
	      {
		SeedSDirs[i][n]=1;
		SeedHighSChans[i][n] = MaxSeedChanPerView[n][i];
		SeedLowSChans[i][n]  = MinSeedChanPerView[n][i];
	      }
	    else
	      {
		SeedSDirs[i][n]=-1;
		SeedHighSChans[i][n] = MinSeedChanPerView[n][i];
		SeedLowSChans[i][n]  = MaxSeedChanPerView[n][i];
	      }
	  }
      }
    
    AverageDir *= (1.0/AverageDir.Mag());
       
    std::vector<double> AngleToAve;
    std::vector<double> ProjAlongAve;
    
    for(size_t i=0; i!=AllSeeds.size(); ++i)
      {
	ProjAlongAve.push_back(VSeedPoss.at(i).Dot(AverageDir));
	AngleToAve.push_back(VSeedDirs.at(i).Dot(AverageDir)/VSeedDirs.at(i).Mag());
      }

    
    // For now we'll only look for 1 long track per cluster combo.
    //  There is room for improvement here if we need it.
    
    std::vector<size_t> TrackParts;
    std::vector<size_t> DoesNotFit;
    std::vector<size_t> AlreadyUsed;
    std::vector<size_t> ThrownByHits;
    std::map<int, bool> SeedNotAvailable;
  
    // This loop keeps running whilst there are still tracks to make
    while((AlreadyUsed.size()+ThrownByHits.size())<AllSeeds.size())
      {
	       
	// This loop keeps running whilst there are still seeds
	//  to check for this track

	
	while((ThrownByHits.size()+TrackParts.size()+DoesNotFit.size()+AlreadyUsed.size())<AllSeeds.size())
	  {

	    // First find the lowest track proj remaining
	    int LowestLeft=-1;
	    for(size_t i=0; i!=AllSeeds.size(); ++i)
	      {
		if(!SeedNotAvailable[i])
		  {
		    if(LowestLeft < 0) LowestLeft = i;
		    if(ProjAlongAve.at(i)<ProjAlongAve.at(LowestLeft))
		      LowestLeft = i;		
		  }
	      }
	
	    std::vector<std::vector<int> > InterpHits(3);
	    
	    if(TrackParts.size()==0) 
	      {
		// If there are no seeds on the track, accept this seed
		//  as the first one to start running with.
		TrackParts.push_back(LowestLeft);
		SeedNotAvailable[LowestLeft]=true;
	
		// Mark the hits that go with this seed as taken
		for(size_t View=0; View!=3; ++View)
		  for(size_t iH=0; iH!=SeedToHitMap[View][LowestLeft].size();++iH)
		    HitStatus[View][SeedToHitMap[View][LowestLeft][iH]] = 2;
	      }
	    else if(TrackParts.size()>0) 
	      {
		int a = TrackParts.at(TrackParts.size()-1);
		int b = LowestLeft;

		double Angle           = AllSeeds.at(a).GetAngle(AllSeeds.at(b));
 		double ProjDisc        = AllSeeds.at(a).GetProjAngleDiscrepancy(AllSeeds.at(b));
		int PointingSign_ab    = AllSeeds.at(a).GetPointingSign(AllSeeds.at(b));
		int PointingSign_ba    = AllSeeds.at(a).GetPointingSign(AllSeeds.at(b));
		int PointingSign_ac    = 0;
		if(TrackParts.size()>1) PointingSign_ac = AllSeeds.at(a).GetPointingSign(AllSeeds.at(TrackParts.size()-2));
		
		bool SeedFits=false;

		// Check seed join conditions (no hits yet)
		bool JoinAngleOK    = fabs(Angle)<fTrackJoinAngle;
		bool PointingSignOK = (PointingSign_ac * PointingSign_ab)<0.5;
		bool ProjDiscOK     = false;
		if(PointingSign_ab > 0)
		  {
		    ProjDiscOK = fabs(ProjDisc) < fTrackJoinAngle;
		  }
		else
		  {
		    ProjDiscOK = fabs(ProjDisc - 3.142) < fTrackJoinAngle;
		  }
		
		if(JoinAngleOK && PointingSignOK && ProjDiscOK)
		  {
		    std::vector<uint32_t> EndChan_a(3), EndChan_b(3), LowChan(3), HighChan(3);
		    std::vector<std::vector<int> > HitsUnderConsideration(3);
  
		    for(size_t n=0; n!=3; ++n)
		      {
			if(PointingSign_ab > 0)
			  EndChan_a[n] = SeedHighSChans[a][n];   
			else
			  EndChan_a[n] = SeedLowSChans[a][n];
			
			if(PointingSign_ba > 0)
			  EndChan_b[n] = SeedHighSChans[b][n];   
			else
			  EndChan_b[n] = SeedLowSChans[b][n];
			
			LowChan[n]  = std::min(EndChan_a[n], EndChan_b[n]) ;
			HighChan[n] = std::max(EndChan_a[n], EndChan_b[n]) ;
		      }
		    
		    std::vector<std::vector<int> > TheseHits(3);
		    
		    std::vector<double> Occupancy = GetOccupancy(AllSeeds.at(a), AllSeeds.at(b), fTrackResolution, AllHits, OrgHits, LowChan, HighChan, HitStatus, TheseHits);
		    
		  
		    /*		    mf::LogVerbatim("BezierTrackerAlgorithm")<<"HighLow "<<a<<" : " << b <<
		      " " <<HighChan[0]-LowChan[0]<<
		      " " <<HighChan[1]-LowChan[1]<<
		      " " <<HighChan[2]-LowChan[2]<<
		      " " <<Occupancy[0]<<
		      " " <<Occupancy[1]<<
		      " " <<Occupancy[2]<<
		      std::endl;
		    */	    

		    if( (Occupancy[0] > fOccupancyThresh[0]) &&
			(Occupancy[1] > fOccupancyThresh[1]) &&
			(Occupancy[2] > fOccupancyThresh[2]) )		    
		      {		    
			SeedFits = true;
			InterpHits=TheseHits;
		      }

		  }
		if(SeedFits) 
		  {
		    // This seed is good! Store it on the track
		    TrackParts.push_back(LowestLeft);
		    SeedNotAvailable[LowestLeft]=true;

		    // Remove the hits that went with this seed from consideration
		    for(size_t View=0; View!=3; ++View)
		      for(size_t iH=0; iH!=SeedToHitMap[View][LowestLeft].size();++iH)
			HitStatus[View][SeedToHitMap[View][LowestLeft][iH]] = 2;
		    // Remove the hits that were interpolated through from consideration
		    for(size_t n=0; n!=3; ++n)
		      for(size_t i=0; i!=InterpHits[n].size(); ++i)
			{
			  HitStatus[n][InterpHits[n][i]] = 1;
			}
		    // Kick out any seeds which are not invalidated by their
		    //  hits being stolen
		    for(size_t j=0; j!=AllSeeds.size(); ++j)
		      {
			if(!SeedNotAvailable[j])
			  {
			    double OverlapFrac[3];
			    for(size_t n=0; n!=3; ++n)
			      {
				int TotalInSeed=0;
				int TotalRejected=0;
				for(size_t iH=0; iH!=SeedToHitMap[n][j].size(); ++iH)
				  {
				    TotalInSeed++;
				    if(HitStatus[n][SeedToHitMap[n][j][iH]]!=0)
				      TotalRejected++;
				  }
				OverlapFrac[n] = float(TotalRejected)/float(TotalInSeed);
			      }
			    for(size_t n=0; n!=3; ++n)
			      {
				size_t n1=(n+1)%3; size_t n2=(n+2)%3;
				if( (OverlapFrac[n]<=OverlapFrac[n1]) &&
				    (OverlapFrac[n]<=OverlapFrac[n2]) )
				  {
				    if(OverlapFrac[n] > fOverlapCut)
				      {
					SeedNotAvailable[j]=true;
					ThrownByHits.push_back(j);
					break;
				      }
				  }
			      }// End view loop (seed tossing)
			  } // End checking seed availability (seed tossing)
		      } // End seed loop (seed tossing)
		  } // End if positive seed fit
		else 
		  {
		    // If seed doesn't fit, mark it as such, and make it
		    //  temporarily unavailable.
		    DoesNotFit.push_back(LowestLeft);
		    SeedNotAvailable[LowestLeft]=true;
		  } 
	      } // End if(not first seed)
	  } // End loop over all remaining seeds
	if(TrackParts.size()>0)
	  {
	    mf::LogVerbatim("BezierTrackerAlgorithm")<<" Prototrack runs through seeds: " ;
	    std::vector<recob::Seed> TheTrackSeeds;
	    for(size_t i=0; i!=TrackParts.size(); ++i)
	      {
		mf::LogVerbatim("BezierTrackerAlgorithm")<<"  " << TrackParts.at(i)<<", ";
		TheTrackSeeds.push_back(AllSeeds.at(TrackParts.at(i)));
		AlreadyUsed.push_back(TrackParts.at(i));
	      }
	    OrganizedByTrack.push_back(TheTrackSeeds);  
	  }

	// Reset the count of unavailable seeds
	SeedNotAvailable.clear();

	// Anything which is either part of a track,
	//  or which has had its hits used, is not fair game.		
	for(size_t j=0; j!=AlreadyUsed.size(); ++j)
	  {
	    SeedNotAvailable[AlreadyUsed[j]]=true;
	  }

	// For all the seeds which didn't fit, see how many
	//  are now ruled out by hit overlaps
	for(size_t i=0; i!=DoesNotFit.size(); ++i)
	  {
	    size_t j=DoesNotFit[i];
	    double OverlapFrac[3];
	    for(size_t n=0; n!=3; ++n)
	      {
		int TotalInSeed=0;
		int TotalRejected=0;
		for(size_t iH=0; iH!=SeedToHitMap[n][j].size(); ++iH)
		  {
		    TotalInSeed++;
		    if(HitStatus[n][SeedToHitMap[n][j][iH]]!=0)
		      TotalRejected++;
		  }
		OverlapFrac[n] = float(TotalRejected)/float(TotalInSeed);
	      }
	    for(size_t n=0; n!=3; ++n)
	      {
		size_t n1=(n+1)%3; size_t n2=(n+2)%3;
		if( (OverlapFrac[n]<=OverlapFrac[n1]) &&
		    (OverlapFrac[n]<=OverlapFrac[n2]) )
		  {
		    if(OverlapFrac[n] > fOverlapCut)
		      {
			SeedNotAvailable[j]=true;
			ThrownByHits.push_back(j);
			break;
		      }
		  }
	      }// End view loop (seed tossing)
	  } // Loop over DoesNotFit
      
	
	for(size_t i=0; i!=ThrownByHits.size(); ++i)
	  SeedNotAvailable[ThrownByHits[i]]=true;
	
	// Ready to go around again
	TrackParts.clear();	
	DoesNotFit.clear();
      }
    
    // Sort vector so most seeded tracks appear first
    std::sort(OrganizedByTrack.begin(), OrganizedByTrack.end(), BTrack_SeedCountComparator);
    return OrganizedByTrack;
  }


  //-------------------------------------------------

  std::vector<double> BezierTrackerAlgorithm::GetOccupancy(recob::Seed& Seed1, recob::Seed& Seed2, double dThresh,  std::vector<art::PtrVector<recob::Hit>*>& AllHits,  std::vector<std::map<uint32_t, std::vector<int> >* >& OrgHits, std::vector<uint32_t>& LowChan, std::vector<uint32_t>& HighChan, std::vector<std::vector<int> >& HitStatus, std::vector<std::vector<int> >& TheseHits)
  {
    art::ServiceHandle<util::DetectorProperties> det;
    std::vector<double> Occupation(3);
    
    int NSteps = 2 * Seed1.GetDistance(Seed2) / dThresh;
    if(NSteps<10) NSteps=10;
    BezierCurveHelper bhlp(NSteps);
  
    std::vector<TVector3> Pts = bhlp.GetBezierPointsNew(Seed2, Seed1, NSteps);
  

    for(size_t n=0; n!=3; ++n)
      {
	
	int NClaimedChannels=0;
	std::vector<bool> ChannelsClaimed(HighChan[n]-LowChan[n]+1, false);
	if((HighChan[n]-LowChan[n])<2) Occupation[n]=1;
	else
	  {
	    for(uint32_t iChan = LowChan[n]+1; iChan < (HighChan[n]-1); ++iChan)
	      {
		for(size_t iH = 0; iH!=OrgHits[n]->operator[](iChan).size(); ++iH)
		  {
		    if(HitStatus[n][OrgHits[n]->operator[](iChan).at(iH)]==0)
		      {
			art::Ptr<recob::Hit> ThisHit = AllHits[n]->at(OrgHits[n]->operator[](iChan)[iH]);
		
			TVector3 WirePoint1 = ( fPitchDir[n]*(fWireZeroOffset[n] + fPitches[n] * ThisHit->WireID().Wire) )+ fXDir * det->ConvertTicksToX(ThisHit->PeakTime(), n, 0, 0);
			// TVector3 WirePoint2 = WirePoint1 + fWireDir[n];
			for(size_t ipt=0; ipt!=Pts.size(); ++ipt)
			  {
			    float d = ((WirePoint1-Pts.at(ipt))-fWireDir[n]*((WirePoint1-Pts.at(ipt)).Dot(fWireDir[n]))).Mag();
			    if(d < dThresh) 
			      {
				TheseHits[n].push_back(OrgHits[n]->operator[](iChan).at(iH));
				if(!ChannelsClaimed.at(iChan-LowChan[n]))
				  {
				    ChannelsClaimed.at(iChan - LowChan[n])=true;
				    ++NClaimedChannels;
				    ipt = Pts.size()-1;
				    continue;
				  }
			      }
			  }
		      }
		  }
	      }
	    Occupation[n] =  double(NClaimedChannels)/double(HighChan[n]-LowChan[n]-2);
	  }
      }
    return Occupation;
    
  }


  //-------------------------------------------------

  void BezierTrackerAlgorithm::GetSeedDirProjected(recob::Seed const& TheSeed, std::vector<double>& WireCoord, std::vector<double>& TimeCoord) 
  {
    art::ServiceHandle<util::DetectorProperties> det;
    
    WireCoord.clear();
    WireCoord.resize(3);
    TimeCoord.clear();
    TimeCoord.resize(3);
    
    double dir[3], err[3];
    
    TheSeed.GetDirection(dir,err);
    TVector3 DirVec(dir[0],dir[1],dir[2]);
    
    for(size_t n=0; n!=3; ++n)
      {
	WireCoord[n] = DirVec.Dot(fPitchDir[n]) / fPitches[n];
	TimeCoord[n] = det->ConvertXToTicks(DirVec.Dot(fXDir),n,0,0);
      }

  }




  //-------------------------------------------------


  void BezierTrackerAlgorithm::CalculateGeometricalElements()
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

  //-----------------------------------------------------

  void BezierTrackerAlgorithm::MakeVertexJoins(std::vector<trkf::BezierTrack>& BTracks, std::vector<recob::Vertex>& Vertices, std::vector<std::vector<int> > Mapping)
  {
    
    std::vector<TVector3> TrackEnd1s;
    std::vector<TVector3> TrackEnd2s;
    std::vector<TVector3> TrackDir1s;
    std::vector<TVector3> TrackDir2s;


    for(size_t i=0; i!=BTracks.size(); ++i)
      {
	TrackEnd1s.push_back( BTracks.at(i).GetTrackPointV(0));
	TrackEnd2s.push_back( BTracks.at(i).GetTrackPointV(1));
	TrackDir1s.push_back( BTracks.at(i).GetTrackDirectionV(0));
	TrackDir2s.push_back( BTracks.at(i).GetTrackDirectionV(1));	
      }
 
    std::vector<std::map<int,int> >     TracksMeetAtPoints;    
    std::vector<std::vector<TVector3> >  MeetPoints;
     
    // The connection[i][j] gives impact parameters for track ends
    //  The sign tells us which end of track j 
    for(size_t i=0; i!=BTracks.size(); ++i)
      {
	for(size_t j=0; j!=BTracks.size(); ++j)
	  {
	    if(i!=j)
	      {
		double impact, disti, distj;
		GetImpact(TrackEnd1s[i],TrackDir1s[i],TrackEnd1s[j],TrackDir1s[j], impact, disti, distj);

		if((impact < fVertexImpactThreshold)
		   && (fabs(disti) < fVertexExtrapDistance )
		   && (fabs(distj) < fVertexExtrapDistance )
		   && (disti < 0 )
		   && (distj < 0 ))
		  {
		    bool ThisTrackAlreadyMet =false;
		    TVector3 ExtrapPoint1 = TrackEnd1s[i] + TrackDir1s[j].Unit()*disti;
		    TVector3 ExtrapPoint2 = TrackEnd1s[j] + TrackDir1s[j].Unit()*distj;
		    
		    for(size_t pt=0; pt!=TracksMeetAtPoints.size(); ++pt)
		      {
			if(TracksMeetAtPoints.at(pt)[i]==-1)
			  {
			    TracksMeetAtPoints.at(pt)[j]=-1;
			    MeetPoints.at(pt).push_back(ExtrapPoint1);
			    MeetPoints.at(pt).push_back(ExtrapPoint2);
			    ThisTrackAlreadyMet=true;
			    
			  }
			else if(TracksMeetAtPoints.at(pt)[j]==-1)
			  {
			    TracksMeetAtPoints.at(pt)[i]=-1;
			    MeetPoints.at(pt).push_back(ExtrapPoint1);
			    MeetPoints.at(pt).push_back(ExtrapPoint2);
			    ThisTrackAlreadyMet=true;
			  }
		      }
		    if(!ThisTrackAlreadyMet)
		      {
			size_t pt = TracksMeetAtPoints.size();
			TracksMeetAtPoints.push_back(std::map<int,int>());
			TracksMeetAtPoints.at(pt)[i]=-1;
			TracksMeetAtPoints.at(pt)[j]=-1;
			MeetPoints.push_back(std::vector<TVector3>());
			MeetPoints.at(pt).push_back(ExtrapPoint1);
			MeetPoints.at(pt).push_back(ExtrapPoint2);
		      }

		  }

		GetImpact(TrackEnd2s[i],TrackDir2s[i],TrackEnd2s[j],TrackDir2s[j], impact, disti, distj);
		if((impact < fVertexImpactThreshold)
		   && (fabs(disti) < fVertexExtrapDistance )
		   && (fabs(distj) < fVertexExtrapDistance )
		   && (disti > 0 )
		   && (distj > 0 ))
		  {
		    bool ThisTrackAlreadyMet =false;
		    TVector3 ExtrapPoint1 = TrackEnd1s[i] + TrackDir2s[j].Unit()*disti;
		    TVector3 ExtrapPoint2 = TrackEnd2s[j] + TrackDir2s[j].Unit()*distj;
		    
		    for(size_t pt=0; pt!=TracksMeetAtPoints.size(); ++pt)
		      {
			if(TracksMeetAtPoints.at(pt)[i]==1)
			  {
			    TracksMeetAtPoints.at(pt)[j]=1;
			    MeetPoints.at(pt).push_back(ExtrapPoint1);
			    MeetPoints.at(pt).push_back(ExtrapPoint2);
			    ThisTrackAlreadyMet=true;
			    
			  }
			else if(TracksMeetAtPoints.at(pt)[j]==1)
			  {
			    TracksMeetAtPoints.at(pt)[i]=1;
			    MeetPoints.at(pt).push_back(ExtrapPoint1);
			    MeetPoints.at(pt).push_back(ExtrapPoint2);
			    ThisTrackAlreadyMet=true;
			  }
		      }
		    if(!ThisTrackAlreadyMet)
		      {
			size_t pt = TracksMeetAtPoints.size();
			TracksMeetAtPoints.push_back(std::map<int,int>());
			TracksMeetAtPoints.at(pt)[i]=1;
			TracksMeetAtPoints.at(pt)[j]=1;
			MeetPoints.push_back(std::vector<TVector3>());
			MeetPoints.at(pt).push_back(ExtrapPoint1);
			MeetPoints.at(pt).push_back(ExtrapPoint2);
		      }
		    
		    
		  }
		
		GetImpact(TrackEnd1s[i],TrackDir1s[i],TrackEnd2s[j],TrackDir2s[j], impact, disti, distj);
		if((impact < fVertexImpactThreshold)
		   && (fabs(disti) < fVertexExtrapDistance )
		   && (fabs(distj) < fVertexExtrapDistance )
		   && (disti < 0 )
		   && (distj > 0 ))
		  {
		    bool ThisTrackAlreadyMet =false;
		    TVector3 ExtrapPoint1 = TrackEnd1s[i] + TrackDir2s[j].Unit()*disti;
		    TVector3 ExtrapPoint2 = TrackEnd2s[j] + TrackDir2s[j].Unit()*distj;
		    
		    for(size_t pt=0; pt!=TracksMeetAtPoints.size(); ++pt)
		      {
			if(TracksMeetAtPoints.at(pt)[i]==-1)
			  {
			    TracksMeetAtPoints.at(pt)[j]=1;
			    MeetPoints.at(pt).push_back(ExtrapPoint1);
			    MeetPoints.at(pt).push_back(ExtrapPoint2);
			    ThisTrackAlreadyMet=true;
			    
			  }
			else if(TracksMeetAtPoints.at(pt)[j]==1)
			  {
			    TracksMeetAtPoints.at(pt)[i]=-1;
			    MeetPoints.at(pt).push_back(ExtrapPoint1);
			    MeetPoints.at(pt).push_back(ExtrapPoint2);
			    ThisTrackAlreadyMet=true;
			  }
		      }
		    if(!ThisTrackAlreadyMet)
		      {
			size_t pt = TracksMeetAtPoints.size();
			TracksMeetAtPoints.push_back(std::map<int,int>());
			TracksMeetAtPoints.at(pt)[i]=-1;
			TracksMeetAtPoints.at(pt)[j]=1;
			MeetPoints.push_back(std::vector<TVector3>());
			MeetPoints.at(pt).push_back(ExtrapPoint1);
			MeetPoints.at(pt).push_back(ExtrapPoint2);
		      }
	  
				    
		  }
		
	      }
	  }
      }


    for(size_t pt=0; pt!=TracksMeetAtPoints.size(); ++pt)
      {
	std::vector<int> TracksThisVertex;
	mf::LogVerbatim("BezierTrackerAlgorithm")<<" Making track adjustments for vertex " <<pt<<std::endl;
	TVector3 FinalMeetPt(0,0,0);
	for(size_t i=0; i!=MeetPoints.at(pt).size(); ++i)
	  {
	    FinalMeetPt += MeetPoints.at(pt).at(i);
	  } 
	FinalMeetPt *= (1./float(MeetPoints.at(pt).size()));
	
	for(size_t i=0; i!=BTracks.size(); ++i)
	  {
	    if(TracksMeetAtPoints.at(pt)[i]==1)
	      {
		TracksThisVertex.push_back(i);
		TVector3 NewDirection = (FinalMeetPt - TrackEnd2s[i])*0.5;
		double NewPos[3];
		double NewDir[3];
		double Err[3];
		
		for(size_t n=0; n!=3; ++n)
		  {
		    NewPos[n]=FinalMeetPt[n]-NewDirection[n];
		    NewDir[n]=NewDirection[n];
		  }
		
		recob::Seed NewSeed(NewPos, NewDir, Err, Err);
		std::vector<recob::Seed> SeedCol = BTracks.at(i).GetSeedVector();
		SeedCol.at(SeedCol.size()-1) = NewSeed;
		BTracks.at(i) = trkf::BezierTrack(SeedCol);
	      }
	    else if(TracksMeetAtPoints.at(pt)[i]==-1)
	      {
		TracksThisVertex.push_back(i);
		TVector3 NewDirection = -(FinalMeetPt - TrackEnd1s[i])*0.5;
		double NewPos[3];
		double NewDir[3];
		double Err[3];
		
		for(size_t n=0; n!=3; ++n)
		  {
		    NewPos[n]=FinalMeetPt[n] + NewDirection[n];
		    NewDir[n]=NewDirection[n];
		  }
		
		recob::Seed NewSeed(NewPos, NewDir, Err, Err);
		std::vector<recob::Seed> SeedCol = BTracks.at(i).GetSeedVector();
		SeedCol.at(0) = NewSeed;
		BTracks.at(i) = trkf::BezierTrack(SeedCol);
	      }

	    
	  }
	Mapping.push_back(TracksThisVertex);
	double VtxPoint[3];
	for(size_t n=0; n!=3; ++n)
	  VtxPoint[n]=FinalMeetPt[n];
	Vertices.push_back(recob::Vertex(VtxPoint, pt));
      }
    
  }


  //------------------------------------------------

  void BezierTrackerAlgorithm::GetImpact(TVector3 t1pt, TVector3 t1dir, TVector3 t2pt, TVector3 t2dir, double& ImpactParam, double& Dist1, double& Dist2)
  {
    double lambda1 = ((t2pt-t1pt).Cross(t2dir)).Dot(t1dir.Cross(t2dir)) / (t1dir.Cross(t2dir)).Mag2();
    double lambda2 = ((t1pt-t2pt).Cross(t1dir)).Dot(t2dir.Cross(t1dir)) / (t2dir.Cross(t1dir)).Mag2();
    TVector3 c = t1pt + t1dir*lambda1 - t2pt - t2dir*lambda2;
    
    ImpactParam = c.Mag();
    Dist1 = lambda1 * t1dir.Mag();
    Dist2 = lambda2 * t2dir.Mag();
  }

  
  //-------------------------------------------
  bool BTrack_SeedCountComparator(std::vector<recob::Seed> const& s1, std::vector<recob::Seed> const& s2) 
  {
    return s1.size()>s2.size();    
  }
}
