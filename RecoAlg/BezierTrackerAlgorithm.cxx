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
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "RecoObjects/BezierTrack.h"
#include "Utilities/AssociationUtil.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 


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

    fMaxKinkDThetaDx   = pset.get<double>("MaxKinkDThetaDx");
    fMaxJumpLengths    = pset.get<double>("MaxJumpLengths");
    fHitDistance       = pset.get<double>("HitDistance");
    fTrackJoinAngle    = pset.get<double>("TrackJoinAngle");
    fTheSeedFinder = new SeedFinderAlgorithm(pset.get<fhicl::ParameterSet>("SeedFinder"));
    
  }


  //------------------------------------------------

  void BezierTrackerAlgorithm::MakeBezierTracksFromSeeds(std::vector<trkf::BezierTrack>& ReturnVector, std::vector<recob::Seed> const& AllSeeds )
  {
    
    mf::LogVerbatim("BezierTrackerAlgorithm")<<"Making bezier tracks from seeds";
    
    std::vector<std::vector<recob::Seed> > OrgSeeds = OrganizeSeedsIntoTracks(AllSeeds);

    mf::LogVerbatim("BezierTrackerAlgorithm")<<"Producing track objects";    

    for(unsigned int i=0; i!=OrgSeeds.size(); i++)
      {
	if(OrgSeeds.at(i).size()>1)
	  ReturnVector.push_back(trkf::BezierTrack(OrgSeeds.at(i)));
	else if(OrgSeeds.at(i).size()==1)
	  {
	    // Throw out any stray seeds which lie on top of another track
	    bool DoNotTrack=false;
	    for(size_t j=0; j!=ReturnVector.size(); ++j)
	      {
		double s,d;
		ReturnVector.at(j).GetClosestApproach(OrgSeeds.at(i).at(0),s,d);
		//		mf::LogVerbatim("BezierTrackerAlgorithm")<<"Checking for overlap " << d<<std::endl;
		if(d<fHitDistance)
		  {
		    mf::LogVerbatim("BezierTrackerAlgorithm") << " Isolated seed track "<<i<<" is subsumed by track " << j << " in seed vector - doesn't deserve  a new track"<<std::endl;
		    DoNotTrack=true;
		  }
	      }
	    if(!DoNotTrack)
	      ReturnVector.push_back(trkf::BezierTrack(OrgSeeds.at(i)));
	  }
      }

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
		  if( ( (End2Directions.at(t2).Angle(End1Directions.at(t1))) < fTrackJoinAngle )
		      && ( (End2Points.at(t2)-End1Points.at(t1) ).Mag()<5.))
		    {
		      mf::LogVerbatim("BezierTrackerAlgorithm")<<" Making track join " << t1<<", " << t2<<std::endl;
		      
		      // Get combined seed collection
		      std::vector<recob::Seed> SeedCol1 = BTracks.at(t1).GetSeedVector();
		      std::vector<recob::Seed> SeedCol2 = BTracks.at(t2).GetSeedVector();
		      SeedCol1.insert(SeedCol1.end(), SeedCol2.begin(), SeedCol2.end());
		      
		      if(t1<t2)
			{
			  // Replace track1 and remove track2
			  
			  ToErase[t2]=true;
			  mf::LogVerbatim("BezierTrackerAlgorithm")
			    << "Marking " << t2 << " for erase"<<std::endl;
			  
			  BTracks.at(t1) = trkf::BezierTrack(SeedCol1);
			  
			  // Add the pointervectors into 1 and remove 2
			  AddPtrVectors(HitVecs.at(t1), HitVecs.at(t2));
			  HitVecs.at(t2).clear();
			}
		      else
			{
			  // Replace track2 and remove track1
			  
			  ToErase[t1]=true;
			  mf::LogVerbatim("BezierTrackerAlgorithm")
			    << "Marking " << t1 << " for erase"<<std::endl;
			  
			  BTracks.at(t2) = trkf::BezierTrack(SeedCol1);
			  
			  // Add the pointervectors into 1 and remove 2
			  AddPtrVectors(HitVecs.at(t2), HitVecs.at(t1));
			  HitVecs.at(t1).clear();
			  
			}
		      
		      
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
					     

  //------------------------------------------------

  void BezierTrackerAlgorithm::MakeBezierTracksFromHits(std::vector<trkf::BezierTrack>& ReturnVector, std::vector<art::Ptr<recob::Hit> > HitVec, std::vector<art::PtrVector<recob::Hit> >& HitsForAssns )
  {
    mf::LogVerbatim("BezierTrackerAlgorithm")<<"Making bezier tracks from hits"<<std::endl;
    
    
    // This vector keeps track of which hits we are still processing  
    art::PtrVector<recob::Hit> HitsToProcess;
    for(size_t i=0; i!=HitVec.size(); ++i)
      {
	HitsToProcess.push_back(HitVec.at(i));
      }

    bool KeepTrying=true;
    while(KeepTrying)
      {
	mf::LogVerbatim("BezierTrackerAlgorithm")<<"Getting space points" <<std::endl;
	
	// Make remaining hits into SPs
	std::vector<recob::SpacePoint> SPVec = 
	  fTheSeedFinder->GetSpacePointsFromHitVector(HitsToProcess);
	    
	// Find seeds in these SPs
	std::vector<std::vector<recob::SpacePoint> > SPUsed;
	mf::LogVerbatim("BezierTrackerAlgorithm")<<"Getting seeds " <<std::endl;
	
	std::vector<recob::Seed> AllSeeds = fTheSeedFinder->FindSeeds(SPVec, SPUsed);
	   	    
	mf::LogVerbatim("BezierTrackerAlgorithm")<<"Organizing seed collections " <<std::endl;
	// Organize these seeds into tracklike collections
	std::vector<std::vector<recob::Seed > > OrgSeeds = OrganizeSeedsIntoTracks(AllSeeds);
	
	// If we didn't find any, continue and quit
	if(OrgSeeds.size()==0) 
	  {
	    KeepTrying=false;
	    continue;
	  }
	    
	mf::LogVerbatim("BezierTrackerAlgorithm")<<"Making tracks from  " << OrgSeeds.size()<<" seed collections"<<std::endl;
	// For each of them make 1 track
	for(unsigned int i=0; i!=OrgSeeds.size(); i++)
	  {
	    
	    // Make a base track using seed positions and directions
	    //  and build a BezierTrack on top of it.		
	    BezierTrack BTrack(OrgSeeds.at(i));
	 	
	    art::PtrVector<recob::Hit>         HitsThisTrack;
	    std::vector<double>                SValues;
	    
	    std::vector<int> HitIDs = DetermineNearbyHits(HitsToProcess, BTrack, SValues);
	    
	    mf::LogVerbatim("BezierTrackerAlgorithm")<<"Found " << HitIDs.size()<<" nearby hits" <<std::endl;
	    
	    for(size_t i=0; i!=HitIDs.size(); i++)
	      HitsThisTrack.push_back(HitsToProcess.at(HitIDs.at(i)));
	    
	    HitsForAssns.push_back(HitsThisTrack);
	    
	    ReturnVector.push_back(BTrack);
	 
	    // Remove the hits we used from the vector and go around
	    mf::LogVerbatim("BezierTrackerAlgorithm")<<"Removing used hits"<<std::endl;
	    for(int i=HitIDs.size()-1; i!=-1; --i)
	      HitsToProcess.erase(HitsToProcess.begin() + HitIDs.at(i));
	  }
	mf::LogVerbatim("BezierTrackerAlgorithm")<<"Uncollected hit size : " << HitsToProcess.size()<<std::endl;
	
	if(HitsToProcess.size()<3) KeepTrying=false;
	
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


  

  //----------------------------------------------------------------------
  // 
  // Organize a vector of seeds into collections which seem to form sensible
  // track objects.
  //

  std::vector<std::vector< recob::Seed > > BezierTrackerAlgorithm::OrganizeSeedsIntoTracks(std::vector<recob::Seed > const& AllSeeds)
  {
    mf::LogVerbatim("BezierTrackerAlgorithm")<<"Organizing seed collection into track collections ";
    
    std::vector<std::vector<recob::Seed > > OrganizedByTrack;
    std::map<int, std::map<int, double> > ConnectionMap;
    
    TVector3 AverageDir;
    std::vector<TVector3> SeedDirs;
    std::vector<TVector3> SeedPoss;

    for(size_t i=0; i!=AllSeeds.size(); ++i)
      {
	double SeedDir[3], SeedPos[3], Err[3];

	AllSeeds.at(i).GetDirection(SeedDir,Err);
	AllSeeds.at(i).GetPoint(SeedPos,Err);
	SeedDirs.push_back(TVector3(SeedDir[0],SeedDir[1],SeedDir[2]));
	SeedPoss.push_back(TVector3(SeedPos[0],SeedPos[1],SeedPos[2]));
	AverageDir+=SeedDirs.at(SeedDirs.size()-1);
      }
    AverageDir *= (1.0/AverageDir.Mag());
       
	
    for(size_t i=0; i!=AllSeeds.size(); ++i)
      {
	for(size_t j=0; j!=AllSeeds.size(); ++j)
	  {
	    if(i==j) continue;
	    double Distance = AllSeeds.at(j).GetDistance(AllSeeds.at(i));
	    double Angle    = AllSeeds.at(j).GetAngle(AllSeeds.at(i));
	    double Length   = std::max(AllSeeds.at(i).GetLength(), AllSeeds.at(j).GetLength());
	    double ProjDisc = AllSeeds.at(i).GetProjAngleDiscrepancy(AllSeeds.at(j));
	    int PointingSign = AllSeeds.at(i).GetPointingSign(AllSeeds.at(j));
	    
	    
	    if( ( Angle<fTrackJoinAngle) && (ProjDisc<fTrackJoinAngle)&& (((Distance-Length)/Length)<fMaxJumpLengths)) 
	      {
		
		ConnectionMap[i][j] = Distance*float(PointingSign);
	       
	      }  
	  }
	
      }
    
    std::vector<double> AngleToAve;
    std::vector<double> ProjAlongAve;
    
    for(size_t i=0; i!=AllSeeds.size(); ++i)
      {
	ProjAlongAve.push_back(SeedPoss.at(i).Dot(AverageDir));
	AngleToAve.push_back(SeedDirs.at(i).Dot(AverageDir)/SeedDirs.at(i).Mag());
      }
  
    
    
    // For now we'll only look for 1 long track per cluster combo.
    //  There is room for improvement here if we need it.
    
    std::vector<size_t> TrackParts;
    std::vector<size_t> DoesNotFit;
    std::vector<size_t> AlreadyUsed;

     

    // This loop keeps running whilst there are still tracks to make
    while(AlreadyUsed.size()<AllSeeds.size())
      {
	       
	// This loop keeps running whilst there are still seeds
	//  to check for this track


	while((TrackParts.size()+DoesNotFit.size()+AlreadyUsed.size())<AllSeeds.size())
	  {
	    // First find the lowest track proj remaining
	    int LowestLeft=-1;
	    for(size_t i=0; i!=AllSeeds.size(); ++i)
	      {
		bool SeedUnavailable=false;
		for(size_t j=0; j!=TrackParts.size(); ++j)
		  if(i==TrackParts.at(j)) SeedUnavailable=true;
		for(size_t j=0; j!=DoesNotFit.size(); ++j)
		  if(i==DoesNotFit.at(j)) SeedUnavailable=true;
		for(size_t j=0; j!=AlreadyUsed.size(); ++j)
		  if(i==AlreadyUsed.at(j)) SeedUnavailable=true;
	    
		if(!SeedUnavailable)
		  {
		    if(LowestLeft < 0) LowestLeft = i;
		    if(ProjAlongAve.at(i)<ProjAlongAve.at(LowestLeft))
		      LowestLeft = i;		
		  }
	      }
	    
	    if(TrackParts.size()==0) TrackParts.push_back(LowestLeft);
	    else if(ConnectionMap[TrackParts.at(TrackParts.size()-1)][LowestLeft]!=0)
	      {
		
		// if track already running, check this comes out the right end
		if(TrackParts.size()>2)
		  {
		    if(ConnectionMap[TrackParts.at(TrackParts.size()-1)][LowestLeft] *
		       ConnectionMap[TrackParts.at(TrackParts.size()-1)][TrackParts.at(TrackParts.size()-2)] < 0)
		      {
			TrackParts.push_back(LowestLeft);
		      }		  
		    else
		      DoesNotFit.push_back(LowestLeft);
		  }
		else
		  TrackParts.push_back(LowestLeft);
	      }
	    else
	      DoesNotFit.push_back(LowestLeft);
	    
	    
	  }
	if(TrackParts.size()>0)
	  {
	    mf::LogVerbatim("BezierTrackerAlgorithm")<<" Prototrack runs through seeds: " ;
	    std::vector<recob::Seed> TheTrackSeeds;
	    for(size_t i=0; i!=TrackParts.size(); ++i)
	      {
		mf::LogVerbatim("BezierTrackerAlgorithm")<<"  " << TrackParts.at(i)<<", ";
		TheTrackSeeds.push_back(AllSeeds.at(TrackParts.at(i)));
	      }
	    OrganizedByTrack.push_back(TheTrackSeeds);  
	  }
	
	for(size_t j=0; j!=TrackParts.size(); ++j)
	  AlreadyUsed.push_back(TrackParts.at(j));
	TrackParts.clear();
	DoesNotFit.clear();
      }
	
    // Sort vector so most seeded tracks appear first
    std::sort(OrganizedByTrack.begin(), OrganizedByTrack.end(), BTrack_SeedCountComparator);
    return OrganizedByTrack;
  }
  


  //-----------------------------------------------
  bool BTrack_SeedCountComparator(std::vector<recob::Seed> const& s1, std::vector<recob::Seed> const& s2)
  {
    return s1.size()>s2.size();
  }


}
