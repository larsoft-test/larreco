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
#include "Utilities/AssociationUtil.h"

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

    fMaxJumpLengths    = pset.get<double>("MaxJumpLengths");
    fHitDistance       = pset.get<double>("HitDistance");
    fTrackJoinAngle    = pset.get<double>("TrackJoinAngle");
    fDirectJoinDistance= pset.get<double>("DirectJoinDistance");
    
    fVertexImpactThreshold = pset.get<double>("VertexImpactThreshold");
    fVertexExtrapDistance = pset.get<double>("VertexExtrapDistance");
    

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
	mf::LogVerbatim("BezierTrackerAlgorithm")<<"Getting seeds " <<std::endl;

	std::vector<art::PtrVector<recob::Hit> > HitCatalogue;	
	std::vector<recob::Seed> AllSeeds = fTheSeedFinder->GetSeedsFromUnSortedHits(HitsToProcess, HitCatalogue);
	   	    
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
