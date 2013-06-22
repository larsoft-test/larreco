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

  BezierTrackerAlgorithm::BezierTrackerAlgorithm(const fhicl::ParameterSet& pset)
  {
    reconfigure(pset);
    fTopTrackID=0;
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

    fTheSeedFinder = new SeedFinderAlgorithm(pset.get<fhicl::ParameterSet>("SeedFinder"));
    
  }


  //------------------------------------------------

  void BezierTrackerAlgorithm::MakeBezierTracksFromSeeds(std::vector<trkf::BezierTrack>& ReturnVector, std::vector<recob::Seed> const& AllSeeds )
  {
    
    mf::LogInfo("BezierTrackerAlgorithm")<<"Making bezier tracks from seeds"<<std::endl;
    
    std::vector<std::vector<recob::Seed> > OrgSeeds = OrganizeSeedsIntoTracks(AllSeeds);

    mf::LogInfo("BezierTrackerAlgorithm")<<"Producing track objects"<<std::endl;    

    for(unsigned int i=0; i!=OrgSeeds.size(); i++)
      {
	mf::LogInfo("BezierTrackerAlgorithm")<<"Seeds in this btrack : " << OrgSeeds.at(i).size()<<std::endl;
	if(OrgSeeds.at(i).size()>0)
	  ReturnVector.push_back(trkf::BezierTrack(OrgSeeds.at(i)));
      }

  }


  //------------------------------------------------

  void BezierTrackerAlgorithm::MakeBezierTracksFromHits(std::vector<trkf::BezierTrack>& ReturnVector, std::vector<art::Ptr<recob::Hit> > HitVec, std::vector<art::PtrVector<recob::Hit> >& HitsForAssns )
  {
    mf::LogInfo("BezierTrackerAlgorithm")<<"Making bezier tracks from hits"<<std::endl;
    
    
    // This vector keeps track of which hits we are still processing  
    art::PtrVector<recob::Hit> HitsToProcess;
    for(size_t i=0; i!=HitVec.size(); ++i)
      {
	HitsToProcess.push_back(HitVec.at(i));
      }

    bool KeepTrying=true;
    while(KeepTrying)
      {
	mf::LogInfo("BezierTrackerAlgorithm")<<"Getting space points" <<std::endl;
	
	// Make remaining hits into SPs
	std::vector<recob::SpacePoint> SPVec = 
	  fTheSeedFinder->GetSpacePointsFromHitVector(HitsToProcess);
	    
	// Find seeds in these SPs
	std::vector<std::vector<recob::SpacePoint> > SPUsed;
	mf::LogInfo("BezierTrackerAlgorithm")<<"Getting seeds " <<std::endl;
	
	std::vector<recob::Seed> AllSeeds = fTheSeedFinder->FindSeeds(SPVec, SPUsed);
	   	    
	mf::LogInfo("BezierTrackerAlgorithm")<<"Organizing seed collections " <<std::endl;
	// Organize these seeds into tracklike collections
	std::vector<std::vector<recob::Seed > > OrgSeeds = OrganizeSeedsIntoTracks(AllSeeds);
	
	// If we didn't find any, continue and quit
	if(OrgSeeds.size()==0) 
	  {
	    KeepTrying=false;
	    continue;
	  }
	    
	mf::LogInfo("BezierTrackerAlgorithm")<<"Making tracks from  " << OrgSeeds.size()<<" seed collections"<<std::endl;
	// For each of them make 1 track
	for(unsigned int i=0; i!=OrgSeeds.size(); i++)
	  {
	    
	    // Make a base track using seed positions and directions
	    //  and build a BezierTrack on top of it.		
	    BezierTrack BTrack(OrgSeeds.at(i));
	 	
	    art::PtrVector<recob::Hit>         HitsThisTrack;
	    std::vector<double>                SValues;
	    
	    std::vector<int> HitIDs = DetermineNearbyHits(HitsToProcess, BTrack, SValues);
	    
	    mf::LogInfo("BezierTrackerAlgorithm")<<"Found " << HitIDs.size()<<" nearby hits" <<std::endl;
	    
	    for(size_t i=0; i!=HitIDs.size(); i++)
	      HitsThisTrack.push_back(HitsToProcess.at(HitIDs.at(i)));
	    
	    HitsForAssns.push_back(HitsThisTrack);
	    
	    ReturnVector.push_back(BTrack);
	 
	    // Remove the hits we used from the vector and go around
	    mf::LogInfo("BezierTrackerAlgorithm")<<"Removing used hits"<<std::endl;
	    for(int i=HitIDs.size()-1; i!=-1; --i)
	      HitsToProcess.erase(HitsToProcess.begin() + HitIDs.at(i));
	  }
	mf::LogInfo("BezierTrackerAlgorithm")<<"Uncollected hit size : " << HitsToProcess.size()<<std::endl;
	
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
  // This is the key function of the bezier tracker.  It takes a set of seeds,
  //  and attempts to organize them into a set of tracks.
  //  Currently the functionality allows for finding 1 long track,
  //  and a set of single-seed leftover parts.  This lends itself
  //  to cluster combination based tracking.  We can revist and improve
  //  this method for multi-tracks in the future.


  std::vector<std::vector< recob::Seed > > BezierTrackerAlgorithm::OrganizeSeedsIntoTracks(std::vector<recob::Seed > const& AllSeeds)
  {
      
    std::vector<std::vector<recob::Seed > > OrganizedByTrack;
    std::map<int, std::map<int, double> > ConnectionMap;
    
    for(size_t i=0; i!=AllSeeds.size(); ++i)
      {
	for(size_t j=0; j!=AllSeeds.size(); ++j)
	  {
	    if(i==j) continue;
	    double Distance = AllSeeds.at(j).GetDistance(AllSeeds.at(i));
	    double Angle    = AllSeeds.at(j).GetAngle(AllSeeds.at(i));
	    double Length   = std::max(AllSeeds.at(i).GetLength(), AllSeeds.at(j).GetLength());
	    int PointingSign = AllSeeds.at(i).GetPointingSign(AllSeeds.at(j));
	     
	    if( ( (Angle/Distance)<fMaxKinkDThetaDx) && (((Distance-Length)/Length)<fMaxJumpLengths)) 
	      {
		
		ConnectionMap[i][j] = (Distance-Length)/Length*float(PointingSign);
	      
	      }      
	  }
      }
    
    std::vector<int> NForwards, NBackwards, BestForward, BestBackward;
    NForwards.resize(AllSeeds.size());
    NBackwards.resize(AllSeeds.size());
    BestForward.resize(AllSeeds.size());
    BestBackward.resize(AllSeeds.size());
    
    std::vector<int> StartCandidates;
    std::vector<int> EndCandidates;
    

    for(size_t i=0; i!=AllSeeds.size(); ++i)	
      {
	for(size_t j=0; j!=AllSeeds.size(); ++j)	
	  {
	    if(ConnectionMap[i][j]>0)
	      {
		NForwards[i]++;
		if( (ConnectionMap[i][j] < ConnectionMap[i][BestForward[i]])
		    || ConnectionMap[i][BestForward[i]]==0)
		  BestForward[i] = j;
	      }
	    else if(ConnectionMap[i][j]<0)
	      {
		NBackwards[i]++;
		if( (ConnectionMap[i][j] > ConnectionMap[i][BestBackward[i]])
		    || ConnectionMap[i][BestBackward[i]]==0)
		  BestBackward[i] = j;
	      }
	  }
	
	if((NForwards[i]>0)&&(NBackwards[i]==0))
	  StartCandidates.push_back(i);
	
	if((NForwards[i]==0)&&(NBackwards[i]>0))
	  EndCandidates.push_back(i);
	
	// If seed is totally unconnected, it is its own track.
	if((NForwards[i]==0)&&(NBackwards[i]==0))
	  {
	    OrganizedByTrack.push_back(std::vector<recob::Seed>() );
	    OrganizedByTrack.at(OrganizedByTrack.size()-1).push_back(AllSeeds.at(i));
		
	  }	
      }
    
    // For now we'll only look for 1 long track per cluster combo.
    //  There is room for improvement here if we need it.
    
    std::vector<int> TrackParts;
    
      
    // If clear start, we're going to go forwards
    if(StartCandidates.size()==1)
      {
	TrackParts.push_back(StartCandidates.at(0));		
	bool KeepGoing=true;
	while(KeepGoing)
	  {
	    if(NForwards[TrackParts.at(TrackParts.size()-1)]>0)
	      {
		int NextPart = BestForward[TrackParts.at(TrackParts.size()-1)];
		//		std::cout<<"going to " << NextPart<<std::endl;
		
		// Escape in case of a loop
		for(size_t p=0; p!=TrackParts.size(); ++p)
		  if(TrackParts.at(p)==NextPart) KeepGoing=false;
		TrackParts.push_back(NextPart);

	      }
	    else KeepGoing=false;   
	  }
      }
    // If no clear start, but a clear end, we're going to go backwards
    else if(EndCandidates.size()==1)
      {
	TrackParts.push_back(EndCandidates.at(0));
	bool KeepGoing=true;
	while(KeepGoing)
	  {
	    if(NBackwards[TrackParts.at(TrackParts.size()-1)]>0)
	      {
		int NextPart = BestBackward[TrackParts.at(TrackParts.size()-1)];
		//		std::cout<<"going to " << NextPart<<std::endl;
		
		// Escape in case of a loop
		for(size_t p=0; p!=TrackParts.size(); ++p)
		  if(TrackParts.at(p)==NextPart) KeepGoing=false;
		TrackParts.push_back(NextPart);
	      }
	    else KeepGoing=false;   
	  }
      }
    
    if(TrackParts.size()>0)
      {
	std::vector<recob::Seed> TheTrackSeeds;
	for(size_t i=0; i!=TrackParts.size(); ++i)
	  {
	    TheTrackSeeds.push_back(AllSeeds.at(TrackParts.at(i)));
	  }
	
	OrganizedByTrack.push_back(TheTrackSeeds);  
      }

     return OrganizedByTrack;
  }
  

}
