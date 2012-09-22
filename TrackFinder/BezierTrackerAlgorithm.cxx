//
// Name: BezierTrackerAlgorithm.cxx
//
// Purpose: Implementation file for module BezierTrackerAlgorithm.
//
// Ben Jones, MIT
//

#include <vector>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TrackFinder/BezierTrackerAlgorithm.h"
#include "TrackFinder/SeedFinderAlgorithm.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "RecoBase/Hit.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Track.h"
#include "RecoBase/Prong.h"
#include "TrackFinder/BezierTrack.h"
#include "Utilities/AssociationUtil.h"


namespace trkf {

  BezierTrackerAlgorithm::BezierTrackerAlgorithm(const fhicl::ParameterSet& pset)
  {
    reconfigure(pset);
    fTopTrackID=0;
  
  }

  BezierTrackerAlgorithm::~BezierTrackerAlgorithm()
  {
  }

  void BezierTrackerAlgorithm::reconfigure(fhicl::ParameterSet const& pset)
  {

    fMaxKinkAngle      = pset.get<double>("MaxKinkAngle");
    fMaxTrackMissAngle = pset.get<double>("MaxTrackMissAngle");
    fMaxJumpDistance   = pset.get<double>("MaxJumpDistance");
    fHitDistance       = pset.get<double>("HitDistance");

    fTheSeedFinder = new SeedFinderAlgorithm(pset.get<fhicl::ParameterSet>("SeedFinder"));
    
  }




  std::vector<trkf::BezierTrack* > BezierTrackerAlgorithm::MakeBezierTracksFromSeeds(std::vector<recob::Seed> const& TrackSeeds )
  {
    
    mf::LogInfo("BezierTrackerAlgorithm")<<"Making bezier tracks from seeds"<<std::endl;
    std::vector<trkf::BezierTrack*> ReturnVector;
    
    std::vector<std::vector<recob::Seed> > OrgSeeds = OrganizeSeedsIntoTracks(TrackSeeds);
	
    mf::LogInfo("BezierTrackerAlgorithm")<<"Bezier Tracker OrgSeeds size : " << OrgSeeds.size()<<std::endl;
    // Loop through track candidate seed collections
    
    for(unsigned int i=0; i!=OrgSeeds.size(); i++)
      {
	mf::LogInfo("BezierTrackerAlgorithm")<<"Seeds in this btrack : " << OrgSeeds.at(i).size()<<std::endl;
	
	// Make a base track using seed positions and directions
	//  and build a BezierTrack on top of it.
	
	BezierTrack *BTrack = ProduceTrackFromSeeds(OrgSeeds.at(i));
	ReturnVector.push_back(BTrack);
      }

    return ReturnVector;

  }




  std::vector<trkf::BezierTrack* > BezierTrackerAlgorithm::MakeBezierTracksFromHits(std::vector<art::Ptr<recob::Hit> > HitVec, std::vector<art::PtrVector<recob::Hit> >& HitsForAssns )
  {
    mf::LogInfo("BezierTrackerAlgorithm")<<"Making bezier tracks from hits"<<std::endl;
  
    std::vector<trkf::BezierTrack*> ReturnVector;
    
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
	
	std::vector<recob::Seed> TrackSeeds = fTheSeedFinder->FindSeeds(SPVec, SPUsed);
	   	    
	mf::LogInfo("BezierTrackerAlgorithm")<<"Organizing seed collections " <<std::endl;
	// Organize these seeds into tracklike collections
	std::vector<std::vector<recob::Seed > > OrgSeeds = OrganizeSeedsIntoTracks(TrackSeeds);
	
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
	    BezierTrack *BTrack = ProduceTrackFromSeeds(OrgSeeds.at(i));
	 	
	    art::PtrVector<recob::Hit>         HitsThisTrack;
	    std::vector<double>                SValues;
	    
	    std::vector<int> HitIDs = DetermineNearbyHits(HitsToProcess, BTrack, fHitDistance, SValues);
	    
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
    
    return ReturnVector;
  }


  



  // Produce a track using the seeds we found
  // Note : we also need to add new seeds at beginning and end to get
  //  entire length of track extremities
  BezierTrack * BezierTrackerAlgorithm::ProduceTrackFromSeeds(std::vector<recob::Seed> const& Seeds)
  {

    std::vector<TVector3> Pos;
    std::vector<TVector3> Dir;
    std::vector<std::vector<double > > dQdx;
    
    double pt[3], dir[3], dummy[3];
    for(unsigned int i=0; i!=Seeds.size(); i++)
      {
	Seeds.at(i).GetPoint(     pt,  dummy );
	Seeds.at(i).GetDirection( dir, dummy );
	/*
	if(i==0)  // put in beginning seed
	  {
	    double ExtraPoint[3];
	    for(int j=0; j!=3; ++j)
	      ExtraPoint[j]=pt[j]+dir[j];
	    TVector3 ThisPos(ExtraPoint[0],ExtraPoint[1],ExtraPoint[2]);
	    TVector3 ThisDir(dir[0],dir[1],dir[2]);
	    
	    Pos.push_back(ThisPos);
	    Dir.push_back(ThisDir);
	  }
	*/
	TVector3 ThisPos(pt[0],  pt[1],  pt[2]  );
	TVector3 ThisDir(dir[0], dir[1], dir[2] );

	Pos.push_back(ThisPos);
	Dir.push_back(ThisDir);
	/*
	if(i==Seeds.size()-1)  // put in end seed
	  {
	    double ExtraPoint[3];
	    for(int j=0; j!=3; ++j)
	      ExtraPoint[j]=pt[j]-dir[j];
	    TVector3 ThisPos(ExtraPoint[0],ExtraPoint[1],ExtraPoint[2]);
	    TVector3 ThisDir(dir[0],dir[1],dir[2]);
	    
	    Pos.push_back(ThisPos);
	    Dir.push_back(ThisDir);
	  }
	*/
      }

    BezierTrack * TheTrack = new BezierTrack(Pos,Dir,dQdx,fTopTrackID++);
    return TheTrack;
  }


  //----------------------------------------------------------------------
  //
  // From a PtrVector of hits, determine which are nearby
  //

  std::vector<int> BezierTrackerAlgorithm::DetermineNearbyHits(art::PtrVector<recob::Hit> Hits, BezierTrack* BTrack, double HitCollectionDistance, std::vector<double>& SValues)
  {
    std::vector<int> ReturnVector;
    std::vector<double> s, distance;
    std::cout<<"Size of hit vector : " <<Hits.size()<<std::endl;
    BTrack->GetClosestApproaches(Hits,SValues, distance);

    for(size_t i=0; i!=Hits.size(); ++i)
      {

	if((distance.at(i)<HitCollectionDistance)&&(SValues.at(i)<1)&&(SValues.at(i)>0)) ReturnVector.push_back(i); 
	//       	std::cout<< i <<"  "<< Hits.at(i)->View()<<"  "<<distance.at(i)<<std::endl;
      }
    return ReturnVector;
  }


  

  //----------------------------------------------------------------------
  std::vector<std::vector< recob::Seed > > BezierTrackerAlgorithm::OrganizeSeedsIntoTracks(std::vector<recob::Seed > TrackSeeds)
  {
    std::vector<std::vector<recob::Seed > > OrganizedByTrack;

    while(TrackSeeds.size()>1)
      {
	// Declare this track
	std::vector<recob::Seed > ThisTrack;
	
        // Add first element
	std::vector<recob::Seed >::iterator ifirst= TrackSeeds.begin();
        ThisTrack.push_back(*ifirst);
        TrackSeeds.erase(ifirst);
	
        for(std::vector<recob::Seed >::iterator it=TrackSeeds.begin();
	    it!=TrackSeeds.end(); it++)
	  {
	    recob::Seed LastSeedAdded = ThisTrack.at(ThisTrack.size()-1);
	    recob::Seed ThisSeed = *it;
	    
            // Does this seed fit with this track?
	    double SeedPos[3],Err[3];

	    ThisSeed.GetPoint(SeedPos,Err);
	    

	    //  float Angle = ThisSeed->GetAngle(*LastSeedAdded);
	    //  float ProjDis = ThisSeed->GetProjAngleDiscrepancy(*LastSeedAdded);
	    //  std::cout<<"BezierTrackerAlgorithm: " << Angle << " " << ProjDis << "  for seed : "; 
	    //  std::cout<<double(SeedPos[0]) << " " << double(SeedPos[1]) << " " << double(SeedPos[2]) << std::endl;
            if((  abs(ThisSeed.GetAngle(LastSeedAdded))                 < fMaxKinkAngle)
	       &&( abs(ThisSeed.GetProjAngleDiscrepancy(LastSeedAdded))  < fMaxTrackMissAngle)
	       &&( abs(ThisSeed.GetDistance(LastSeedAdded)<fMaxJumpDistance)))
                {
                  // if so, add it into the track, erase it from the stack
                  // and start looping at the beginning again
                  ThisTrack.push_back(ThisSeed);
                  TrackSeeds.erase(it);
		  it=TrackSeeds.begin();
                }
          }
        // We ran out of seeds on the stack. store this track and go again
	if(ThisTrack.size()>1)
	  OrganizedByTrack.push_back(ThisTrack);
      }   
    return OrganizedByTrack;
  }


}
