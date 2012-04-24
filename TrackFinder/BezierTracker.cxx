//
// Name: BezierTracker.cxx
//
// Purpose: Implementation file for module BezierTracker.
//
// Ben Jones, MIT
//

#include <vector>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TrackFinder/BezierTracker.h"
#include "TrackFinder/SeedFinder.h"
#include "Geometry/geo.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "RecoBase/Hit.h"
#include "RecoBase/Seed.h"
#include "RecoBase/Track.h"
#include "RecoBase/Prong.h"
#include "TrackFinder/BezierTrack.h"
#include "TrackFinder/SeedFinder.h"
#include "Utilities/AssociationUtil.h"


namespace trkf {

  BezierTracker::BezierTracker(const fhicl::ParameterSet& pset)
  {
    reconfigure(pset);
    produces< std::vector<recob::Track> >();
    produces< std::vector<recob::Hit> >();
    produces< art::Assns<recob::Track, recob::Hit> >();
    fTopTrackID=0;
  
  }

  BezierTracker::~BezierTracker()
  {
  }

  void BezierTracker::reconfigure(fhicl::ParameterSet const& pset)
  {
    fSeedModuleLabel   = pset.get<std::string>("SeedModuleLabel");
    fHitModuleLabel    = pset.get<std::string>("HitModuleLabel");
    fMaxKinkAngle      = pset.get<double>("MaxKinkAngle");
    fMaxTrackMissAngle = pset.get<double>("MaxTrackMissAngle");
    fMaxJumpDistance   = pset.get<double>("MaxJumpDistance");
    fHitDistance       = pset.get<double>("HitDistance");
    fTrackMode         = pset.get<double>("TrackMode");
    
    if((fTrackMode==2)||(fTrackMode==3)); 
      {
	fTheSeedFinder = new SeedFinder(pset.get<fhicl::ParameterSet>("SeedFinder"));
	
      }  
}

  void BezierTracker::beginJob()
  {}


  void BezierTracker::produce(art::Event& evt)
  {
 
    // Extract hits PtrVector from event

    art::Handle< std::vector<recob::Hit> > hith;
    evt.getByLabel(fHitModuleLabel, hith);

    art::PtrVector<recob::Hit> HitVec;
    for(unsigned int i=0; i < hith->size(); ++i)
      {
	art::Ptr<recob::Hit> hit(hith,i);
	HitVec.push_back(hit);
      }
    

    // Declare products to store

    std::auto_ptr< std::vector<recob::Track > > btracks ( new std::vector<recob::Track>);
    std::auto_ptr< std::vector<recob::Hit > > leftovers ( new std::vector<recob::Hit>);
    std::auto_ptr< art::Assns<recob::Track, recob::Hit > > assn( new art::Assns<recob::Track, recob::Hit>);
   
    if(fTrackMode==1)
      {
	std::map<int,bool> UsedHitIDs;

	// Load a vector of track seeds from the event
      	art::Handle< std::vector<recob::Seed> > seedh;
	evt.getByLabel(fSeedModuleLabel, seedh);  
	
	std::vector<art::Ptr<recob::Seed> > TrackSeeds;
	if(seedh->size()>0)
	  for(unsigned int iseed=0; iseed!=seedh->size(); iseed++)
	    {
	      art::Ptr<recob::Seed> theseed(seedh, iseed);
	      TrackSeeds.push_back(theseed);
	    }
	
	// Organize these seeds into track candidates based on
	// proximity and pointing
	
	std::vector<std::vector<art::Ptr<recob::Seed> > > OrgSeeds = OrganizeSeedsIntoTracks(TrackSeeds);
	
	std::cout<<"Bezier Tracker OrgSeeds size : " << OrgSeeds.size()<<std::endl;
	// Loop through track candidate seed collections
	
	for(unsigned int i=0; i!=OrgSeeds.size(); i++)
	  {
	    std::cout<<"Seeds in this btrack : " << OrgSeeds.at(i).size()<<std::endl;
	    
	    // Make a base track using seed positions and directions
	    //  and build a BezierTrack on top of it.
	    
	    BezierTrack *BTrack = ProduceTrackFromSeeds(OrgSeeds.at(i));
	    
	    // Collect hits in the vicinity of this track
	    
	    std::cout<<"Determining nearby hits for bezier track"<<std::endl;
	    std::vector<double> SValues;
	    std::vector<int> HitIDs = DetermineNearbyHits(HitVec, BTrack, fHitDistance, SValues);
	    art::PtrVector<recob::Hit> HitsToAssoc;	
	    
	    std::cout<<"Finding hits "<<std::endl;
	    std::vector<double> SValuesFordQdx;
	    for(size_t i=0; i!=HitIDs.size(); i++)
	      {
		SValuesFordQdx.push_back(SValues.at(i));
		UsedHitIDs[HitIDs.at(i)]=true;
		HitsToAssoc.push_back(HitVec.at(HitIDs.at(i)));
	      }
	    
	    // Using the hits, fill in dQdx in the track
	    std::cout<<"Calculating dQdx"<<std::endl;
	    BTrack->CalculatedQdx(HitsToAssoc,SValuesFordQdx);
	    BTrack->FillMySpacePoints(1000);
	    
	    //Put the base object into the auto_ptr and make associations
	    std::cout<<"Inserting into vector and creating assn"<<std::endl;
	   
	    recob::Track TheTrack = BTrack->GetBaseTrack();
	    
	    

	    btracks->push_back(TheTrack);       
	    util::CreateAssn(*this, evt, *(btracks.get()), HitsToAssoc, *(assn.get())); 
	    
	    for(size_t i=0; i!=HitVec.size(); i++)
	      {
		if(UsedHitIDs[i]!=true) leftovers->push_back(*HitVec.at(i));
	      }
	  }
      }

   
    else if(fTrackMode==2)
      {
	std::cout<<"Ready to try to do iterative tracking"<<std::endl;

	art::PtrVector<recob::Hit> HitsToProcess = HitVec;

	// In this mode we have to figure it all out first
	//  and store late
	std::vector<art::PtrVector<recob::Hit > > HitsForBTracks;
	std::vector<trkf::BezierTrack*> TracksToStore;
	
	bool KeepTrying=true;
	while(KeepTrying)
	  {
	    std::cout<<"Getting space points" <<std::endl;
	    
	    // Make remaining hits into SPs
	    std::vector<recob::SpacePoint> SPVec = 
	      fTheSeedFinder->GetSpacePointsFromHitVector(HitsToProcess);
	    
	    // Find seeds in these SPs
	    std::vector<std::vector<recob::SpacePoint> > SPUsed;
	    std::cout<<"Getting seeds " <<std::endl;
	    std::vector<recob::Seed*> TrackSeeds = fTheSeedFinder->FindAsManySeedsAsPossible(SPVec, SPUsed);
	    std::cout<<"Beginning iterative refitting " <<std::endl;
	    for(size_t i=0; i!=TrackSeeds.size(); i++)
	      {
		//		fTheSeedFinder->RefitSeed(TrackSeeds.at(i),SPUsed.at(i));
	      }
	    
	    std::cout<<"Organizing seed collections " <<std::endl;
	    // Organize these seeds into tracklike collections
	    std::vector<std::vector<recob::Seed* > > OrgSeeds = OrganizeSeedsIntoTracks(TrackSeeds);
	    
	    // If we didn't find any, continue and quit
	    if(OrgSeeds.size()==0) 
	      {
		KeepTrying=false;
		continue;
	      }
	    
	    std::cout<<"Making tracks from  " << OrgSeeds.size()<<" seed collections"<<std::endl;
	    // For each of them make 1 track
	    for(unsigned int i=0; i!=OrgSeeds.size(); i++)
	      {
		    
		// Make a base track using seed positions and directions
		//  and build a BezierTrack on top of it.		
		BezierTrack *BTrack = ProduceTrackFromSeeds(OrgSeeds.at(i));
		
		// Collect hits in the vicinity of this track
		std::vector<double> SValues;
		std::vector<int> HitIDs = DetermineNearbyHits(HitsToProcess, BTrack, fHitDistance, SValues);
	   
		art::PtrVector<recob::Hit> HitsToAssoc;	
		HitsToAssoc.clear();
		std::cout<<"seeking hits for associations"<<std::endl;
		std::vector<double> SValuesFordQdx;
		std::cout<<"Found " << HitIDs.size()<<" nearby hits" <<std::endl;
		for(size_t i=0; i!=HitIDs.size(); i++)
		  {
		    HitsToAssoc.push_back(HitsToProcess.at(HitIDs.at(i)));
		    SValuesFordQdx.push_back(SValues.at(i));
		  }
		
		HitsForBTracks.push_back(HitsToAssoc);
	    
		// Using the hits, make SPs and fill dQdx in the track
		std::cout<<"filling spacepoints"<<std::endl;
		BTrack->FillMySpacePoints(100);
		std::cout<<"Calulating dqdx"<<std::endl;
		BTrack->CalculatedQdx(HitsToAssoc,SValuesFordQdx);
	    
		//Put the base object into a vector for storage later
		TracksToStore.push_back(BTrack); 
			
		// Remove the hits we used from the vector and go around
		std::cout<<"Removing used hits"<<std::endl;
		for(int  i=HitIDs.size()-1; i!=-1; --i)
		  {
		    HitsToProcess.erase(HitsToProcess.begin() + HitIDs.at(i));
		  }
		std::cout<<"Uncollected hit size : " << HitsToProcess.size()<<std::endl;
		if(HitsToProcess.size()<3) KeepTrying=false;
		
	      }	    
	  }

	// Associate hits to the tracks we found
	for(size_t i=0; i!=TracksToStore.size(); i++)
	  {
	    recob::Track TheTrack = TracksToStore.at(i)->GetBaseTrack();
	    TheTrack.SetID(TracksToStore.at(i)->ID());
	    btracks->push_back(TheTrack);
	    util::CreateAssn(*this, evt, *(btracks.get()), HitsForBTracks.at(i), *(assn.get())); 
	  }
	for(size_t i=0; i!=HitsToProcess.size(); i++)
	  {
	    leftovers->push_back(*HitsToProcess.at(i));
	  }

      }

    else if(fTrackMode==3)
      {
	std::cout<<"Ready to try to do iterative tracking"<<std::endl;

	art::PtrVector<recob::Hit> HitsToProcess = HitVec;

	std::vector<recob::Seed*> AllTheSeeds;

	bool KeepTrying=true;
	while(KeepTrying)
	  {
	    
	    // Make remaining hits into SPs
	    std::vector<recob::SpacePoint> SPVec = 
	      fTheSeedFinder->GetSpacePointsFromHitVector(HitsToProcess);
	    
	    std::vector<std::vector<recob::SpacePoint> > SPUsed;

	    // Find seeds in these SPs
	    std::vector<recob::Seed*> TrackSeeds = fTheSeedFinder->FindAsManySeedsAsPossible(SPVec, SPUsed);
	    std::cout<<"Filling vector"<<std::endl;
	    AllTheSeeds.insert(AllTheSeeds.end(), TrackSeeds.begin(), TrackSeeds.end());
	    
	   
	    // Organize these seeds into tracklike collections
	    std::vector<std::vector<recob::Seed* > > OrgSeeds = OrganizeSeedsIntoTracks(TrackSeeds);
	    
	    // If we didn't find any, continue and quit
	    if(OrgSeeds.size()==0) 
	      {
		KeepTrying=false;
		continue;
	      }

	    // For each of them make 1 track
	    for(unsigned int i=0; i!=OrgSeeds.size(); i++)
	      {
		    
		// Make a base track using seed positions and directions
		//  and build a BezierTrack on top of it.		
		BezierTrack *BTrack = ProduceTrackFromSeeds(OrgSeeds.at(i));
		
		// Collect hits in the vicinity of this track	    
		std::vector<double> SValues;
		std::vector<int> HitIDs = DetermineNearbyHits(HitsToProcess, BTrack, fHitDistance, SValues);
	   
		// Remove the hits we used from the vector and go around
		for(int  i=HitIDs.size()-1; i!=-1; --i)
		  {
		    HitsToProcess.erase(HitsToProcess.begin() + HitIDs.at(i));
		  }
		if(HitsToProcess.size()<3) KeepTrying=false;
		delete BTrack;
		
	      }    
	  }
	std::cout<<"Produced all seeds"<<std::endl;

	// Now we have fully seeded, try to make best tracks we can

	std::vector<std::vector<recob::Seed* > > OrgSeeds = OrganizeSeedsIntoTracks(AllTheSeeds);
	
	// Loop through track candidate seed collections
	for(unsigned int i=0; i!=OrgSeeds.size(); i++)
	  {
	    // Make a base track using seed positions and directions
	    //  and build a BezierTrack on top of it.
	    BezierTrack *BTrack = ProduceTrackFromSeeds(OrgSeeds.at(i));
	    
	    // Collect hits in the vicinity of this track
	    std::vector<double> SValues;
	    std::vector<int> HitIDs = DetermineNearbyHits(HitVec, BTrack, fHitDistance,SValues);
	    art::PtrVector<recob::Hit> HitsToAssoc;	
	    
	    std::vector<double> SValuesFordQdx;
	    for(size_t i=0; i!=HitIDs.size(); i++)
	      {
		SValuesFordQdx.push_back(SValues.at(i));
		HitsToAssoc.push_back(HitVec.at(HitIDs.at(i)));
	      }
	    
	    // Using the hits, fill in dQdx in the track
	    std::cout<<"Calculating dQdx"<<std::endl;
	    BTrack->CalculatedQdx(HitsToAssoc,SValuesFordQdx);
	    BTrack->FillMySpacePoints(1000);
	    
	    //Put the base object into the auto_ptr and make associations
	    recob::Track TheTrack = BTrack->GetBaseTrack();
	 
	    btracks->push_back(TheTrack);
	    
	    util::CreateAssn(*this, evt, *(btracks.get()), HitsToAssoc, *(assn.get())); 
	    
	  }
	for(size_t i=0; i!=HitsToProcess.size(); i++)
	  {
	    leftovers->push_back(*HitsToProcess.at(i));
	  }

      }




    // Store the fruits of our labour in the event
    std::cout<<"Storing in evt"<<std::endl;
    evt.put(btracks);
    evt.put(assn);
    evt.put(leftovers);
  }



  // Produce a track using the seeds we found
  // Note : we also need to add new seeds at beginning and end to get
  //  entire length of track extremities
  BezierTrack * BezierTracker::ProduceTrackFromSeeds(std::vector<art::Ptr<recob::Seed> > Seeds)
  {

    std::vector<TVector3> Pos;
    std::vector<TVector3> Dir;
    std::vector<std::vector<double > > dQdx;
    
    double pt[3], dir[3], dummy[3];
    for(unsigned int i=0; i!=Seeds.size(); i++)
      {
	Seeds.at(i)->GetPoint(     pt,  dummy );
	Seeds.at(i)->GetDirection( dir, dummy );

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
	
	TVector3 ThisPos(pt[0],  pt[1],  pt[2]  );
	TVector3 ThisDir(dir[0], dir[1], dir[2] );

	Pos.push_back(ThisPos);
	Dir.push_back(ThisDir);

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

      }
    BezierTrack * TheTrack = new BezierTrack(Pos,Dir,dQdx);
    TheTrack->SetID(fTopTrackID++);
    return TheTrack;
  }

  BezierTrack * BezierTracker::ProduceTrackFromSeeds(std::vector<recob::Seed* > Seeds)
  {

    std::vector<TVector3> Pos;
    std::vector<TVector3> Dir;
    std::vector<std::vector<double > > dQdx;
    
    double pt[3], dir[3], dummy[3];
    for(unsigned int i=0; i!=Seeds.size(); i++)
      {
	Seeds.at(i)->GetPoint(     pt,  dummy );
	Seeds.at(i)->GetDirection( dir, dummy );

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
	
	TVector3 ThisPos(pt[0],  pt[1],  pt[2]  );
	TVector3 ThisDir(dir[0], dir[1], dir[2] );

	Pos.push_back(ThisPos);
	Dir.push_back(ThisDir);

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

      }
    BezierTrack * TheTrack = new BezierTrack(Pos,Dir,dQdx);
    TheTrack->SetID(fTopTrackID++);
    return TheTrack;
  }


  //
  // From a PtrVector of hits, determine which are nearby
  //

  std::vector<int> BezierTracker::DetermineNearbyHits(art::PtrVector<recob::Hit> Hits, BezierTrack* BTrack, double HitCollectionDistance, std::vector<double>& SValues)
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


  

  std::vector<std::vector<art::Ptr< recob::Seed> > > BezierTracker::OrganizeSeedsIntoTracks(std::vector<art::Ptr<recob::Seed> > TrackSeeds)
  {
    std::vector<std::vector<art::Ptr<recob::Seed> > > OrganizedByTrack;

    while(TrackSeeds.size()>1)
      {
	// Declare this track
	std::vector<art::Ptr<recob::Seed> > ThisTrack;
	ThisTrack.clear();

        // Add first element
	std::vector<art::Ptr<recob::Seed> >::iterator ifirst= TrackSeeds.begin();
        ThisTrack.push_back(*ifirst);
        TrackSeeds.erase(ifirst);
	
        for(std::vector<art::Ptr<recob::Seed> >::iterator it=TrackSeeds.begin();
	    it!=TrackSeeds.end(); it++)
	  {
	    art::Ptr<recob::Seed> LastSeedAdded = ThisTrack.at(ThisTrack.size()-1);
	    art::Ptr<recob::Seed> ThisSeed = *it;
	    
            // Does this seed fit with this track?
	    double SeedPos[3],Err[3];
	    float Angle = ThisSeed->GetAngle(*LastSeedAdded);
	    float ProjDis = ThisSeed->GetProjAngleDiscrepancy(*LastSeedAdded);
	    ThisSeed->GetPoint(SeedPos,Err);
	    
	    std::cout<<"BezierTracker: " << Angle << " " << ProjDis << "  for seed : "; 
	    std::cout<<double(SeedPos[0]) << " " << double(SeedPos[1]) << " " << double(SeedPos[2]) << std::endl;
            if((  abs(ThisSeed->GetAngle(*LastSeedAdded))                 < fMaxKinkAngle)
	       &&( abs(ThisSeed->GetProjAngleDiscrepancy(*LastSeedAdded))  < fMaxTrackMissAngle)
	       &&( abs(ThisSeed->GetDistance(*LastSeedAdded)<fMaxJumpDistance)))
                {
                  // if so, add it into the track, erase it from the stack
                  // and start looping at the beginning again
                  ThisTrack.push_back(ThisSeed);
                  TrackSeeds.erase(it);
		  it=TrackSeeds.begin();
                }
          }
        // We ran out of seeds on the stack. store this track and go again
	if(ThisTrack.size()>2)
	  OrganizedByTrack.push_back(ThisTrack);
      }   
    return OrganizedByTrack;
  }




  std::vector<std::vector< recob::Seed* > > BezierTracker::OrganizeSeedsIntoTracks(std::vector<recob::Seed* > TrackSeeds)
  {
    std::vector<std::vector<recob::Seed* > > OrganizedByTrack;

    while(TrackSeeds.size()>1)
      {
	// Declare this track
	std::vector<recob::Seed* > ThisTrack;
	
        // Add first element
	std::vector<recob::Seed* >::iterator ifirst= TrackSeeds.begin();
        ThisTrack.push_back(*ifirst);
        TrackSeeds.erase(ifirst);
	
        for(std::vector<recob::Seed* >::iterator it=TrackSeeds.begin();
	    it!=TrackSeeds.end(); it++)
	  {
	    recob::Seed* LastSeedAdded = ThisTrack.at(ThisTrack.size()-1);
	    recob::Seed* ThisSeed = *it;
	    
            // Does this seed fit with this track?
	    double SeedPos[3],Err[3];
	    float Angle = ThisSeed->GetAngle(*LastSeedAdded);
	    float ProjDis = ThisSeed->GetProjAngleDiscrepancy(*LastSeedAdded);
	    ThisSeed->GetPoint(SeedPos,Err);
	    
	    std::cout<<"BezierTracker: " << Angle << " " << ProjDis << "  for seed : "; 
	    std::cout<<double(SeedPos[0]) << " " << double(SeedPos[1]) << " " << double(SeedPos[2]) << std::endl;
            if((  abs(ThisSeed->GetAngle(*LastSeedAdded))                 < fMaxKinkAngle)
	       &&( abs(ThisSeed->GetProjAngleDiscrepancy(*LastSeedAdded))  < fMaxTrackMissAngle)
	       &&( abs(ThisSeed->GetDistance(*LastSeedAdded)<fMaxJumpDistance)))
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


  void BezierTracker::endJob()
  {

  }
}
