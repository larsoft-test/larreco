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
    
    if(fTrackMode==2) fTheSeedFinder = new SeedFinder(pset);
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
    std::auto_ptr< art::Assns<recob::Track, recob::Hit > > assn( new art::Assns<recob::Track, recob::Hit>);
   

    if(fTrackMode==1)
      {
	
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
	
	
	// Loop through track candidate seed collections
	
	for(unsigned int i=0; i!=OrgSeeds.size(); i++)
	  {
	    std::cout<<"Seeds in this btrack : " << OrgSeeds.at(i).size()<<std::endl;
	    
	    // Make a base track using seed positions and directions
	    //  and build a BezierTrack on top of it.
	    
	    BezierTrack *BTrack = ProduceTrackFromSeeds(OrgSeeds.at(i));
	    
	    // Collect hits in the vicinity of this track
	    
	    std::cout<<"Determining nearby hits for bezier track"<<std::endl;
	    std::vector<int> HitIDs = DetermineNearbyHits(HitVec, BTrack, fHitDistance);
	    art::PtrVector<recob::Hit> HitsToAssoc;	
	    
	    std::cout<<"Finding hits "<<std::endl;
	    for(size_t i=0; i!=HitIDs.size(); i++)
	      {
		HitsToAssoc.push_back(HitVec.at(HitIDs.at(i)));
	      }
	    
	    // Using the hits, fill in dQdx in the track
	    std::cout<<"Calculating dQdx"<<std::endl;
	    BTrack->CalculatedQdx(HitsToAssoc);
	    
	    
	    //Put the base object into the auto_ptr and make associations
	    std::cout<<"Inserting into vector and creating assn"<<std::endl;
	    recob::Track TheTrack = BTrack->GetBaseTrack();
	    
	    btracks->push_back(TheTrack);       
	    util::CreateAssn(*this, evt, *(btracks.get()), HitsToAssoc, *(assn.get())); 
	    
	  }
      }
    else if(fTrackMode==2)
      {
	std::cout<<"Ready to try to do iterative tracking"<<std::endl;

      }

    // Store the fruits of our labour in the event
    std::cout<<"Storing in evt"<<std::endl;
    evt.put(btracks);
    evt.put(assn);
  }

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
	
	TVector3 ThisPos(pt[0],  pt[1],  pt[2]  );
	TVector3 ThisDir(dir[0], dir[1], dir[2] );

	Pos.push_back(ThisPos);
	Dir.push_back(ThisDir);
      }
    BezierTrack * TheTrack = new BezierTrack(Pos,Dir,dQdx);
    TheTrack->SetID(fTopTrackID++);
    return TheTrack;
  }


  //
  // From a PtrVector of hits, determine which are nearby
  //

  std::vector<int> BezierTracker::DetermineNearbyHits(art::PtrVector<recob::Hit> Hits, BezierTrack* BTrack, double HitCollectionDistance)
  {
    std::vector<int> ReturnVector;
    double s, distance;
    std::cout<<"Size of hit vector : " <<Hits.size()<<std::endl;
    for(size_t i=0; i!=Hits.size(); ++i)
      {
	//	if((i%100)==0) std::cout<< i <<"  ";
	BTrack->GetClosestApproach(Hits.at(i),s, distance);
	if(distance<HitCollectionDistance) ReturnVector.push_back(i); 
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
	    //  float Angle = ThisSeed->GetAngle(*LastSeedAdded);
	    //  float ProjDis = ThisSeed->GetProjAngleDiscrepancy(*LastSeedAdded);
	    //  std::cout<<"BezierTracker: " << Angle<< " " <<ProjDis<<std::endl;
	    
            if((  abs(ThisSeed->GetAngle(*LastSeedAdded))                 < fMaxKinkAngle)
	       &&( abs(ThisSeed->GetProjAngleDiscrepancy(*LastSeedAdded))  < fMaxTrackMissAngle)
	       &&( abs(ThisSeed->GetDistance(*LastSeedAdded)<fMaxJumpDistance)));
                {
                  // if so, add it into the track, erase it from the stack
                  // and start looping at the beginning again
                  ThisTrack.push_back(ThisSeed);
                  TrackSeeds.erase(it);
		  it=TrackSeeds.begin();
                }
          }
        // We ran out of seeds on the stack. store this track and go again
	
	OrganizedByTrack.push_back(ThisTrack);
      }   
    return OrganizedByTrack;
  }


  void BezierTracker::endJob()
  {

  }
}
