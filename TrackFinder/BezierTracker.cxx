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

  }

  void BezierTracker::beginJob()
  {}


  void BezierTracker::produce(art::Event& evt)
  {
 
    // Extract seeds from event
    art::Handle< std::vector<recob::Seed> > seedh;
    evt.getByLabel(fSeedModuleLabel, seedh);  
  

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
   
    
    std::vector<art::Ptr<recob::Seed> > TrackSeeds;
    if(seedh->size()>0)
      for(unsigned int iseed=0; iseed!=seedh->size(); iseed++)
	{
	  art::Ptr<recob::Seed> theseed(seedh, iseed);
	  TrackSeeds.push_back(theseed);
	}
    
    std::vector<std::vector<art::Ptr<recob::Seed> > > OrgSeeds = OrganizeSeedsIntoTracks(TrackSeeds);
   
    for(unsigned int i=0; i!=OrgSeeds.size(); i++)
      {
	std::cout<<"Seeds in this btrack : " << OrgSeeds.at(i).size();
	recob::Track BezTrackBase = ProduceBaseTrack(OrgSeeds.at(i));
	BezierTrack BTrack(BezTrackBase);
	
	std::vector<int> HitsToAssoc = DetermineNearbyHits(HitVec, &BTrack, fHitDistance);
	btracks->push_back(BezTrackBase);

	art::PtrVector<recob::Hit> HitAssocCol;
	
	for(size_t i=0; i!=HitsToAssoc.size(); i++)
	  {
	    HitAssocCol.push_back(HitVec.at(HitsToAssoc.at(i)));
	  }
	util::CreateAssn(*this, evt, *(btracks.get()), HitAssocCol, *(assn.get())); 
	
      }
    evt.put(btracks);
    evt.put(assn);
  }

  recob::Track BezierTracker::ProduceBaseTrack(std::vector<art::Ptr<recob::Seed> > Seeds)
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
    return recob::Track(Pos,Dir,dQdx);
    
  }


  //
  // From a PtrVector of hits, determine which are nearby
  //

  std::vector<int> BezierTracker::DetermineNearbyHits(art::PtrVector<recob::Hit> Hits, BezierTrack* BTrack, double HitCollectionDistance)
  {
    std::vector<int> ReturnVector;
    double s, distance;
    for(size_t i=0; i!=Hits.size(); i++)
      {
	BTrack->GetClosestApproach(Hits.at(i),s, distance);
	std::cout<<"Hit found with distance " <<distance;
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
