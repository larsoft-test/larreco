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
#include "TrackFinder/SpacePointService.h"
#include "Geometry/geo.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "RecoBase/Hit.h"
#include "RecoBase/Seed.h"
#include "RecoBase/BezierTrackBase.h"
#include "RecoBase/Track.h"
#include "RecoBase/Prong.h"
#include "Utilities/AssociationUtil.h"

namespace trkf {

  BezierTracker::BezierTracker(const fhicl::ParameterSet& pset)
  {
    reconfigure(pset);
    produces<std::vector<recob::BezierTrackBase> >();
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

  }

  void BezierTracker::beginJob()
  {}


  void BezierTracker::produce(art::Event& evt)
  {
 
    //   int EventNumber = evt.id().event();

    art::Handle< std::vector<recob::Seed> > seedh;
    evt.getByLabel(fSeedModuleLabel, seedh);  
  
    art::Handle< std::vector<recob::Hit> > hith;
    evt.getByLabel(fHitModuleLabel, hith);
    
    std::auto_ptr<std::vector<recob::BezierTrackBase> > btracks(new std::vector<recob::BezierTrackBase>);

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
	recob::BezierTrackBase BezTrackBase = ProduceBaseTrack(OrgSeeds.at(i));
       	btracks->push_back(BezTrackBase);
      }
    evt.put(btracks);
    
  }

  recob::BezierTrackBase BezierTracker::ProduceBaseTrack(std::vector<art::Ptr<recob::Seed> > Seeds)
  {
    std::vector<double> PtX, PtY, PtZ, DirX, DirY, DirZ;
    double pt[3], dir[3], dummy[3];
    for(unsigned int i=0; i!=Seeds.size(); i++)
      {
	Seeds.at(i)->GetPoint(     pt,  dummy );
	Seeds.at(i)->GetDirection( dir, dummy );
	
	PtX.push_back(pt[0]);
	PtY.push_back(pt[1]);
	PtZ.push_back(pt[2]);

	DirX.push_back(dir[0]);
	DirY.push_back(dir[1]);
	DirZ.push_back(dir[2]);
      }
    return recob::BezierTrackBase(PtX, PtY, PtZ, DirX, DirY, DirZ, fTopTrackID++);
    
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
	    float Angle = ThisSeed->GetAngle(*LastSeedAdded);
	    float ProjDis = ThisSeed->GetProjAngleDiscrepancy(*LastSeedAdded);
	    std::cout<<"BezierTracker: " << Angle<< " " <<ProjDis<<std::endl;
	    
            if((  abs(ThisSeed->GetAngle(*LastSeedAdded))                 < fMaxKinkAngle)
	       &&( abs(ThisSeed->GetProjAngleDiscrepancy(*LastSeedAdded))  < fMaxTrackMissAngle))
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
