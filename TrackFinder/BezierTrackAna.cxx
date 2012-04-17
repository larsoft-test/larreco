//
// Name: BezierTrackAna.cxx
//
// Purpose: Implementation file for module BezierTrackAna.
//
// Created: April 2012, bjpjones@mit.edu
//

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Utilities/DetectorProperties.h"
#include "TrackFinder/BezierTrackAna.h"
#include "TrackFinder/BezierTrack.h"
#include "Geometry/geo.h"
#include "Utilities/LArProperties.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "MCCheater/BackTracker.h"

#include "TTree.h"

namespace trkf {

  BezierTrackAna::BezierTrackAna(const fhicl::ParameterSet& pset) :
    //
    // Purpose: Constructor.
    //
    // Arguments: pset - Module parameters.
    //
    fBezierTrackModuleLabel(pset.get<std::string>("BezierTrackModuleLabel")),
    fBooked(false)
  {


  }

  BezierTrackAna::~BezierTrackAna()
  //
  // Purpose: Destructor.
  //
  {}

  void BezierTrackAna::bookHistograms(bool mc)
  //
  // Purpose: Book histograms.
  //
  {
    if(!fBooked) {
      fBooked = true;

      art::ServiceHandle<geo::Geometry> geom;
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory dir = tfs->mkdir("bezana", "BezierTrackAna histograms");

      //      if(mc) {
      //      }
     
      fHitS = tfs->make<TH1D>("HitS","HitS",100,0,1);
      fHitDistance = tfs->make<TH1D>("HitDistance","HitDistance",100,0,1);
      

      fTree          = (TTree*)tfs->make<TTree>("EventInfo","EventInfo");
      
      fTree->Branch("Length",        &fLength,       "Length/F");
      fTree->Branch("RMSCurvature",  &fRMSCurvature, "RMSCurvature/F");
      fTree->Branch("EventNumber",   &fNumEvent,     "EventNumber/I");
      fTree->Branch("NHits",         &fNHits,        "NHits/I");
      fTree->Branch("NHitsU",        &fNHitsU,       "NHitsU/I");
      fTree->Branch("NHitsV",        &fNHitsV,       "NHitsV/I");
      fTree->Branch("NHitsW",        &fNHitsW,       "NHitsW/I");
      fTree->Branch("AverageS",      &fAverageS,     "AverageS/I");
      fTree->Branch("AverageDistance",   &fAverageDistance,   "AverageDistance/I");
      fTree->Branch("dQdxU",         &fdQdxU,       "dQdxU/F");
      fTree->Branch("dQdxV",         &fdQdxV,       "dQdxV/F");
      fTree->Branch("dQdxW",         &fdQdxW,       "dQdxW/F");
      fTree->Branch("ChargeU",       &fChargeU,     "ChargeU/F");
      fTree->Branch("ChargeV",       &fChargeV,     "ChargeV/F");
      fTree->Branch("ChargeW",       &fChargeW,     "ChargeW/F");
      fTree->Branch("Segments",      &fSegments,    "Segments/I");
      
	
    }
    
  } 

  void BezierTrackAna::analyze(const art::Event& evt)
  //
  // Purpose: Analyze method.
  //
  // Arguments: event - Art event.
  //
  {
    ++fNumEvent;

    // Make sure histograms are booked.

    bool mc = !evt.isRealData();
    bookHistograms(mc);

    // Get Services.

    art::Handle< std::vector<recob::Track> > btbh;
    evt.getByLabel(fBezierTrackModuleLabel, btbh);   
    
    art::PtrVector<recob::Track> trackvec;   
    for(unsigned int i= 0; i<btbh->size(); ++i)
      {
	art::Ptr<recob::Track> trkptr(btbh,i);
	trackvec.push_back(trkptr);
      }
    
    
    
    
    std::vector<trkf::BezierTrack> BTracks;
    BTracks.clear();
    for(size_t i = 0; i < trackvec.size(); ++i){
      BTracks.push_back(trkf::BezierTrack( *(trackvec.at(i)) ));
    }
    
  
  if(BTracks.size()>0)
    for(size_t i=0; i!=BTracks.size(); ++i)
      {
	fLength       = BTracks.at(i).GetLength();
	fRMSCurvature = BTracks.at(i).GetRMSCurvature();
	

	art::PtrVector<recob::Hit> hits = util::FindManyP<recob::Hit>(trackvec, evt, fBezierTrackModuleLabel,i); 
	
	fNHitsU=fNHitsV=fNHitsW=fNHits=0;
	fAverageS = fAverageDistance = 0;
 
	for(art::PtrVector<recob::Hit>::const_iterator it=hits.begin();
	    it!=hits.end(); ++it)
	  {
	    fNHits++;
	    if((*it)->View()==geo::kU) fNHitsU++;
	    if((*it)->View()==geo::kV) fNHitsV++;
	    if((*it)->View()==geo::kW) fNHitsW++;
	    double HitDistance, S;
	    
	    BTracks.at(i).GetClosestApproach(*it, S, HitDistance);
	    
	    fHitS->Fill(S);
	    fHitDistance->Fill(HitDistance);
	    
	    fAverageS+=S;
	    fAverageDistance+=HitDistance;
	  }
	

	fAverageS        /= fNHits;
	fAverageDistance /= fNHits;
	
	fdQdxU = BTracks.at(i).GetViewdQdx(geo::kU);
  	fdQdxV = BTracks.at(i).GetViewdQdx(geo::kV);
  	fdQdxW = BTracks.at(i).GetViewdQdx(geo::kW);

  	fChargeU = BTracks.at(i).GetTotalCharge(geo::kU);
  	fChargeV = BTracks.at(i).GetTotalCharge(geo::kV);
  	fChargeW = BTracks.at(i).GetTotalCharge(geo::kW);
	
  	fSegments = BTracks.at(i).NSegments();
	
	
	
	fTree->Fill();
	std::cout<<"L : " <<fLength<<"  C : " <<fRMSCurvature<<std::endl;
	
	}
  }  
}
