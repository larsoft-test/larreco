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
     
      fRMSCurvatures = dir.make<TH1D>("RMSCurvature","RMSCurvature",100,0,0.01);
      fLengths       = dir.make<TH1D>("Length",      "Length",      100,0,0.01);
      fTree          = (TTree*)tfs->make<TTree>("EventInfo","EventInfo");
      
      fTree->Branch("Length",        &fLength,       "Length/F");
      fTree->Branch("RMSCurvature",  &fRMSCurvature, "RMSCurvature/F");
      fTree->Branch("EventNumber",   &fNumEvent,  "EventNumber/I");
      
	
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
    std::vector<trkf::BezierTrack> BTracks;
    BTracks.clear();
    for(size_t i = 0; i < btbh->size(); ++i){
      art::Ptr<recob::Track> btb(btbh, i);
      BTracks.push_back(trkf::BezierTrack(*btb));
    }

  
  if(BTracks.size()>0)
    for(size_t i=0; i!=BTracks.size(); ++i)
      {
	fLength       = BTracks.at(i).GetLength();
	fRMSCurvature = BTracks.at(i).GetRMSCurvature();
	std::cout<<"Getting RMS Curvature"<<std::endl;
	fRMSCurvatures->Fill(fRMSCurvature);
	std::cout<<"Getting Length"<<std::endl;
	fLengths->Fill(fLength);
	fTree->Fill();
	std::cout<<"L : " <<fLength<<"  C : " <<fRMSCurvature<<std::endl;
      }
  }  
}
