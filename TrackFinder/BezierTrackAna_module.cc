#ifndef BEZIERTRACKANA_H
#define BEZIERTRACKANA_H

//
// Name: BezierTrackAna.h
//
// Purpose: Header file for module BezierTrackAna.
//
// Created: April 2012, bjpjones@mit.edu
//

#include "art/Framework/Core/EDAnalyzer.h"
#include "TH1F.h"
#include "TH2D.h"
class TTree;

namespace trkf {

  class BezierTrackAna : public art::EDAnalyzer
  {
  public:
 
    // Constructors, destructor

    explicit BezierTrackAna(fhicl::ParameterSet const& pset);
    virtual ~BezierTrackAna();

    // Book histograms.

    void bookHistograms(bool mc);

    // Overrides.

    void analyze(const art::Event& evt);

  private:

    // Fcl Attributes.

    std::string fBezierTrackModuleLabel;
    bool fBooked;

    // Histograms
    
    //TH1D * fRMSCurvatures;
    //TH1D * fLengths;
    
    TH1D * fhHitDistance;
    TH1D * fhHitS;
    TH1D * fhdQdxU;
    TH1D * fhdQdxV;
    TH1D * fhdQdxW;
    TH2D * fhdQdxVW;
    TH1D * fhCurv;
    
    
   
    //Entries for TTree
    TTree*  fTree;
    Float_t fLength;
    Float_t fRMSCurvature;
    Int_t   fNumEvent;
    Int_t   fNHitsU;
    Int_t   fNHitsV;
    Int_t   fNHitsW;
    Int_t   fNHits;
    Float_t fAverageS;
    Float_t fAverageDistance;
    
    Float_t fdQdxU;
    Float_t fdQdxV;
    Float_t fdQdxW;
    
    Float_t fChargeU;
    Float_t fChargeV;
    Float_t fChargeW;
    
    Int_t fSegments;


  };
}

#endif 



#include "art/Framework/Core/ModuleMacros.h" 

namespace trkf {
  DEFINE_ART_MODULE(BezierTrackAna);
}


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
#include "TrackFinder/BezierTrack.h"
#include "Geometry/Geometry.h"
#include "Utilities/LArProperties.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"

#include "TTree.h"

#include "TH2D.h"
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
     
      fhHitS = tfs->make<TH1D>("HitS","HitS",100,0,1);
      fhHitDistance = tfs->make<TH1D>("HitDistance","HitDistance",100, 0.0, 1.0);
      fhdQdxU = tfs->make<TH1D>("dQdxU", "dQdxU", 100, 0.0, 1.0);
      fhdQdxV = tfs->make<TH1D>("dQdxV", "dQdxV", 100, 0.0, 1.0);
      fhdQdxW = tfs->make<TH1D>("dQdxW", "dQdxW", 100, 0.0, 1.0);
      
      fhdQdxVW =  tfs->make<TH2D>("dQdxVW", "dQdxWV", 100, 600, 1200,100,600,1200);
      fhCurv  = tfs->make<TH1D>("Curv",  "Curv",  1000, 0.0, 1.0);

      fTree          = (TTree*)tfs->make<TTree>("EventInfo","EventInfo");
      
      fTree->Branch("Length",        &fLength,       "Length/F");
      fTree->Branch("RMSCurvature",  &fRMSCurvature, "RMSCurvature/F");
      fTree->Branch("EventNumber",   &fNumEvent,     "EventNumber/I");
      fTree->Branch("NHits",         &fNHits,        "NHits/I");
      fTree->Branch("NHitsU",        &fNHitsU,       "NHitsU/I");
      fTree->Branch("NHitsV",        &fNHitsV,       "NHitsV/I");
      fTree->Branch("NHitsW",        &fNHitsW,       "NHitsW/I");
      fTree->Branch("AverageS",      &fAverageS,     "AverageS/F");
      fTree->Branch("AverageDistance",   &fAverageDistance,   "AverageDistance/F");
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
    
    art::FindManyP<recob::Hit> fmh(btbh, evt, fBezierTrackModuleLabel);
    
    
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
	

	std::vector< art::Ptr<recob::Hit> > hits = fmh.at(i); 
	
	fNHitsU=fNHitsV=fNHitsW=fNHits=0;
	fAverageS = fAverageDistance = 0;
	
	std::vector<double> SVals; std::vector<double> DVals;
	BTracks.at(i).GetClosestApproaches(hits, SVals, DVals);

	for(size_t j=0; j!=SVals.size();++j)
	  {
	    fhHitS->Fill(SVals.at(j));
	    fhHitDistance->Fill(DVals.at(j));
	    
	    fAverageS+=SVals.at(j);
	    fAverageDistance+=DVals.at(j);

	  }
	/*
	for(art::PtrVector<recob::Hit>::const_iterator it=hits.begin();
	    it!=hits.end(); ++it)
	  {
	    fNHits++;
	    if((*it)->View()==geo::kU) fNHitsU++;
	    if((*it)->View()==geo::kV) fNHitsV++;
	    if((*it)->View()==geo::kZ) fNHitsW++;
	    double HitDistance, S;
	    
	    
	   
	  }
	*/
	int jDivs = 100;
	for(int j=1; j!=jDivs-1; j++)
	  {
	    float Point = float(j)/jDivs + 0.001;
	    std::cout<<"btrkana" << Point<< " " << BTracks.at(i).GetCurvature(float(j)/jDivs)<<std::endl;
	    //	    fhdQdxU->Fill( Point, BTracks.at(i).GetdQdx(float(j)/jDivs, geo::kU));
	    //	    fhdQdxV->Fill( Point, BTracks.at(i).GetdQdx(float(j)/jDivs, geo::kV));
	    //	    fhdQdxW->Fill( Point, BTracks.at(i).GetdQdx(float(j)/jDivs, geo::kZ));
	    //	    fhdQdxVW->Fill( BTracks.at(i).GetdQdx(float(j)/jDivs, geo::kZ), BTracks.at(i).GetdQdx(float(j)/jDivs, geo::kV));
	    fhCurv->Fill(  Point, BTracks.at(i).GetCurvature(float(j)/jDivs));
			 
	  }
	
	//	fAverageS        /= fNHits;
	//fAverageDistance /= fNHits;
	
	//	fdQdxU = BTracks.at(i).GetViewdQdx(geo::kU);
	// 	fdQdxV = BTracks.at(i).GetViewdQdx(geo::kV);
	// 	fdQdxW = BTracks.at(i).GetViewdQdx(geo::kZ);

	//  	fChargeU = BTracks.at(i).GetTotalCharge(geo::kU);
	//	fChargeV = BTracks.at(i).GetTotalCharge(geo::kV);
	//	fChargeW = BTracks.at(i).GetTotalCharge(geo::kZ);
	
  	fSegments = BTracks.at(i).NSegments();
	
	
	
	fTree->Fill();
	std::cout<<"L : " <<fLength<<"  C : " <<fRMSCurvature<<std::endl;
	
	}
  }  
}
