////////////////////////////////////////////////////////////////////////
// $Id$
//
// MagDriftAna class designed to make histograms
//
// brebel@fnal.gov
//
//
////////////////////////////////////////////////////////////////////////
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

// LArSoft includes
#include "TrackFinder/MagDriftAna.h"
#include "Geometry/geo.h"
#include "MagneticField/mag.h"
#include "Utilities/LArProperties.h"
#include "MCCheater/BackTracker.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"

// ROOT includes
#include "TMath.h"
#include "TH2.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLine.h"

// C++ includes
#include <algorithm>
#include <sstream>
#include <fstream>
#include <bitset>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace hit{

  //-------------------------------------------------
  MagDriftAna::MagDriftAna(fhicl::ParameterSet const& pset)
    :initDone(false)
    ,fDirCosY(0.0),fDirCosZ(0.0)
    ,fChargeXpos(),fChargeYpos(),fChargeZpos()
    ,fHitZpos()
    ,fDriftDeltaZ(),fDeltaZoverX(),fDeltaZvsX()
    ,fDriftDeltaZAway(),fDeltaZoverXAway()
  {
    this->reconfigure(pset);
  }

  //-------------------------------------------------
  MagDriftAna::~MagDriftAna()
  {
  }

  void MagDriftAna::reconfigure(fhicl::ParameterSet const& p)
  {
    fFFTHitFinderModuleLabel  = p.get< std::string >("HitsModuleLabel");
//     fTrackFinderModuleLabel   = p.get< std::string >("TracksModuleLabel");
    fLArG4ModuleLabel         = p.get< std::string >("LArGeantModuleLabel");
    return;
  }
  //-------------------------------------------------
  void MagDriftAna::ensureHists() {
    if (initDone) return; // Bail if we've already done this.
    initDone = true; // Insure that we bail later on

    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    // Find magetic field related corrections
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<mag::MagneticField> MagField;
    double Efield         = larprop->Efield();  
    double Temperature    = larprop->Temperature();  
    double DriftVelocity  = larprop->DriftVelocity(Efield,Temperature)/1000.;
    
    if (MagField->UseField()) {
      fDirCosY = -DriftVelocity * MagField->FieldAtPoint().z() / Efield;
      fDirCosZ = +DriftVelocity * MagField->FieldAtPoint().y() / Efield;
      std::cerr << "Drift ratios: "
		<< "dY/dX = " << fDirCosY << ", " 
		<< "dZ/dX = " << fDirCosZ 
		<< std::endl;
    } else {
      std::cerr << "Why are you running magdriftana without a magnetic field?"
		<< std::endl;
      std::cerr << "\t" <<  "Results are going to be nonesense."
		<< std::endl;
    }
      
    // geometry data.
    art::ServiceHandle<geo::Geometry> geom; 
    // assumes all TPCs are the same
    double width  = 2 * geom->TPC(0).HalfWidth(); 
    double halfHeight = geom->TPC(0).HalfHeight(); 
    double length = geom->TPC(0).Length(); 

    double zScale = std::max(fDirCosZ,2e-4);

    // Assumes microboone dimensions. Ideally we'd fix this later...
    fChargeXpos  = tfs->make<TH1D>("hChargeXpos",
				   "MC X charge depositions; X (cm); Events",
				   101, 0.0, width);
    fChargeYpos  = tfs->make<TH1D>("hChargeYpos",
				   "MC Y charge depositions; Y (cm); Events",
				   101, -halfHeight, halfHeight);
    fChargeZpos  = tfs->make<TH1D>("hChargeZpos",
				   "MC Z charge depositions; Z (cm); Events",
				   101, 0.0, length);
    fHitZpos     = tfs->make<TH1D>("hHitZpos",
				   "Z charge collection; Z (cm); Events",
				   101, 0.0, length);

    fDriftDeltaZ = tfs->make<TH1D>("hDriftDeltaZ",
				   "Z drift of charge; delta Z (cm); Events",
				   101, -5*zScale*width, 5*zScale*width);
    fDeltaZoverX = tfs->make<TH1D>("hDeltaZoverX",
				   "Z drift of charge; delta Z/X; Events",
				   101, -20*zScale, 20*zScale);
    fDeltaZvsX = tfs->make<TH2D>("hDeltaZvsX",
				 "delta Z vs X; X (cm); delta Z (cm), Events",
				 51, 0.0, width,
				 51, -20*zScale, 20*zScale);

    // Some stats only when the Xdrift is large (more than 3/4)
    fDriftDeltaZAway = 
      tfs->make<TH1D>("hDriftDeltaZAway",
		      "Z drift of charge (long drift); delta Z (cm); Events",
		      101, -5*zScale*width, 5*zScale*width);
    fDeltaZoverXAway = 
      tfs->make<TH1D>("hDeltaZoverXAway",
		      "Z drift of charge (long drift); delta Z/X; Events",
		      101, -20*zScale, 20*zScale);

    return;

  }

  //-------------------------------------------------
  void MagDriftAna::endJob() 
  {
    // Add a line on the deltaZ/X graph to denote the calculated value
    // of the drift ration
    TLine *l = new TLine(fDirCosZ, 0,
			 fDirCosZ, 1.05*fDeltaZoverX->GetMaximum());
    l->SetLineColor(kRed);
    l->SetLineStyle(kDotted);
    fDeltaZoverX->GetListOfFunctions()->Add(l);

    // I know this looks like a memory leak, but each historgram needs
    // it's own copy of the line to prevent double freeing by the
    // framework...
    ////////////
    l = new TLine(fDirCosZ, 0,
		  fDirCosZ, 1.05*fDeltaZoverX->GetMaximum());
    ///////////
    //
    l->SetLineColor(kRed);
    l->SetLineStyle(kDotted);
    fDeltaZoverX->GetListOfFunctions()->Add(l);
    fDeltaZoverXAway->GetListOfFunctions()->Add(l);
  }

  //-------------------------------------------------
  void MagDriftAna::analyze(const art::Event& evt)
  {

    if (evt.isRealData()) {
      throw cet::exception("MagDriftAna: ") << "Not for use on Data yet... " 
					    << "\n";
    }
    
    ensureHists();

    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fFFTHitFinderModuleLabel,hitHandle);

//     art::Handle< std::vector<recob::Track> > trackHandle;
//     evt.getByLabel(fTrackFinderModuleLabel,trackHandle);

    art::ServiceHandle<geo::Geometry> geom; 

    // We're going to want to compare the reconstructed Z with the
    // simulted Z. For that purpose we use the simultion backtracking.
    //

    art::ServiceHandle<cheat::BackTracker> bt;

    //    art::PtrVector<recob::Hit> hits;
    std::vector< art::Ptr<recob::Hit> > hits;
    art::fill_ptr_vector(hits, hitHandle);
//     std::vector< art::Ptr<recob::Track> > tracks;
//     art::fill_ptr_vector(tracks, trackHandle);
    
    unsigned int p(0); // plane
    unsigned int w(0); // wire
    unsigned int t(0); // TPC 
    unsigned int c(0); // cryostat
    unsigned int channel(0);

    //++++++++++
    // Loop over the hits (data) and fill histos
    //++++++++++
    for ( std::vector< art::Ptr<recob::Hit> >::iterator itr = hits.begin();
	  itr != hits.end();
	  ++itr) {
      // Find the DAQ channel associated with the hit
      channel=(*itr)->Wire()->RawDigit()->Channel();
      c=0;t=0;p=0;w=0;
      geom->ChannelToWire(channel,c,t,p,w);
      // By assumption the drift occurs only in the z-direction, so
      // we can get all the info we need from the z-measug plane.
      if (p != (geom->Nplanes()-1)) continue;
      
      // Charge collected at the wire 
      //
      // Exactly once for each recob::Hit
      double w0pos[3] = {0.};
      geom->TPC(t).Plane(p).Wire(0).GetCenter(w0pos);
      double HitZpos = w0pos[2] + w * geom->TPC(t).WirePitch();
      double Charge = (*itr)->Charge();
      fHitZpos->Fill(HitZpos,Charge);
	
      // Charge deposition in the detector
      std::vector<double> xyz = bt->HitToXYZ(*itr);
      fChargeXpos->Fill(xyz[0],Charge);
      fChargeYpos->Fill(xyz[1],Charge);
      double ChargeZpos = xyz[2];
      fChargeZpos->Fill(ChargeZpos,Charge);

// 	// Delta-Y
// 	fDriftDeltaY->Fill(HitYpos-ChargeYpos,Charge);
// 	// Delta Y correlation with X
// 	fDeltaYoverX->Fill((HitYpos-ChargeYpos)/xyz[0],Charge);
// 	fDeltaYvsX->Fill(xyz[0],HitYpos-ChargeZpos,Charge);

      // Delta-Z
      //
      // Compares the collected z position from the wire to the
      // simulated z position
      double DeltaZ = HitZpos-ChargeZpos;
      fDriftDeltaZ->Fill(DeltaZ,Charge);
      // Delta Z correlation with X
      fDeltaZoverX->Fill(DeltaZ/xyz[0],Charge);
      fDeltaZvsX->Fill(xyz[0],DeltaZ,Charge);
      // The X related histograms know the dimensions of the
      // detector, so we use them to set the "away" limit
      if (xyz[0] >  (fChargeYpos->GetXaxis()->GetXmax() * 0.80) ) {  
	fDriftDeltaZAway->Fill(DeltaZ,Charge);
	fDeltaZoverXAway->Fill(DeltaZ/xyz[0],Charge);
      } 

    } // loop on Hits
    

    return;
  }//end analyze method

}//end namespace
