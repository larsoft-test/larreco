////////////////////////////////////////////////////////////////////////
//
// GausHitFinderAna class designed to make histograms
//
// jaasaadi@syr.edu
//
//  This algorithm is designed to analyze hits created on wires after 
//  deconvolution and can be used with FFTHitFinder or GausHitFinder
//
// Note: This has been based (stolen) from the FFTHitFinderAna thus
// there is still some unneeded hold overs that will get cleaned up later
//
// To use this simply include in your analyzers:
// gaushitfinderana: @local::gaus_hitfinderana
////////////////////////////////////////////////////////////////////////
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

// LArSoft includes
#include "HitFinder/GausHitFinderAna.h"
#include "Geometry/geo.h"
#include "MCCheater/BackTracker.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"

// ROOT includes
#include <TMath.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TFile.h>

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
  GausHitFinderAna::GausHitFinderAna(fhicl::ParameterSet const& pset) 
  {
    this->reconfigure(pset);
  }

  //-------------------------------------------------
  GausHitFinderAna::~GausHitFinderAna()
  {
  }

  void GausHitFinderAna::reconfigure(fhicl::ParameterSet const& p)
  {
    fGausHitFinderModuleLabel = p.get< std::string >("HitsModuleLabel");
    fLArG4ModuleLabel        = p.get< std::string >("LArGeantModuleLabel");
    return;
  }
  //-------------------------------------------------
  void GausHitFinderAna::beginJob() 
  {
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
    /*fNp0 = 9000;
    fNp1 = 9000;
    fNp2 = 9000;*/

    fHTree = tfs->make<TTree>("HTree","HTree");
    // ------------------------------------
    // --- Defining Peak Time Variables ---
    // ------------------------------------
    /*fPeakTime0 = new Float_t[9000];
    fPeakTime1 = new Float_t[9000];
    fPeakTime2 = new Float_t[9000];
    
    
    fWirep0 = new Int_t[fNp0];
    fWirep1 = new Int_t[fNp1];
    fWirep2 = new Int_t[fNp2];
    fChgp0 = new Float_t[fNp0];
    fChgp1 = new Float_t[fNp1];
    fChgp2 = new Float_t[fNp2];
    fXYZp0 = new Float_t[fNp0*3];
    fXYZp1 = new Float_t[fNp1*3];
    fXYZp2 = new Float_t[fNp2*3];

    fMCPdg0 = new Int_t[fNp0];
    fMCPdg1 = new Int_t[fNp1];
    fMCPdg2 = new Int_t[fNp2];
    fMCTId0 = new Int_t[fNp0];
    fMCTId1 = new Int_t[fNp1];
    fMCTId2 = new Int_t[fNp2];
    fMCE0 = new Float_t[fNp0];
    fMCE1 = new Float_t[fNp1];
    fMCE2 = new Float_t[fNp2];*/

    fHTree->Branch("HEvt", &fEvt, "HEvt/I");
    fHTree->Branch("HRun", &fRun, "HRun/I");
    fHTree->Branch("NHits", &fnhits, "NHits/I");
    fHTree->Branch("SinglePulseEvent", &fSingleHit, "SingleHitEvent/I");
    fHTree->Branch("MulitPulseEvent", &fMultiHit, "MulitPulseEvent/I");
    
    
    fHTree->Branch("WireNumbern1", &fWiren1, "WireNumbern1/I");
    fHTree->Branch("NOnePulseHit", &fnOnePulseHits, "NOnePulseHit/I");
    fHTree->Branch("GOFMulti1", &fgoodoffitn1, "GOFMulti1/F");
    fHTree->Branch("ChargeMulti1", &fChargen1, "ChargeMulti1/F");
    fHTree->Branch("SigmaChargeMulti1", &fSigmaChargen1, "SigmaChargeMulti1/F");
    fHTree->Branch("WidthMulti1", &fWidthn1, "WidthMulti1/F");
    fHTree->Branch("PeakPosMulti1", &fPeakn1, "PeakPosMulti1/F");
    fHTree->Branch("PeakUncertMulti1", &fPeakUncertn1, "PeakUncertMulti1/F");
    fHTree->Branch("StartPosMulti1", &fStartTimen1, "StartPosMulti1/F");
    fHTree->Branch("StartPosUncertMulti1", &fStartTimeUncertn1, "StartPosUncertMulti1/F");
    fHTree->Branch("EndPosMulti1", &fEndTimen1, "EndPosMulti1/F");
    fHTree->Branch("EndPosUncertMulti1", &fEndTimeUncertn1, "EndPosUncertMulti1/F");
    
    
    fHTree->Branch("WireNumbernGT1", &fWirenGT1, "WireNumbernGT1/I");
    fHTree->Branch("NMultiPulseHit", &fmulitPulseHits, "NMultiPulseHit/I");
    fHTree->Branch("GOFMultiGT1", &fgoodoffitnGT1, "GOFMultiGT1/F");
    fHTree->Branch("ChargeMultiGT1", &fChargenGT1, "ChargeMultiGT1/F");
    fHTree->Branch("SigmaChargeMultiGT1", &fSigmaChargenGT1, "SigmaChargeMultiGT1/F");
    fHTree->Branch("WidthMultiGT1", &fWidthnGT1, "WidthMultiGT1/F");
    fHTree->Branch("PeakPosMultiGT1", &fPeaknGT1, "PeakPosMultiGT1/F");
    fHTree->Branch("PeakUncertMultiGT1", &fPeakUncertnGT1, "PeakUncertMultiGT1/F");
    fHTree->Branch("StartPosMultiGT1", &fStartTimenGT1, "StartPosMultiGT1/F");
    fHTree->Branch("StartPosUncertMultiGT1", &fStartTimeUncertnGT1, "StartPosUncertMultiGT1/F");
    fHTree->Branch("EndPosMultiGT1", &fEndTimenGT1, "EndPosMultiGT1/F");
    fHTree->Branch("EndPosUncertMultiGT1", &fEndTimeUncertnGT1, "EndPosUncertMultiGT1/F");
    
    /*fHTree->Branch("Hit0PeakTime", &fNp0, "Hit0PeakTime/I");
    fHTree->Branch("HNp1", &fNp1, "HNp1/I");
    fHTree->Branch("HNp2", &fNp2, "HNp2/I");
    fHTree->Branch("HN3p0", &fN3p0, "HN3p0/I");
    fHTree->Branch("HN3p1", &fN3p1, "HN3p1/I");
    fHTree->Branch("HN3p2", &fN3p2, "HN3p2/I");
    fHTree->Branch("Htp0", fPeakTime0, "Htp0[Hit0PeakTime]/F");
    fHTree->Branch("Htp1", fPeakTime1, "Htp1[HNp1]/F");
    fHTree->Branch("Htp2", fPeakTime2, "Htp2[HNp2]/F");
    fHTree->Branch("Hwp0", fWirep0, "Hwp0[Hit0PeakTime]/I");
    fHTree->Branch("Hwp1", fWirep1, "Hwp1[HNp1]/I");
    fHTree->Branch("Hwp2", fWirep2, "Hwp2[HNp2]/I");
    fHTree->Branch("Hchgp0", fChgp0, "Hchgp0[Hit0PeakTime]/F");
    fHTree->Branch("Hchgp1", fChgp1, "Hchgp1[HNp1]/F");
    fHTree->Branch("Hchgp2", fChgp2, "Hchgp2[HNp2]/F");
    fHTree->Branch("HMCXYZp0", fXYZp0, "HMCXYZp0[HN3p0]/F");
    fHTree->Branch("HMCXYZp1", fXYZp1, "HMCXYZp1[HN3p1]/F");
    fHTree->Branch("HMCXYZp2", fXYZp2, "HMCXYZp2[HN3p2]/F");
    fHTree->Branch("HMCPdgp0", fMCPdg0, "HMCPdgp0[Hit0PeakTime]/I");
    fHTree->Branch("HMCPdgp1", fMCPdg1, "HMCPdgp1[HNp1]/I");
    fHTree->Branch("HMCPdgp2", fMCPdg2, "HMCPdgp2[HNp2]/I");
    fHTree->Branch("HMCTIdp0", fMCTId0, "HMCTIdp0[Hit0PeakTime]/I");
    fHTree->Branch("HMCTIdp1", fMCTId1, "HMCTIdp1[HNp1]/I");
    fHTree->Branch("HMCTIdp2", fMCTId2, "HMCTIdp2[HNp2]/I");
    fHTree->Branch("HMCEp0", fMCE0, "HMCEp0[Hit0PeakTime]/F");
    fHTree->Branch("HMCEp1", fMCE1, "HMCEp1[HNp1]/F");
    fHTree->Branch("HMCEp2", fMCE2, "HMCEp2[HNp2]/F");*/

  
    return;

  }

  //-------------------------------------------------
  void GausHitFinderAna::analyze(const art::Event& evt)
  {

  // ##############################################
  // ### Outputting Run Number and Event Number ###
  // ##############################################
  //std::cout << "run    : " << evt.run() <<" event  : "<<evt.id().event() << std::endl;
  int NSinglePulseEvents = 0 , NMultiPulseEvents = 0;
  
  int SinglePulse = 0, Multipulse = 0;
  fRun = evt.run();
  fEvt = evt.id().event();
  
  // ####################################
  // ### Getting Geometry Information ###
  // ####################################
  art::ServiceHandle<geo::Geometry> geom; 

  // ##################################################
  // ### Getting the Reconstructed Hits (hitHandle) ###
  // ##################################################
  art::Handle< std::vector<recob::Hit> > hitHandle;
  evt.getByLabel(fGausHitFinderModuleLabel,hitHandle);
  //art::PtrVector<recob::Hit> allhits;
  
  // #########################################
  // ### Putting Hits into a vector (hits) ###
  // #########################################
  std::vector< art::Ptr<recob::Hit> > hits;
  art::fill_ptr_vector(hits, hitHandle);
  
  unsigned int channel = 0, c = 0, t = 0, p = 0, w = 0;

  std::cout<<std::endl;
  std::cout<<"Number of Hits in the Event = "<<hitHandle->size()<<std::endl;
  fnhits = hitHandle->size();

  // #########################
  // ### Looping over Hits ###
  // #########################
  for(int nHits = 0; nHits< hitHandle->size(); nHits++)
  	{
	//std::cout<<"Hit = "<<nHits<<std::endl;
	
	// === Finding Channel associated with the hit ===
	art::Ptr<recob::Hit> hit(hitHandle, nHits);
	channel= hit->Wire()->RawDigit()->Channel();
	//std::cout<<"channel = "<<channel<<std::endl;
	// === Going from the Channel to the wire location ===
	// (Note:3/16/12 Channel to wire function now reads
	// (Channedl,cryostat,tpc,plane,wire)
	geom->ChannelToWire(channel,c,t,p,w);
	
	
	// ##################################################
	// ### Looking at "Hits" with a multiplicity == 1 ###
	// ##################################################
	if(hit->Multiplicity() == 1)
		{
		NSinglePulseEvents++;
		SinglePulse = 1;
		
		fSingleHit         = SinglePulse;
		fWiren1            = w;
		fgoodoffitn1       = hit->GoodnessOfFit();
		fChargen1          = hit->Charge();
		fSigmaChargen1     = hit->SigmaCharge();
		fWidthn1           = (hit->EndTime() - hit->PeakTime());
		
		fPeakn1            = hit->PeakTime();	
		fPeakUncertn1      = hit->SigmaPeakTime();
		fStartTimen1       = hit->StartTime();
		fStartTimeUncertn1 = hit->SigmaStartTime();
		fEndTimen1         = hit->EndTime();
		fEndTimeUncertn1   = hit->SigmaEndTime();
		
		}//<---End Hit Multiplicity == 1
	
	// ##################################################
	// ### Looking at "Hits" with a multiplicity == 1 ###
	// ##################################################
	if(hit->Multiplicity() > 1)
		{
		Multipulse = 1;
		NMultiPulseEvents++;
		
		fMultiHit            = Multipulse;
		fWirenGT1            = w;
		fgoodoffitnGT1       = hit->GoodnessOfFit();
		fChargenGT1          = hit->Charge();
		fSigmaChargenGT1     = hit->SigmaCharge();
		fWidthnGT1           = (hit->EndTime() - hit->PeakTime());
		
		fPeaknGT1            = hit->PeakTime();
		fPeakUncertnGT1      = hit->SigmaPeakTime();
		fStartTimenGT1       = hit->StartTime();
		fStartTimeUncertnGT1 = hit->SigmaStartTime();
		fEndTimenGT1         = hit->EndTime();
		fEndTimeUncertnGT1   = hit->SigmaEndTime();
		
		
		}//<---End Hit Multiplicity > 1
		
	
	/*std::cout<<"c = "<<c<<" t = "<<t<<" p = "<<p<<" w = "<<w<<std::endl;
	std::cout<<"Start Time        = "<<	hit->StartTime()	<<	std::endl;
	std::cout<<"Sigma Start Time  = "<<	hit->SigmaStartTime()	<<	std::endl;
	std::cout<<"End Time          = "<<	hit->EndTime()		<<	std::endl;
	std::cout<<"Sigma End Time    = "<<	hit->SigmaEndTime()	<<	std::endl;
	std::cout<<"Peak Time         = "<<	hit->PeakTime()		<<	std::endl;
	std::cout<<"Sigma Peak Time   = "<<	hit->SigmaPeakTime()	<<	std::endl;
	std::cout<<"Multiplicity      = "<<	hit->Multiplicity()	<<	std::endl;
	std::cout<<"Charge            = "<<	hit->Charge()		<<	std::endl;
	std::cout<<"Sigma Charge      = "<<	hit->SigmaCharge()	<<	std::endl;
	std::cout<<"Goodness Fit      = "<<	hit->GoodnessOfFit()	<<	std::endl;
	std::cout<<std::endl;*/
	
	
	
	fHTree->Fill();
	Multipulse = 0;
	SinglePulse = 0;
  	}//<---End Loop over hits
    fnOnePulseHits = NSinglePulseEvents;
    fmulitPulseHits = NMultiPulseEvents;
    fHTree->Fill();
    return;
    
  }//end analyze method
  
}//end namespace

/*
    if (evt.isRealData()){
      throw cet::exception("HitFinderAna: ") << "Not for use on Data yet... " << "\n";
    }
    
    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fFFTHitFinderModuleLabel,hitHandle);

    sim::ParticleList _particleList = sim::SimListUtils::GetParticleList(evt, fLArG4ModuleLabel);
    std::vector<const sim::SimChannel*> sccol;
    evt.getView(fLArG4ModuleLabel, sccol);

    std::cout << _particleList << std::endl;

    //    art::PtrVector<recob::Hit> hits;
    std::vector< art::Ptr<recob::Hit> > hits;
    art::fill_ptr_vector(hits, hitHandle);
    
    art::ServiceHandle<geo::Geometry> geom;  
  
    unsigned int p(0),w(0), t(0), cs(0), channel(0);
    for(unsigned int cstat = 0; cstat < geom->Ncryostats(); ++cstat){
      for(unsigned int tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc){
	//      for(unsigned int plane=0; plane<geom->Nplanes(tpc); plane++){
	//	for(unsigned int i = 0; i< hitcol->size(); ++i){
	fNp0=0;       fN3p0=0;
	fNp1=0;       fN3p1=0;
	fNp2=0;       fN3p2=0;
	
	//now make a vector where each channel in the detector is an entry
	std::vector<const sim::SimChannel*> scs(geom->Nchannels(),0);
	for(size_t i = 0; i < sccol.size(); ++i) scs[sccol[i]->Channel()] = sccol[i];
	
	std::vector< art::Ptr<recob::Hit> >::iterator itr = hits.begin();
	while(itr != hits.end()) {
	  
	  //art::Ptr<recob::Hit> hit(hitHandle, i);
	  channel=(*itr)->Wire()->RawDigit()->Channel();
	  cs=0;t=0;p=0;w=0;
	  geom->ChannelToWire(channel, cs, t, p, w);
	  
	  fRun = evt.run();
	  fEvt = evt.id().event();
	  

	  if (!scs[channel]) {itr++;continue;}
	  
	  std::vector<cheat::TrackIDE> trackides = cheat::BackTracker::HitToTrackID(*(scs[(*itr)->Channel()]), *itr);
	  std::vector<cheat::TrackIDE>::iterator idesitr = trackides.begin();
	  std::vector<double> xyz = cheat::BackTracker::HitToXYZ(*(scs[(*itr)->Channel()]),*itr);
	  
	  
	  if (p==0 && fNp0<9000){
	    fPeakTime0[fNp0] = (*itr)->PeakTime();
	    fWirep0[fNp0] = w;
	    fChgp0[fNp0] = (*itr)->Charge();
	    
	    for (unsigned int kk=0;kk<3;kk++){
	      fXYZp0[fNp0*3+kk] = xyz[kk];
	    }
	    
	    
	    while( idesitr != trackides.end() ){
	      fMCTId0[fNp0] = (*idesitr).trackID;
	      if (_particleList.find((*idesitr).trackID) != _particleList.end()){
		const sim::Particle* particle = _particleList.at( (*idesitr).trackID);
		fMCPdg0[fNp0] = particle->PdgCode();
		fMCE0[fNp0] = particle->E();
	      }
	      idesitr++;
	    }
	    
	    fNp0++;
	  }
	  
	  else if (p==1 && fNp1<9000){
	    fPeakTime1[fNp1] = (*itr)->PeakTime();
	    fWirep1[fNp1] = w;
	    fChgp1[fNp1] = (*itr)->Charge();
	    
	    for (unsigned int kk=0;kk<3;kk++){
	      fXYZp1[fNp1*3+kk] = xyz[kk];
	    }
	    
	    while( idesitr != trackides.end() ){
	      fMCTId1[fNp1] = (*idesitr).trackID;
	      if (_particleList.find((*idesitr).trackID) != _particleList.end()){
		const sim::Particle* particle = _particleList.at( (*idesitr).trackID);
		fMCPdg1[fNp1] = particle->PdgCode();
		fMCE1[fNp1] = particle->E();
	      }
	      idesitr++;
	    }
	    fNp1++;
	  }
	  
	  else if (p==2  && fNp2<9000){
	    fPeakTime2[fNp2] = (*itr)->PeakTime();
	    fWirep2[fNp2] = w;
	    fChgp2[fNp2] = (*itr)->Charge();
	    
	    for (unsigned int kk=0;kk<3;kk++){
	      fXYZp2[fNp2*3+kk] = xyz[kk];
	    }
	    
	    while( idesitr != trackides.end()){
	      fMCTId2[fNp2] = (*idesitr).trackID;
	      if (_particleList.find((*idesitr).trackID) != _particleList.end() ){
		const sim::Particle* particle = _particleList.at( (*idesitr).trackID);
		fMCPdg2[fNp2] = particle->PdgCode();
		fMCE2[fNp2] = particle->E();
	      }
	      idesitr++;
	    }
	    fNp2++;
	  }
	  
	  fN3p0 = 3* fNp0;
	  fN3p1 = 3* fNp1;
	  fN3p2 = 3* fNp2;
	  
	  fHTree->Fill();
	  itr++;
	} // loop on Hits
	//      }
      } //  loop on NTPCs
    } // loop on cryostats
	*/

