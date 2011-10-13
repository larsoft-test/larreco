////////////////////////////////////////////////////////////////////////
//
// FinalStateParticleFilter class:
// Algoritm to produce a filtered event file having
// events with user-defined final state particles 
//
// saima@ksu.edu
//
////////////////////////////////////////////////////////////////////////

//Framework Includes
#include "art/Framework/Core/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Core/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 


//Larsoft Includes
#include "FinalStateParticleFilter.h"
#include "TDatabasePDG.h"
#include "SimulationBase/MCNeutrino.h"
#include "SimulationBase/simbase.h"

namespace filt{

  //-------------------------------------------------
  FinalStateParticleFilter::FinalStateParticleFilter(fhicl::ParameterSet const & pset)  
  {   
    this->reconfigure(pset);
  }

  //-------------------------------------------------
  FinalStateParticleFilter::~FinalStateParticleFilter()
  {
  }
  
  //-------------------------------------------------
  void FinalStateParticleFilter::reconfigure(fhicl::ParameterSet const& p)
  {
    fGenieModuleLabel = p.get< std::string      >("GenieModuleLabel");
    fPDG              = p.get< std::vector<int> >("PDG");
  } 

  //-------------------------------------------------
  void FinalStateParticleFilter::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    fPDGEvents = tfs->make<TH1D>("fPDGEvents", "Number of Selected of Events", 3, 0, 3); //counts the #events that go to the filtered file
    fTotalEvents = tfs->make<TH1D>("fTotalEvents", "Total Events", 3, 0, 3); //counts the initial #events in the unfiltered root input file
  }

  //-------------------------------------------------
  bool FinalStateParticleFilter::filter(art::Event &evt)
  { 
    
    //const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
    
    art::Handle< std::vector<simb::MCTruth> > mclist;
    evt.getByLabel(fGenieModuleLabel,mclist);
    art::Ptr<simb::MCTruth> mc(mclist,0);
    
    fTotalEvents->Fill(1);

    std::vector<int> finalstateparticles;

    //get a vector of final state particles   
    for(int i = 0; i < mc->NParticles(); ++i){
      simb::MCParticle part(mc->GetParticle(i));
      if(part.StatusCode()== 1)
	finalstateparticles.push_back(part.PdgCode());
    }

    return isSubset(fPDG, finalstateparticles); 

  } // bool  
} // namespace
    
//------------------------------------------------   

bool filt::FinalStateParticleFilter::isSubset(std::vector<int>& a, std::vector<int>& b)
{
  for (std::vector<int>::iterator i = a.begin(); i != a.end(); i++){
    bool found = false;
    for (std::vector<int>::iterator j = b.begin(); j != b.end(); j++){
      if (*i == *j){
	found = true;
	break;
      }
    }
      
    if (!found){
      return false;
    }
  }
  std::cout << "this is a selected event " << std::endl;
  return true;
}

//--------------------------------------------------
