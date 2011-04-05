////////////////////////////////////////////////////////////////////////
//
// DBSCAN analyzer
//
// \author kinga.partyka@yale.edu
// echurch@fnal.gov
// 
////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional>

#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TDatabasePDG.h"
#include "TSystem.h"

#include "art/Framework/Core/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Core/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 


#include "ClusterFinder/KingaClusterAna.h"
#include "Geometry/geo.h"
#include "SimulationBase/simbase.h"
#include "Simulation/sim.h"
#include "Simulation/SimListUtils.h"
#include "RecoBase/recobase.h"
#include "RawData/RawDigit.h"


 
//-------------------------------------------------
cluster::KingaClusterAna::KingaClusterAna(fhicl::ParameterSet const& pset) : 
  fKingaModuleLabel         (pset.get< std::string >("KingaModuleLabel")        ),
  //fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel")        ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")     ),
  //fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")         ),
  fClusterFinderModuleLabel (pset.get< std::string >("ClusterFinderModuleLabel"))
{



}

//-------------------------------------------------
cluster::KingaClusterAna::~KingaClusterAna()
{

}

void cluster::KingaClusterAna::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;
Mu_theta=tfs->make<TH1F>("Mu_theta","Muon theta angle", 360,0 ,360);
Mu_phi=tfs->make<TH1F>("Mu_phi","Muon phi angle", 360,0 ,360);
Mu_phi_oneside=tfs->make<TH1F>("Mu_phi_oneside","Muon phi angle", 360,0 ,360);

pion_theta=tfs->make<TH1F>("pion_theta","Pion theta angle", 360,0 ,360);
pion_phi=tfs->make<TH1F>("pion_phi","Pion phi angle", 360,0 ,360);
pion_phi_oneside=tfs->make<TH1F>("pion_phi_oneside","Pion phi angle", 360,0 ,360);
  
}

void cluster::KingaClusterAna::analyze(const art::Event& evt)
{
  std::cout<<"Hello, You are in KingaClusterAna::analyze"<<std::endl;
  std::cout << "run    : " << evt.run() << std::endl;
  //std::cout << "subrun : " << evt.subRun() << std::endl; // Doesn't compile w. or w.o. id().
  std::cout << "event  : " << evt.id().event() << std::endl;
  //----------------------------------------------------------------

  /* This is basically a module for studying MC efficiency/purity. Kick out now if not MC. EC, 8-Oct-2010 */
  if (evt.isRealData()) 
    {
      std::cout<<"**** KingaClusterAna: Bailing. Don't call this module if you're not MC. "<<std::endl;
      return;
    }

  art::Handle< std::vector<recob::Cluster>  > kingaListHandle;
  evt.getByLabel(fKingaModuleLabel,kingaListHandle);
  
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitsModuleLabel,hitListHandle);
  
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterFinderModuleLabel,clusterListHandle);
 
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);
  
  // sim::ParticleList _particleList = sim::SimListUtils::GetParticleList(evt, fLArG4ModuleLabel);
  
   //................................................................

   art::PtrVector<simb::MCTruth> mclist;
   for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
      mclist.push_back(mctparticle);
    } 
 
 double mu_theta_true, mu_phi_true;
 double pion_theta_true, pion_phi_true;
 
 for( unsigned int i = 0; i < mclist.size(); ++i ){
    art::Ptr<simb::MCTruth> mc(mclist[i]);
    
   // std::cout<<"Number of Primaries= "<<mclist->NumberOfPrimaries()<<std::endl;
    
    for(int j = 0; j < mc->NParticles(); ++j){
    simb::MCParticle part(mc->GetParticle(j));
 //std::cout<<"Process="<<part.Process()<<std::endl;
    std::cout<<"pdg= "<<part.PdgCode()<<" ,Process="<<part.Process()<<" StatusCode= "<<part.StatusCode()<<std::endl;
    
    if((part.PdgCode()==13 || part.PdgCode()==-13) && (part.StatusCode()==1)){
    
    std::cout<<"we have a mu"<<std::endl;
    
    mu_theta_true=(TMath::ACos(part.Pz()/sqrt(pow(part.Px(),2)+pow(part.Py(),2)+pow(part.Pz(),2))))*(180/TMath::Pi());
    std::cout<<"mu_theta_true= "<<mu_theta_true<<std::endl;
    
    mu_phi_true=(TMath::Pi()+TMath::ATan2(-part.Py(),-part.Px()))*(180/TMath::Pi());
 std::cout<<"mu_phi_true= "<<mu_phi_true<<std::endl;
 
 Mu_theta->Fill(mu_theta_true);
 Mu_phi->Fill(mu_phi_true);
 if(mu_phi_true<=180){
 Mu_phi_oneside->Fill(mu_phi_true);}
 if(mu_phi_true>180){
 Mu_phi_oneside->Fill(360-mu_phi_true);}
 
 
    }
    
    
    
    //.........................................................
     if((part.PdgCode()==211 || part.PdgCode()==-211 || part.PdgCode()==111) && (part.StatusCode()==1)){
    
    std::cout<<"we have a pion"<<std::endl;
    
    pion_theta_true=(TMath::ACos(part.Pz()/sqrt(pow(part.Px(),2)+pow(part.Py(),2)+pow(part.Pz(),2))))*(180/TMath::Pi());
    std::cout<<"pion_theta_true= "<<pion_theta_true<<std::endl;
    
    pion_phi_true=(TMath::Pi()+TMath::ATan2(-part.Py(),-part.Px()))*(180/TMath::Pi());
  std::cout<<"pion_phi_true= "<<pion_phi_true<<std::endl;
 
  pion_theta->Fill(pion_theta_true);
  pion_phi->Fill(pion_phi_true);
  if(pion_phi_true<=180){
  pion_phi_oneside->Fill(pion_phi_true);}
  if(pion_phi_true>180){
  pion_phi_oneside->Fill(360-pion_phi_true);}
 
 
    }
    //................................................................
    
    
    
 }
 }
 
 
 
 
 
 
 
 
 //................................................................

  art::PtrVector<recob::Hit> hits;
  for (unsigned int ii = 0; ii <  hitListHandle->size(); ++ii)
    {
      art::Ptr<recob::Hit> hitHolder(hitListHandle,ii);
      hits.push_back(hitHolder);
    }

std::cout<<"trying to get  DBSCAN Clusters"<<std::endl;
  art::PtrVector<recob::Cluster> clusters;
  for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
      clusters.push_back(clusterHolder);
    }
std::cout<<"GOT DBSCAN Clusters"<<std::endl;
std::cout<<"trying to get KingaClusters"<<std::endl;
art::PtrVector<recob::Cluster> kingaclusters;
  for (unsigned int ii = 0; ii <  kingaListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> kingaHolder(kingaListHandle,ii);
      kingaclusters.push_back(kingaHolder);
    }
    
    std::cout<<"GOT KingaClusters"<<std::endl;
    
  std::cout<<"in Efficiency, kingaclusters.size()= "<<kingaclusters.size()<<std::endl;
  std::cout<<"in Efficiency, dbscanclusters.size()= "<<clusters.size()<<std::endl;

  
  
  
  
  art::ServiceHandle<geo::Geometry> geom;  
 
 // if(kingaclusters.size()!=0 && hits.size()!=0)// {
//     for(unsigned int plane=0;plane<geom->Nplanes();++plane){
//       geo::View_t view = geom->Plane(plane).View();
//          
//       for(unsigned int j=0; j<clusters.size();++j) 
// 	
// 	{
// 	 
// 	  if( kingaclusters[j]->View() == view){
// 	    art::PtrVector<recob::Hit> _hits; 
// 	    
// 	   
// 	    _hits=kingaclusters[j]->Hits();
// 	    
// 	    if(_hits.size()!=0){ //need this b/c of plane
// 	      
// 	      for(unsigned int i = 0; i < _hits.size(); ++i) {
// 		
// 		
// 		
// 		
// 		double XTime=_hits[i]->PeakTime();
// 		
// 		
// 		unsigned int channel = _hits[i]->Wire()->RawDigit()->Channel();
// 		
// 		
// 		
// 		
// 	   
// 
// 
// 		
// 	   
// 	      }//for hits
// 	
// 	    
// 	  
// 	     
// 	  
// 	  
// 	   
// 	      
// 	  
// 	  
// 	  
// 	    }//non-zero hits
// 	  }//end if cluster is in correct view
// 	 
// 	
// 	}//for each cluster
//       
//      
//      
//       
//      
//  
//  
//     
// 
//     }//for each plane
//   
//       
//   
//   }

 
}

