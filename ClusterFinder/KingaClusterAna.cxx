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
#include <math.h>

#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include "TDatabasePDG.h"
#include "TSystem.h"
#include "TMath.h"

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
#include "Simulation/LArVoxelCalculator.h"
#include "Simulation/LArVoxelData.h"
#include "Simulation/LArVoxelID.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"
#include "Utilities/LArProperties.h"


 
//-------------------------------------------------
cluster::KingaClusterAna::KingaClusterAna(fhicl::ParameterSet const& pset) : 
  fKingaModuleLabel         (pset.get< std::string >("KingaModuleLabel")        ),
  //fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel")        ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")     ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")     ),
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

Energy_in_Sphere=tfs->make<TH1F>("Energy_in_Sphere","Energy contained within a sphere around a vertex", 1000,0 ,1000);

M_Delta_plus_plus=tfs->make<TH1F>("M_Delta_plus_plus","Inv Mass of delta++", 1000,0 ,10);
M_Delta_plus_plus2=tfs->make<TH1F>("M_Delta_plus_plus2","Inv Mass of delta++", 1000,0 ,10);

M_Delta_plus_plus_Mother=tfs->make<TH1F>("M_Delta_plus_plus_Mother","Inv Mass of delta++ by using Mother()", 1000,0 ,5);

Ind_eng_rectangle=tfs->make<TH1F>("Ind_eng_rectangle","Energy contained within a rectangle", 2000,0 ,80000);
Coll_eng_rectangle=tfs->make<TH1F>("Coll_eng_rectangle","Energy contained within a rectangle", 2000,0 ,80000);

Ind_eng_rectangle2=tfs->make<TH1F>("Ind_eng_rectangle2","Energy contained within a rectangle", 2000,0 ,80000);
Coll_eng_rectangle2=tfs->make<TH1F>("Coll_eng_rectangle2","Energy contained within a rectangle", 2000,0 ,80000);

Ind_eng_rectangle3=tfs->make<TH1F>("Ind_eng_rectangle3","Energy contained within a rectangle", 2000,0 ,80000);
Coll_eng_rectangle3=tfs->make<TH1F>("Coll_eng_rectangle3","Energy contained within a rectangle", 2000,0 ,80000);


Number_protons=tfs->make<TH1F>("Number_protons","Number of protons with StatusCode=1 for each Event", 15,0 ,15);

Vertex_x=tfs->make<TH1F>("Vertex_x","vertex X coordinate value", 60,0 ,60);
Vertex_y=tfs->make<TH1F>("Vertex_y","vertex Y coordinate value", 60,-30 ,30);
Vertex_z=tfs->make<TH1F>("Vertex_z","vertex Z coordinate value", 100,0 ,100);

  
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
  
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);
  
  art::Handle< std::vector<sim::LArVoxelData> > vxlistHandle;
  evt.getByLabel(fLArG4ModuleLabel,vxlistHandle);
  
  sim::ParticleList _particleList = sim::SimListUtils::GetParticleList(evt, fLArG4ModuleLabel);
  
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitsModuleLabel,hitListHandle);
  
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterFinderModuleLabel,clusterListHandle);
 
  
  
  for ( sim::ParticleList::const_iterator i = _particleList.begin();
	i != _particleList.end(); ++i )
    {
       const sim::Particle* particle = (*i).second;
      int pdgcode=particle->PdgCode();
        std::cout<<" QQQQ.... pdg= "<<pdgcode<<std::endl;
        if(pdgcode==2224){std::cout<<"Victory"<<std::endl;}
      //int trackID = (*i).first;
      //std::cout<<"trackID= "<<trackID<<std::endl;
     // mc_trackids.push_back(trackID);      
    }
  
   //................................................................
double vertex [3] = { 0, 0, 0 };
int have_p=0;
int have_pion=0;
double E_p, E_pion;
int event_has_pi0=0;
int event_has_pi_plus=0;
double INV_MASS=0,INV_MASS_Mother=0;
double MC_Total_Eng=0;

   art::PtrVector<simb::MCTruth> mclist;
   for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
      mclist.push_back(mctparticle);
    } 
    
    
    
    
    
  //..........................................................
   
    //First determine what kind of event we are dealing with, and select only the ones you want:
 //..........................................................
 int no_protons=0;
 int delta_index=-999;
 
 for( unsigned int i = 0; i < mclist.size(); ++i ){
    art::Ptr<simb::MCTruth> mc(mclist[i]);
     for(int j = 0; j < mc->NParticles(); ++j){
    simb::MCParticle part(mc->GetParticle(j));
    
    if(part.PdgCode()==2212 && part.StatusCode()==1){
    //std::cout<<"and its mother = "<<part.Mother()<<std::endl;
    no_protons++;
    }
    
     if(part.PdgCode()==111 && part.StatusCode()==1){
    event_has_pi0=1;
    }
    
    if(part.PdgCode()==211 && part.StatusCode()==1){
    event_has_pi_plus=1;
    }
    
    //now let's check for delta++ and get its daughters (should be a proton and pi+)
    if(part.PdgCode()==2224){
    //pick it's trackID, remember that here it was set to be negative and they start from -1 instead of 0 so we need to account for this in order to agree with Daughter() from genie :
    
    delta_index=abs(part.TrackId())-1;
    std::cout<<"delta_index= "<<delta_index<<std::endl;
    std::cout<<"Have a delta++ and it has this many daughters: "<<part.NumberDaughters()<<std::endl;
     for(int k=0; k<part.NumberDaughters(); k++){
    std::cout<<"No of Daughters= "<<part.NumberDaughters()<<std::endl;
//      
//      
     }//for each daughter
    
    
    
    
    }
    
    
    
    }
 }
 
 Number_protons->Fill(no_protons);
 
  //..........................................................

 double mu_theta_true, mu_phi_true;
 double pion_theta_true, pion_phi_true;
 
 std::cout<<"mclist.size() = "<<mclist.size()<<std::endl;
 
 for( unsigned int i = 0; i < mclist.size(); ++i ){
    art::Ptr<simb::MCTruth> mc(mclist[i]);
    
    simb::MCParticle neut(mc->GetParticle(i));
    
    vertex[0] =neut.Vx();
    vertex[1] =neut.Vy();
    vertex[2] =neut.Vz();
    std::cout<<"neut.Vx()= "<<neut.Vx()<<" ,y= "<<neut.Vy()<<" ,z= "<<neut.Vz()<<std::endl;
    
    Vertex_x->Fill(neut.Vx());
    Vertex_y->Fill(neut.Vy());
    Vertex_z->Fill(neut.Vz());
   
    
  
      for(int j = 0; j < mc->NParticles(); ++j){
    simb::MCParticle part(mc->GetParticle(j));
 //std::cout<<"Process="<<part.Process()<<std::endl;
    std::cout<<"pdg= "<<part.PdgCode()<<" ,Process="<<part.Process()<<" StatusCode= "<<part.StatusCode()<<" mass= "<<part.Mass()<<" p= "<<part.P()<<" E= "<<part.E()<<" trackID= "<<part.TrackId()<<" ND= "<<part.NumberDaughters()<<" Mother= "<<part.Mother()<<std::endl;
    //std::cout<<"P()= "<<part.P()<<" sqrt(pow(part.Px(),2)+pow(part.Py(),2.)+pow(part.Pz(),2.))= "<<sqrt(pow(part.Px(),2.)+pow(part.Py(),2.)+pow(part.Pz(),2.))<<std::endl;
    
    //std::cout<<"Mass()= "<<part.Mass()<<" sqrt(pow(part.E(),2.)-pow(part.P(),2.))= "<<sqrt(pow(part.E(),2.)-pow(part.P(),2.))<<std::endl;
    MC_Total_Eng+=part.E();
   
//now we know which proton and pion come from delta++ so let's pick them to do invariant mass:
  if(part.Mother()==delta_index){
std::cout<<"FOUND DECAY PRODUCT OF DELTA++, PDG= "<<part.PdgCode()<<" mass= "<<part.Mass()<<" sqrt(pow(part.E(),2.)-pow(part.P(),2.))= "<<sqrt(pow(part.E(),2.)-pow(part.P(),2.))<<std::endl;
   INV_MASS_Mother+=sqrt(pow(part.E(),2.)-pow(part.P(),2.));




   }
    
   
    
    if( event_has_pi_plus==1){
    //std::cout<<"THIS EVENT HAS PI+ YAY!!!!!!!!!!!"<<std::endl;
    if(part.PdgCode()!=13 && part.PdgCode()!=-13 && part.PdgCode()!=2112 && part.PdgCode()<10000 && part.StatusCode()==1){
    std::cout<<"including pdg= "<<part.PdgCode()<<std::endl;
INV_MASS+=part.Mass();

}


    if((part.PdgCode()==2212) && (part.StatusCode()==1)){
    std::cout<<"we have a proton with E="<<part.E()<<std::endl;
    have_p++;
    E_p=part.E();
    }
    
    if((part.PdgCode()==13 || part.PdgCode()==-13) && (part.StatusCode()==1)){
   
    std::cout<<"we have a mu with E="<<part.E()<<std::endl;
    
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
    
    
    //pi- captures so only do pi+
    //.........................................................
     if((part.PdgCode()==211 ) && (part.StatusCode()==1)){
    have_pion++;
    E_pion=part.E();
    std::cout<<"we have a pion with E= "<<part.E()<<std::endl;
    
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
    //Lets draw a rectangle around a vertex, run dbscan and then loop through all the hits. If the fall within the rectangle then add the energy from them. This might be another way to distinguish RES from COH events:
    
    
    
    
    
    
    
    
    
    
    
    //................................................................

    }//event_has_pi0==0
    
 }
 
 }
 
 if( event_has_pi_plus==1){
 std::cout<<"INVARIANT MASS (BY USING STATUS CODES)= "<<INV_MASS<<" GeV"<<std::endl;
 M_Delta_plus_plus2->Fill(INV_MASS);
 
 if(INV_MASS_Mother!=0) {M_Delta_plus_plus_Mother->Fill(INV_MASS_Mother);}
 
 }
 
 
 //if we have one proton and one pion we can reconstruct the resonance for delta++
 // if(have_p+have_pion==2 && event_has_pi0==0){
//  std::cout<<"E_Delta_plus_plus= "<<E_p+E_pion<<std::endl;
//  E_Delta_plus_plus->Fill(E_p+E_pion);
//  
//  }
 
 
 
 
 //................................................................
 // Let's draw a sphere around a vertex and find out how much energy is contained in it. Should be a big difference for COH vs RES interactions.
 
 
  
  std::vector<int>::iterator it,it2;
  std::vector<int> trackIDs;
  double Inv_Mass=0;
  INV_MASS_Mother=0;
  double Eng_211=0, Eng_211_2=0;
  double Eng_13=0;
  double Tot_Energy_Event=0;
  double Eng_else=0;
  
  
  
  if( event_has_pi_plus==1){
  
  
  double Energy=0;
  double fActivityRadius=5;
  
  
  //std::cout<<"vxlistHandle->size()= "<<vxlistHandle->size()<<std::endl;
  for(unsigned int i = 0; i < vxlistHandle->size(); ++i){
    // Get the reference to the LArVoxelID in the LArVoxelList.
    art::Ptr<sim::LArVoxelData> voxel(vxlistHandle, i);
      
    int numberParticles = voxel->NumberParticles();
	      
      //std::cout<<"numberParticles "<<numberParticles<<std::endl;
      
      
      
       
    for ( int i = 0; i != numberParticles; ++i )
      {
      

      int trackID=voxel->TrackID(i);
   const sim::Particle* particle = _particleList.at(trackID);
   
  
   //if trackID does not yet exist in trackIDs vector add it
  // it=std::find(trackIDs.begin(),trackIDs.end(),trackID);
   //if(it==trackIDs.end()){ //if new trackID
   //trackIDs.push_back(trackID);
   //double p=sqrt(pow(particle->Px(),2.)+pow(particle->Py(),2.)+pow(particle->Pz(),2.));
   int pdg = particle->PdgCode();
   
  // std::cout<<"FROM VOXELS : pdg= "<<pdg<<" Process= "<<particle->Process()<<" trackID= "<<trackID<<" p="<<p<<" E= "<<particle->E()<<std::endl;
   
   if(pdg==2224){
   std::cout<<"delta++ from G4 $$$$$$$$$$$$$$$$$$$$$"<<std::endl;
   }
   
   if( _particleList.IsPrimary( trackID )){
  // std::cout<<"pdg= "<<pdg<<" trackID= "<<trackID<<" Primary"<<" with Energy= "<<particle->E()<<" and p= "<<particle->P()<<" P(3mom)="<<p<<" sqrt(pow(part.E(),2.)-pow(p,2.))= "<<sqrt(pow(particle->E(),2.)-pow(p,2.))<<" Mass= "<<particle->Mass()<< std::endl;
  
//   int mother = particle->Mother();
//      particle = _particleList.at( trackID );
//      int MPDG=particle->PdgCode();
    //std::cout<<"Primary,Mother is now (if -1 then no mother)= "<<mother<<std::endl; 
//      
   }
   
    if(! _particleList.IsPrimary( trackID )){
   //std::cout<<"pdg= "<<pdg<<" trackID= "<<trackID<<" NOT Primary"<<" with Energy= "<<particle->E()<<" and p= "<<particle->P()<<" P(3mom)="<<p<<" sqrt(pow(part.E(),2.)-pow(p,2.))= "<<sqrt(pow(particle->E(),2.)-pow(p,2.))<<" Mass= "<<particle->Mass()<<std::endl;
   //std::cout<<"BEFORE *********************"<<std::endl;
   while(! _particleList.IsPrimary( trackID ))
   {
     trackID = particle->Mother();
     particle = _particleList.at( trackID );
     int MPDG=particle->PdgCode();
     if(MPDG==2224){std::cout<<"ATTENTION, FOUND D++"<<std::endl;}
std::cout<<"MPDG is now = "<<MPDG<<std::endl;   
   
   }
      //std::cout<<"AFTER *********************"<<std::endl;

   
   
   
   }
   
   
   // if( _particleList.IsPrimary( trackID ) && (pdg==13 || pdg==-13)){
//  //  std::cout<<"IT'S A PRIMARY MUON!!! DON'T CONSIDER"<<std::endl;
//    continue;
//    }
   
   
   
  //  if(pdg!=13 && pdg!=-13){
//    std::cout<<"pdg= "<<pdg<<" Inv_MAss= "<<particle->Mass()<<std::endl;
//    Inv_Mass+=particle->Mass();
//    }
  // }
  
  
  while ( (! _particleList.IsPrimary( trackID )) && (((_particleList.at(particle->Mother()))->PdgCode())==pdg))
			    {
			    trackID = particle->Mother();
			    particle = _particleList.at( trackID );
			    
			    }
  //std::cout<<"The trackID is now = "<<trackID<<std::endl;
  it2=std::find(trackIDs.begin(),trackIDs.end(),trackID);
   if(it2==trackIDs.end()){ //if new trackID
  trackIDs.push_back(trackID);
//  std::cout<<"ADDING INV_MASS of pdg= "<<pdg<<" with trackID= "<<trackID<<" and Mass= "<<particle->Mass()<<std::endl;
  
  Inv_Mass+=particle->Mass();
  if(pdg==211){ Eng_211+=voxel->Energy(i)*1000.;}
  }
  //}//only new trackID
  
//      if(_particleList.IsPrimary( trackID )){

      if(pdg==13 || pdg==-13){ Eng_13+=voxel->Energy(i)*1000.;}
      if(pdg==211){ Eng_211_2+=voxel->Energy(i)*1000.;}
      if(pdg!=211 && pdg!=13){ Eng_else+=voxel->Energy(i)*1000.;}
	  Tot_Energy_Event+= (voxel->Energy(i)*1000.);
	 if(sqrt(pow(TMath::Abs(voxel->VoxelID().X()-vertex[0]),2)+pow(TMath::Abs(voxel->VoxelID().Y()-vertex[1]),2)+pow(TMath::Abs(voxel->VoxelID().Z()-vertex[2]),2))<fActivityRadius){
	 
		
		 Energy+= (voxel->Energy(i)*1000.);
		
		 }
		//}//if primary 
 
 }
 
 }
 
  std::cout<<"*** TOTAL ENERGY in the EVENT from MC "<<MC_Total_Eng<<std::endl;	

 std::cout<<"*** TOTAL ENERGY in the EVENT "<<Tot_Energy_Event<<std::endl;	
 std::cout<<"*** TOTAL 211 ENERGY= "<<Eng_211<<std::endl;
 std::cout<<"*** TOTAL 211_2 ENERGY= "<<Eng_211_2<<std::endl;
 std::cout<<"*** TOTAL 13 ENERGY= "<<Eng_13<<std::endl;
 std::cout<<"*** TOTAL else ENERGY= "<<Eng_else<<std::endl;
 std::cout<<"Total Inv_Mass= "<<Inv_Mass<<std::endl;
 M_Delta_plus_plus->Fill(Inv_Mass);

 std::cout<<"ENERGY= "<<Energy<<std::endl;
  Energy_in_Sphere->Fill(Energy);
  
  
  }// event_has_pi0==0 exclude RES CCpi0 events
 
 //................................................................

 //  art::PtrVector<recob::Hit> hits;
//   for (unsigned int ii = 0; ii <  hitListHandle->size(); ++ii)
//     {
//       art::Ptr<recob::Hit> hitHolder(hitListHandle,ii);
//       hits.push_back(hitHolder);
//     }
// 
// std::cout<<"trying to get  DBSCAN Clusters"<<std::endl;
//   art::PtrVector<recob::Cluster> clusters;
//   for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii)
//     {
//       art::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
//       clusters.push_back(clusterHolder);
//     }
// std::cout<<"GOT DBSCAN Clusters"<<std::endl;
// std::cout<<"trying to get KingaClusters"<<std::endl;
// art::PtrVector<recob::Cluster> kingaclusters;
//   for (unsigned int ii = 0; ii <  kingaListHandle->size(); ++ii)
//     {
//       art::Ptr<recob::Cluster> kingaHolder(kingaListHandle,ii);
//       kingaclusters.push_back(kingaHolder);
//     }
//     
//     std::cout<<"GOT KingaClusters"<<std::endl;
//     
//   std::cout<<"in Efficiency, kingaclusters.size()= "<<kingaclusters.size()<<std::endl;
//   std::cout<<"in Efficiency, dbscanclusters.size()= "<<clusters.size()<<std::endl;

  
  
  
  
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









// Let's get dbscan's hits and see how much energy is deposited around a rectangle centered at the vertex. Should see a difference for RES vs COH events:





art::ServiceHandle<util::LArProperties> larp;
double drifttick=(vertex[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./.198);
std::cout<<"drifttick= "<<drifttick<<std::endl;
ftime_vertex.push_back(drifttick);
art::PtrVector<recob::Cluster> clusIn;
art::PtrVector<recob::Hit> hits;
ftimetick      =  0.198; 
fdriftvelocity =  0.157; 
double ind_hits_eng_rectangle=0, ind_hits_eng_rectangle2=0, ind_hits_eng_rectangle3=0;
double coll_hits_eng_rectangle=0, coll_hits_eng_rectangle2=0,coll_hits_eng_rectangle3=0;
//double MCvertex [3];


 for(unsigned int ii = 0; ii < clusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster(clusterListHandle, ii);
      clusIn.push_back(cluster);
    }
std::cout<<"No of DBSCAN clusters= "<<clusIn.size()<<std::endl;



unsigned int p(0),w(0), channel(0), wirevertex, channel2;

std::cout<<"No of planes = "<<geom->Nplanes()<<std::endl;
  for(unsigned int plane = 0; plane < geom->Nplanes(); plane++) {



if(plane==0){
	vertex[0]=.3;//force time coordinate to be closer to induction plane 
	}
      else{
	vertex[0]=-.3;//force time coordinate to be closer to collection plane
     }
      channel2 = geom->NearestChannel(vertex);
      geom->ChannelToWire(channel2,plane,wirevertex); 
      std::cout<<"wirevertex= "<<wirevertex<<std::endl;
      fwire_vertex.push_back(wirevertex);




    for(unsigned int j=0; j<clusIn.size();++j) {
   
    hits=clusIn[j]->Hits();
    for(unsigned int i = 0; i< hits.size(); ++i){
    channel=hits[i]->Wire()->RawDigit()->Channel();
    geom->ChannelToWire(channel,p,w);

    if(p == plane){
    allhits.push_back(hits[i]);
    //std::cout<<"plane= "<<plane<<" wire= "<<w<<" time= "<<hits[i]->PeakTime()<<std::endl; 
   // std::cout<<"w= "<<w<<"  fwire_vertex["<<plane<<"] = "<<fwire_vertex[plane]<<" fabs(w*0.4-fwire_vertex[plane]*0.4)= "<<fabs(w*0.4-fwire_vertex[plane]*0.4)<<std::endl;
    
    //std::cout<<"allhits[i]->StartTime()+allhits[i]->EndTime())/2.= "<<(allhits[i]->StartTime()+allhits[i]->EndTime())/2.<<" PeakTime= "<<allhits[i]->PeakTime()<<" ftime_vertex[0]= "<<ftime_vertex[0]<<" fabs(((allhits[i]->StartTime()+allhits[i]->EndTime())/2.)-ftime_vertex[0])* ftimetick *fdriftvelocity= "<<fabs(((allhits[i]->StartTime()+allhits[i]->EndTime())/2.)-ftime_vertex[0])* ftimetick *fdriftvelocity<<std::endl;
    
    
    if(p==plane && p==0 && fabs(w*0.4-fwire_vertex[plane]*0.4)<2 && fabs(((allhits[i]->StartTime()+allhits[i]->EndTime())/2.)-ftime_vertex[0])* ftimetick *fdriftvelocity<5){
    std::cout<<"filling Ind_eng_rectangle"<<std::endl;
    ind_hits_eng_rectangle+=allhits[i]->Charge();

    }
    
     if(p==plane && p==1 && fabs(w*0.4-fwire_vertex[plane]*0.4)<2 && fabs(((allhits[i]->StartTime()+allhits[i]->EndTime())/2.)-ftime_vertex[0])* ftimetick *fdriftvelocity<5){
    std::cout<<"filling Coll_eng_rectangle"<<std::endl;
    coll_hits_eng_rectangle+=allhits[i]->Charge();
    
    }
    //................................
    
    if(p==plane && p==0 && fabs(w*0.4-fwire_vertex[plane]*0.4)<3 && fabs(((allhits[i]->StartTime()+allhits[i]->EndTime())/2.)-ftime_vertex[0])* ftimetick *fdriftvelocity<6){
    std::cout<<"filling Ind_eng_rectangle2"<<std::endl;
    ind_hits_eng_rectangle2+=allhits[i]->Charge();

    }
    
     if(p==plane && p==1 && fabs(w*0.4-fwire_vertex[plane]*0.4)<3 && fabs(((allhits[i]->StartTime()+allhits[i]->EndTime())/2.)-ftime_vertex[0])* ftimetick *fdriftvelocity<6){
    std::cout<<"filling Coll_eng_rectangle2"<<std::endl;
    coll_hits_eng_rectangle2+=allhits[i]->Charge();
    
    }
    
    
     //................................
    
    if(p==plane && p==0 && fabs(w*0.4-fwire_vertex[plane]*0.4)<6 && fabs(((allhits[i]->StartTime()+allhits[i]->EndTime())/2.)-ftime_vertex[0])* ftimetick *fdriftvelocity<7){
    std::cout<<"filling Ind_eng_rectangle2"<<std::endl;
    ind_hits_eng_rectangle3+=allhits[i]->Charge();

    }
    
     if(p==plane && p==1 && fabs(w*0.4-fwire_vertex[plane]*0.4)<6 && fabs(((allhits[i]->StartTime()+allhits[i]->EndTime())/2.)-ftime_vertex[0])* ftimetick *fdriftvelocity<7){
    std::cout<<"filling Coll_eng_rectangle2"<<std::endl;
    coll_hits_eng_rectangle3+=allhits[i]->Charge();
    
    }
    
    
    
    
    
    }//if p==plane
   
    }//loop through hits for each cluster
   // std::cout<<"hits.size()= "<<hits.size()<<std::endl;
    }//loop through each cluster
   
    std::cout<<"allhits.size()="<<allhits.size()<<std::endl;
   
   
    //Now we have hits for the plane that we are on right now, so let's do some work:
    if(allhits.size()>10){
if(plane==0){
std::cout<<"Final ind_hits_eng_rectangle= "<<ind_hits_eng_rectangle<<std::endl;
Ind_eng_rectangle->Fill(ind_hits_eng_rectangle);

std::cout<<"Final ind_hits_eng_rectangle2= "<<ind_hits_eng_rectangle2<<std::endl;
Ind_eng_rectangle2->Fill(ind_hits_eng_rectangle2);

std::cout<<"Final ind_hits_eng_rectangle3= "<<ind_hits_eng_rectangle3<<std::endl;
Ind_eng_rectangle3->Fill(ind_hits_eng_rectangle3);
}
    
if(plane==1){
std::cout<<"Final coll_hits_eng_rectangle= "<<coll_hits_eng_rectangle<<std::endl;
 Coll_eng_rectangle->Fill(coll_hits_eng_rectangle);
 
 std::cout<<"Final coll_hits_eng_rectangle2= "<<coll_hits_eng_rectangle2<<std::endl;
 Coll_eng_rectangle2->Fill(coll_hits_eng_rectangle2);
 
 std::cout<<"Final coll_hits_eng_rectangle3= "<<coll_hits_eng_rectangle3<<std::endl;
 Coll_eng_rectangle3->Fill(coll_hits_eng_rectangle3);

}
}//if allhits.size()>10
allhits.clear();

}//for each plane


std::cout<<"ORIG,Final ind_hits_eng_rectangle= "<<ind_hits_eng_rectangle<<std::endl;
std::cout<<"ORIG,Final coll_hits_eng_rectangle= "<<coll_hits_eng_rectangle<<std::endl;
   
    ind_hits_eng_rectangle=0;
    coll_hits_eng_rectangle=0;
    ind_hits_eng_rectangle2=0;
    coll_hits_eng_rectangle2=0;
    ind_hits_eng_rectangle3=0;
    coll_hits_eng_rectangle3=0;
    allhits.clear();
    ftime_vertex.clear();
    fwire_vertex.clear();












 
}

