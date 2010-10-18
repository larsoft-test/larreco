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
#include "Geometry/geo.h"
#include "DBSCANfinderAna.h"

#include "TDatabasePDG.h"
#include "TSystem.h"


#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/ParameterSet/interface/ParameterSet.h" 
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h" 
#include "DataFormats/Common/interface/Handle.h" 
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ptr.h" 
#include "DataFormats/Common/interface/PtrVector.h" 
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "FWCore/Services/interface/TFileService.h" 
#include "FWCore/Framework/interface/TFileDirectory.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h" 

#include "Geometry/Geometry.h"
#include "SimulationBase/MCFlux.h"
#include "SimulationBase/MCNeutrino.h"
#include "SimulationBase/MCTruth.h"
#include "RawData/RawDigit.h" 
#include "Simulation/Particle.h"
#include "Simulation/ParticleList.h"
#include "Simulation/LArVoxelList.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include <Simulation/LArVoxelID.h>
#include <Simulation/Electrons.h>
#include <Simulation/SimDigit.h>
#include "Geometry/WireGeo.h"
#include <SimulationBase/MCParticle.h>


//-------------------------------------------------
cluster::DBSCANfinderAna::DBSCANfinderAna(edm::ParameterSet const& pset) : 
  fDigitModuleLabel    (pset.getParameter< std::string >("DigitModuleLabel")),
  fHitModuleLabel    (pset.getParameter< std::string >("HitModuleLabel")),
  fLArG4ModuleLabel    (pset.getParameter< std::string >("LArG4ModuleLabel")),
  fDetSimModuleLabel    (pset.getParameter< std::string >("DetSimModuleLabel")),
  fGenieGenModuleLabel    (pset.getParameter< std::string >("GenieGenModuleLabel")),
  fClusterFinderModuleLabel    (pset.getParameter< std::string >("ClusterFinderModuleLabel"))
{



}

//-------------------------------------------------
cluster::DBSCANfinderAna::~DBSCANfinderAna()
{

}

void cluster::DBSCANfinderAna::beginJob(edm::EventSetup const&)
{
     // get access to the TFile service
     edm::Service<edm::TFileService> tfs;

  fNoParticles_pdg_per_event = tfs->make<TH1F>("fNoParticles_pdg_per_event","Average # of Particles per cluster for each event", 500,0 ,5);
    fNoParticles_pdg=tfs->make<TH1F>("fNoParticles_pdg","Number of Particles in a Cluster for each cluster", 500,0 ,5);
    fNoParticles_trackid=tfs->make<TH1F>("fNoParticles_trackid","Number of different TrackIDs in a Cluster", 300,0 ,30);

fNoParticles_trackid_mother=tfs->make<TH1F>("fNoParticles_trackid_mother","Number of different TrackIDs in a Cluster(using mother)for each cluster", 300,0 ,30);

fNoParticles_trackid_per_event=tfs->make<TH1F>("fNoParticles_trackid_per_event","Avg Number of different TrackIDs per Cluster per event", 300,0 ,30);
fCl_for_Muon=tfs->make<TH1F>("fCl_for_Muon","Number of Clusters for Muon per plane (pdg)", 1500,0 ,15);
//  fCl_for_Electron=tfs->make<TH1F>("fCl_for_Electron","Number of Clusters for Electron  (pdg)", 1500,0 ,15);
//  fCl_for_Positron=tfs->make<TH1F>("fCl_for_Positron","Number of Clusters for Positron", 1500,0 ,15);
//  fCl_for_Pion_111=tfs->make<TH1F>("fCl_for_Pion_111","Number of Clusters for Pion (111)", 1500,0 ,15);
//  fCl_for_Pion_211=tfs->make<TH1F>("fCl_for_Pion_211","Number of Clusters for Pion (211)", 1500,0 ,15);
//  fCl_for_Pion_m211=tfs->make<TH1F>("fCl_for_Pion_m211","Number of Clusters for Pion (-211)", 1500,0 ,15);
// fCl_for_Proton=tfs->make<TH1F>("fCl_for_Proton","Number of Clusters for Proton", 1500,0 ,15);

fNoClustersInEvent=tfs->make<TH1F>("fNoClustersInEvent","Number of Clusters in an Event", 5000,0 ,50);

 fPercentNoise=tfs->make<TH1F>("fPercentNoise","% of hits that were marked as Noise by DBSCAN",2500,0 ,25);

fno_of_clusters_per_track=tfs->make<TH1F>("fno_of_clusters_per_track","Number of Clusters per TrackID per plane", 1500,0 ,15);

fPercent_lost_muon_hits=tfs->make<TH1F>("fPercent_lost_muon_hits","Number of muon hits excluded by dbscan in % (per Event)", 10000,0 ,100);
fPercent_lost_electron_hits=tfs->make<TH1F>("fPercent_lost_electron_hits","Number of electron hits excluded by dbscan in % (per Event)", 10000,0 ,100);
fPercent_lost_positron_hits=tfs->make<TH1F>("fPercent_lost_positron_hits","Number of positron hits excluded by dbscan in % (per Event)", 10000,0 ,100);
fPercent_lost_111_hits=tfs->make<TH1F>("fPercent_lost_111_hits","Number of pion(111) hits excluded by dbscan in % (per Event)", 10000,0 ,100);
fPercent_lost_211_hits=tfs->make<TH1F>("fPercent_lost_211_hits","Number of pion(211) hits excluded by dbscan in % (per Event)", 10000,0 ,100);
fPercent_lost_m211_hits=tfs->make<TH1F>("fPercent_lost_m211_hits","Number of pion(-211) hits excluded by dbscan in % (per Event)", 10000,0 ,100);
fPercent_lost_2212_hits=tfs->make<TH1F>("fPercent_lost_2212_hits","Number of proton hits excluded by dbscan in % (per Event)", 10000,0 ,100);
fPercent_lost_2112_hits=tfs->make<TH1F>("fPercent_lost_2112_hits","Number of neutron hits excluded by dbscan in % (per Event)", 10000,0 ,100);

fPercent_lost_muon_energy=tfs->make<TH1F>("fPercent_lost_muon_energy"," muon energy excluded by dbscan in % (per Event)", 10000,0 ,100);
fPercent_lost_electron_energy=tfs->make<TH1F>("fPercent_lost_electron_energy","electron energy excluded by dbscan in % (per Event)", 10000,0 ,100);
fPercent_lost_positron_energy=tfs->make<TH1F>("fPercent_lost_positron_energy"," positron energy excluded by dbscan in % (per Event)", 10000,0 ,100);
fPercent_lost_111_energy=tfs->make<TH1F>("fPercent_lost_111_energy","pion(111) energy excluded by dbscan in % (per Event)", 10000,0 ,100);
fPercent_lost_211_energy=tfs->make<TH1F>("fPercent_lost_211_energy","pion(211) energy excluded by dbscan in % (per Event)", 10000,0 ,100);
fPercent_lost_m211_energy=tfs->make<TH1F>("fPercent_lost_m211_energy"," pion(-211) energy excluded by dbscan in % (per Event)", 10000,0 ,100);
fPercent_lost_2212_energy=tfs->make<TH1F>("fPercent_lost_2212_energy","proton energy excluded by dbscan in % (per Event)", 10000,0 ,100);
fPercent_lost_2112_energy=tfs->make<TH1F>("fPercent_lost_2112_energy","neutron energy excluded by dbscan in % (per Event)", 10000,0 ,100);

fEnergy=tfs->make<TH1F>("fEnergy","energy for each voxel", 100000,0 ,0.0005);

 fbrian_in = tfs->make<TH2F>("fbrian_in", ";# Electrons deposited; # Electrons detected by hitfinder", 1000,     0, 10000000, 1000, 0, 10000000);
fbrian_coll = tfs->make<TH2F>("fbrian_coll", ";# Electrons deposited; # Electrons detected by hitfinder", 1000,     0, 10000000, 1000, 0, 10000000);

}

void cluster::DBSCANfinderAna::analyze(const edm::Event& evt,  edm::EventSetup const&)
{
  
  std::cout << "run    : " << evt.run() << std::endl;
  //std::cout << "subrun : " << evt.subRun() << std::endl; // Doesn't compile w. or w.o. id().
  std::cout << "event  : " << evt.id().event() << std::endl;
  //----------------------------------------------------------------

  /* This is basically a module for studying MC efficiency/purity. Kick out now if not MC. EC, 8-Oct-2010 */
  if (evt.isRealData()) 
    {
      std::cout<<"**** DBSCANfinderAna: Bailing. Don't call this module if you're not MC. "<<std::endl;
      exit (1);
    }
  //  edm::Handle< edm::View <std::vector<raw::RawDigit> > > rdListHandle;
  edm::Handle< edm::View <std::vector<sim::SimDigit> > > rdListHandle;
  evt.getByLabel(fDigitModuleLabel,rdListHandle);
  edm::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitModuleLabel,hitListHandle);
  edm::Handle< std::vector<sim::ParticleList> > partListHandle;
  evt.getByLabel(fLArG4ModuleLabel,partListHandle);
  edm::Handle< std::vector<sim::LArVoxelList> > voxelListHandle;
  evt.getByLabel(fLArG4ModuleLabel,voxelListHandle);
  edm::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);
  edm::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterFinderModuleLabel,clusterListHandle);
  edm::Handle< std::vector<recob::Wire> > wireListHandle;
  evt.getByLabel(fDetSimModuleLabel,wireListHandle);



  //  std::cout<<"****simdigit.size()= "<<simdigit.size()<<std::endl;
  //   for(int i=0; i<simdigit.size();i++){
  //     std::cout<<"***numofele: "<<simdigit[i]->NumberOfElectrons()<<"  ";}
  
  // std::cout<<"simdigit->getelecetrons: "<<*simdigit[0]->GetElectrons(0)<<"  ";
  
  //----------------------------------------------------------------

  //----------------------------------------------------------------   

  /*
  edm::PtrVector<raw::RawDigit> rawdigits;
  for (unsigned int ii = 0; ii <  rdListHandle->size(); ++ii)
    {
      edm::Ptr<raw::RawDigit> rawdigit(rdListHandle,ii);
      rawdigits.push_back(rawdigit);
    }
  */

  edm::Ptr<sim::ParticleList> check_particleList;
  edm::Ptr<sim::ParticleList> particleList;
  edm::PtrVector<sim::ParticleList> _particleList;
  for (unsigned int ii = 0; ii <  partListHandle->size(); ++ii)
    {
      edm::Ptr<sim::ParticleList> partlist(partListHandle,ii);
      _particleList.push_back(partlist);
    }
  
  std::vector<int> mc_trackids;
  
  check_particleList = _particleList[0];
  //std::cout<<"checking trackID: ";
  for ( sim::ParticleList::const_iterator i = check_particleList->begin();
	i != check_particleList->end(); ++i )
    {
      // const sim::Particle* particle = (*i).second;
      //int pdgcode=particle->PdgCode();
      //  std::cout<<"pdg= "<<pdgcode<<" ";
      int trackID = (*i).first;
      //std::cout<<"trackID= "<<trackID<<std::endl;
      mc_trackids.push_back(trackID);      
    }
  
  // std::cout<<std::endl;
   // std::cout<<" Track ID list from MC: "<<std::endl;
 //  for(unsigned int i=0;i<mc_trackids.size();++i){
  
//      std::cout<<mc_trackids[i]<<" ";
    
//   }
  // std::cout<<"I have in total "<<mc_trackids.size()<<" different tracks"<<std::endl;
  
  //----------------------------------------------------------------
  
  edm::PtrVector<sim::LArVoxelList> larVoxelList;
  for (unsigned int ii = 0; ii <  voxelListHandle->size(); ++ii)
    {
      edm::Ptr<sim::LArVoxelList> voxellist(voxelListHandle,ii);
      larVoxelList.push_back(voxellist);
    }
  // There's probably only one LArVoxelList per event, but ART-nee'-FMWK... blah blah blah
  edm::Ptr<sim::LArVoxelList> voxelList(voxelListHandle,larVoxelList.size()-1);
 
  
  //---------------------------------------------------------------- 
  edm::PtrVector<simb::MCTruth> mclist;
  for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      edm::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
      mclist.push_back(mctparticle);
    }

  edm::PtrVector<recob::Hit> hits;
  for (unsigned int ii = 0; ii <  hitListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Hit> hitHolder(hitListHandle,ii);
      hits.push_back(hitHolder);
    }


  edm::PtrVector<recob::Cluster> clusters;
  for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
      clusters.push_back(clusterHolder);
    }

  std::cout<<"in Efficiency, clusters.size()= "<<clusters.size()<<std::endl;
  //---------------------------------------------------------------
  edm::PtrVector<recob::Wire> wirelist;
  for (unsigned int ii = 0; ii <  wireListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Wire> wireHolder(wireListHandle,ii);
      wirelist.push_back(wireHolder);
    }
  
  
  //...........................................................................
  // How many different particles do we have in a cluster???
  //  -- I will answer that by 3 methods
  //  1) count different TrackIDs in each cluster
  //  2)  count different TrackIDs in each cluster with the use of "Mother" to get rid of the ones that just randomly got their TrackIDs changed by Geant4 
  //  3) count different pdg codes in each cluster <--probably not very good b/c if you have 2 "different" electrons in a cluster it's going to count them as one.
  //..........................................................................
  // How many clusters it takes to contain a particle???
  // take each TrackID and count in how many clusters it appears
  //.........................................................................
  
  
  double no_of_particles_in_cluster=0;
  double sum_vec_trackid=0;
  double no_of_clusters=0;
  double total_no_hits_in_clusters=0;
  //unsigned int plane=0;
  //  edm::Ptr<sim::SimDigit > _rawdigit;
  //edm::Ptr<sim::SimDigit > _rawdigit2;
  const sim::Electrons* _electrons=0;
  const sim::Electrons* electrons=0;
  std::vector<int> vec_pdg;
  std::vector<int> vec_trackid,vec_trackid_mother, vec_trackid_mother_en;
  std::vector<int> all_trackids;
  std::vector<int> ids;
  std::vector<int>::iterator it,it2,it3,it4,it5,it6,it7,it8;
  int no_cl_for_muon=0;
  int no_cl_for_electron=0;
  int no_cl_for_positron=0;
  int no_cl_for_pion_111=0;
  int no_cl_for_pion_211=0;
  int no_cl_for_pion_m211=0;
  int no_cl_for_proton=0;
  double noCluster=0;
  //  int muon=0,electron=0,positron=0,pion=0;
  double _hit_13=0,_hit_11=0,_hit_m_11=0,_hit_111=0,_hit_22=0,_hit_211=0,_hit_m211=0,_hit_2212=0,_hit_2112=0;
  double _en_13=0,_en_11=0,_en_m11=0,_en_111=0,_en_22=0,_en_211=0,_en_m211=0,_en_2212=0,_en_2112=0;
  std::vector<double> diff_vec;

  edm::Service<geo::Geometry> geom;  
  if(clusters.size()!=0){
    for(unsigned int plane=0;plane<geom->Nplanes();++plane){
      
      for(unsigned int j=0; j<clusters.size();++j) {
       
	//	std::cout<<"I AM ON PLANE #"<<plane<<std::endl;
	edm::Ptr<const recob::Cluster> clusTmp(clusterListHandle,j);
	edm::PtrVector<recob::Hit> _hits (clusTmp->Hits(plane,-1));
	//	std::cout<<"in the "<<j<<" cluster loop, hits' size= "<<_hits.size()<<"****************************************"<<std::endl;
		if(_hits.size()!=0){ //need this b/c of plane
	  
	  for(unsigned int i = 0; i < _hits.size(); ++i) {
	    //      // geo->ChannelToWire(hits[i]->Wire()->RawDigit()->Channel(), p, w);
	    //  std::cout<<"channel: "<<_hits[i]->Wire()->RawDigit()->Channel()<<"  time= "<<(_hits[i]->StartTime()+_hits[i]->EndTime())/2.<<" X time= "<<_hits[i]-> CrossingTime()<<std::endl;
	    
	    double XTime=_hits[i]-> CrossingTime();
	    const sim::SimDigit* simdigit = static_cast<const sim::SimDigit*>(_hits[i]->Wire()->RawDigit().get());
	    int numberOfElectrons = simdigit->NumberOfElectrons();
	    // std::cout<<"# of elec: "<<numberOfElectrons<<"  ";
	    if(simdigit==0){std::cout<<"simdigit=0 !!!!!!!!!!"<<std::endl;}


	    for ( int i = 0; i != numberOfElectrons; ++i )
	      {
		//std::cout<<"in the loop"<<std::endl;
		_electrons = simdigit->GetElectrons(i);
		double ArrivalTime=(_electrons->ArrivalT())/200;
		double diff=XTime-ArrivalTime;
	   
		//	std::cout<<"e's ArrivalT = "<<ArrivalTime<<" diff= "<<diff<<std::endl;
		diff_vec.push_back(diff);

		/* Below should be generalized for N planes. EC, 5-Oct-2010. */

		if(plane==0)  { 
		  if((diff<22)&&(diff>13))
		  {
		    electrons = simdigit->GetElectrons(i);
		    // double _ArrivalTime=(electrons->ArrivalT())/200;
		    // double _diff=XTime-_ArrivalTime;
		    // std::cout<<"PLANE 0,diff= "<<_diff<<std::endl;
		  }
		}
		if(plane==1)  { 
		  if((diff<36)&&(diff>27))
		    {
		      electrons = simdigit->GetElectrons(i);
		      //double _ArrivalTime=(electrons->ArrivalT())/200;
		      //double _diff=XTime-_ArrivalTime;
		      // std::cout<<"PLANE 1,diff= "<<_diff<<std::endl;
		    }
		}


	      }//for
	    //   std::cout<<"MIN of DIFF= "<<*min_element(diff_vec.begin(),diff_vec.end())<<" MAX of DIFF= "<<*max_element(diff_vec.begin(),diff_vec.end())<<std::endl;
	    diff_vec.clear();
	    if(electrons!=0){
	      const sim::LArVoxelID voxelIDval = electrons->Voxel()->VoxelID(); 
	      const sim::LArVoxelID* voxelID = &voxelIDval;
	      if(voxelID==0){std::cout<<"voxelID =0!!!!!!!!! ************"<<std::endl;}
	      // sim::LArVoxelData& voxelData =const_cast<sim::LArVoxelList*> (voxelList)->at(*voxelID);
	      const sim::LArVoxelData& voxelData = voxelList->at(*voxelID);
	    
	    
	      int numberParticles = voxelData.NumberParticles();
	      //std::cout<<"numberParticles= "<<numberParticles<<std::endl;
	      for ( int i = 0; i != numberParticles; ++i )
		{
		  int trackID = voxelData.TrackID(i);
		  //	std::cout<<"trackID= "<<trackID<<std::endl;
		  //
		  double energy=voxelData.Energy(i);
		  // std::cout<<"energy= "<<energy<<std::endl;
		  // fEnergy->Fill(energy);
		  
		  vec_trackid.push_back(trackID);
		  for( unsigned int i=0; i<_particleList.size(); ++i )
		    {
		      // const sim::ParticleList* particleList = _particleList[i];
		      particleList = _particleList[i];
		      

		      // 	double energyTrackID=voxelData.Energy[trackID];
//  	std::cout<<"ENERGY OF PRIMARY TRACKID= "<<energyTrackID<<std::endl;
		    
		      const sim::Particle* particle = particleList->at( trackID );
		      int pdg = particle->PdgCode();
		      
		      double energy2=voxelData.Energy(i);
		      // std::cout<<"energy2= "<<energy2<<std::endl;

		      // std::cout<<"part eng= "<<particle->E()<<std::endl;
		      if(pdg==13){_hit_13++;
			_en_13+=energy2;}
		      if(pdg==11){_hit_11++;
			_en_11+=energy2;
			// std::cout<<"in clus: _en_11="<<_en_11<<std::endl;
		      }
		      if(pdg==-11){_hit_m_11++;
			_en_m11+=energy2;}
		      if(pdg==111){_hit_111++;
			_en_111+=energy2;}
		      if(pdg==22){_hit_22++;
			_en_22+=energy2;}
		      if(pdg==211){_hit_211++;
			_en_211+=energy2;}
		      if(pdg==-211){_hit_m211++;
			_en_m211+=energy2;}
		      if(pdg==2212){_hit_2212++;
			_en_2212+=energy2;}
		      if(pdg==2112){_hit_2112++;
			_en_2112+=energy2;}
		    
		      //std::cout<<"True PDG= "<<pdg<<std::endl;
		      vec_pdg.push_back(pdg);
		      // std::cout<<"_en_11= "<<_en_11<<std::endl;
		    //while particle is not a primary particle and going up in a chain of trackIDs is not going to change its pdg code, go up the chain.
		      while ( (! particleList->IsPrimary( trackID )) && (((particleList->at(particle->Mother()))->PdgCode())==pdg))
			{
			  trackID = particle->Mother();
			  //std::cout<<"((NOt a PRIMARY ORIGINALLY!!! ) trackID= "<<trackID<<std::endl;
			  particle = particleList->at( trackID );
			  pdg= particle->PdgCode();
			  //	 std::cout<<"(NOt a PRIMARY ORIGINALLY!!! ) The PDG from HIT is: "<<pdg<<std::endl;
			  
			}
		    
		      // std::cout<<"The PDG from HIT is: "<<pdg<<std::endl;
		      //std::cout<<"after mother trackid= "<<trackID<<std::endl;
		      vec_trackid_mother.push_back(trackID);
		      if(energy>(7e-5)){ vec_trackid_mother_en.push_back(trackID);}
		    }
		
		}
	    }//if electrons!=0
	    ////////////////////////////////////////////
	    
	    // int numberPrimaries = particleList->NumberOfPrimaries();
// 	     std::cout<<"no of PRIMARIES: "<<numberPrimaries<<std::endl;
// 	    for ( int i = 0; i != numberPrimaries; ++i )
// 	    {
// 	    const sim::Particle* primaryParticle = particleList->Primary(i);
// 			int trackID = primaryParticle->TrackId();
// 		std::cout<<"from PRIMARY trackID= "<<trackID<<std::endl;
// 	double energy=voxelData.Energy[trackID];
// 	std::cout<<"ENERGY OF PRIMARY TRACKID= "<<energy<<std::endl;
// 	    }
	    
	    //////////////////////////////////////////////

	    _electrons=0;
	    electrons=0;
	   
	  }//for hits
	
	  //  std::cout<<"vec_pdg("<<vec_pdg.size()<<")= " ;
	  for(unsigned int i=0;i<vec_pdg.size();++i){
	    
	    // std::cout<<vec_pdg[i]<<" ";
	    
	  }
	  //std::cout<<std::endl;
	  //  std::cout<<"vec_trackid("<<vec_trackid.size()<<")= ";
// 	  for(unsigned int i=0;i<vec_trackid.size();++i){
	    
// 	     std::cout<<vec_trackid[i]<<" ";
	    
// 	  }
	  
	  // std::cout<<"vec_trackid_mother("<<vec_trackid_mother.size()<<")= ";
// 	  for(unsigned int i=0;i<vec_trackid_mother.size();++i){
	    
// 	     std::cout<<vec_trackid_mother[i]<<" ";
	    
// 	  }
	  
	  
	  it=find(vec_pdg.begin(),vec_pdg.end(),13);
	  if(it!=vec_pdg.end()){
	    // std::cout<<"matched found at position="<<int(it-vec_pdg.begin())<<std::endl;
	    no_cl_for_muon++;	  
	  }
	  else{
	    // std::cout<<"no match!"<<std::endl;
	  }
	  
	  it2=find(vec_pdg.begin(),vec_pdg.end(),11);
	  if(it2!=vec_pdg.end()){
	    no_cl_for_electron++;
	  }
	  
	  it3=find(vec_pdg.begin(),vec_pdg.end(),-11);
	  if(it3!=vec_pdg.end()){
	    no_cl_for_positron++;
	  }
	  
	  it4=find(vec_pdg.begin(),vec_pdg.end(),111);
	  if(it4!=vec_pdg.end()){
	    no_cl_for_pion_111++;
	  }
	  it6=find(vec_pdg.begin(),vec_pdg.end(),211);
	  if(it6!=vec_pdg.end()){
	    no_cl_for_pion_211++;
	  }
	  it7=find(vec_pdg.begin(),vec_pdg.end(),-211);
	  if(it7!=vec_pdg.end()){
	    no_cl_for_pion_m211++;
	  }
	  it8=find(vec_pdg.begin(),vec_pdg.end(),2212);
	  if(it8!=vec_pdg.end()){
	    no_cl_for_proton++;
	  }
	  
	  // std::cout<<std::endl;
	  //std::cout<<"numberParticles= "<<numberParticles<<std::endl;
	  //  std::cout<<"size of vec_pdg= "<<vec_pdg.size()<<std::endl;
	  sort( vec_pdg.begin(), vec_pdg.end() );
	  vec_pdg.erase( unique( vec_pdg.begin(), vec_pdg.end() ), vec_pdg.end() );
	  // std::cout<<" NO OF PARTICLES IN THIS CLUSTER IS: "<<vec_pdg.size()<<std::endl;
	   //  std::cout<<"They are: ";
// 	  for(unsigned int i=0;i<vec_pdg.size();++i){
//  	     std::cout<<vec_pdg[i]<<" ";
//  	  }
//  	  std::cout<<std::endl;
	  
	  //same for vec_trackid:
	  
	  
	  sort( vec_trackid.begin(), vec_trackid.end() );
	  vec_trackid.erase( unique( vec_trackid.begin(), vec_trackid.end() ), vec_trackid.end() );
	  //  std::cout<<" NO OF DIFFERENT TRACKIDS IN THIS CLUSTER IS: "<<vec_trackid.size()<<std::endl;
	  // std::cout<<"They are: ";
	  for(unsigned int i=0;i<vec_trackid.size();++i){
	    //   std::cout<<vec_trackid[i]<<" ";
	    all_trackids.push_back(vec_trackid[i]);
	    //mytracklist.push_back(vec_trackid[i]);
	  }
	  // std::cout<<std::endl;
	  //...............................................................
	  //Also Make vec_trackid_mother unique:
	  
	  sort( vec_trackid_mother.begin(), vec_trackid_mother.end() );
	  vec_trackid_mother.erase( unique( vec_trackid_mother.begin(), vec_trackid_mother.end() ), vec_trackid_mother.end() );
	  // std::cout<<" NO OF DIFFERENT TRACKIDS_MOTHER IN THIS CLUSTER IS: "<<vec_trackid_mother.size()<<std::endl;
	  // std::cout<<"They are: ";
	  for(unsigned int i=0;i<vec_trackid_mother.size();++i){
	    // std::cout<<vec_trackid_mother[i]<<" ";
	    
	  }
	  // std::cout<<std::endl;
	  
	  //........................................................................
	  
	  //Also Make vec_trackid_mother_en unique:

	  sort( vec_trackid_mother_en.begin(), vec_trackid_mother_en.end() );
	  vec_trackid_mother_en.erase( unique( vec_trackid_mother_en.begin(), vec_trackid_mother_en.end() ), vec_trackid_mother_en.end() );
	  // std::cout<<" NO OF DIFFERENT TRACKIDS_MOTHER_en IN THIS CLUSTER IS: "<<vec_trackid_mother_en.size()<<std::endl;
		  // std::cout<<"They are: ";
	  for(unsigned int i=0;i<vec_trackid_mother_en.size();++i){
	    // std::cout<<vec_trackid_mother_en[i]<<" ";
	    
	  }
	  //  std::cout<<std::endl;
	  
	  //........................................................................
	  
	  
	  // Q: How many clusters it takes to contain a certain particle?
	  
	  for(unsigned int i=0;i<mc_trackids.size();++i){
	    it5=find(vec_trackid.begin(),vec_trackid.end(),mc_trackids[i]);
	    if(it5!=vec_trackid.end()){
	      // std::cout<<"found match for: "<<mc_trackids[i]<<" at position: "<<it5-vec_trackid.begin()<<std::endl;
	      ids.push_back(mc_trackids[i]);//then make it unique 
	      noCluster++;
	    }
	  }
	  
	  // std::cout<<"noCluster = "<< noCluster<<std::endl;
	 
	  

	  
	  fNoParticles_pdg->Fill(vec_pdg.size());
	  // std::cout<<"LOOK ->> vec_trackid.size()= "<<vec_trackid.size()<<std::endl;
	  fNoParticles_trackid->Fill(vec_trackid.size());
	  fNoParticles_trackid_mother->Fill(vec_trackid_mother.size());
	  no_of_clusters++;
	  // std::cout<<"MUON IS CONTAINED IN "<<no_cl_for_muon<<" CLUSTERS"<<std::endl;
	  
	  no_of_particles_in_cluster+=vec_pdg.size();
	  sum_vec_trackid+=vec_trackid_mother.size();
	  total_no_hits_in_clusters+=_hits.size();
	  
	  vec_pdg.clear();
	  vec_trackid.clear();
	  vec_trackid_mother.clear();
	  vec_trackid_mother_en.clear();
	  	}//non-zero hits
	
      }//for each cluster
      
      // std::cout<<"sum_vec_trackid= "<<sum_vec_trackid<<std::endl;
      
      
sort( all_trackids.begin(), all_trackids.end() );
 all_trackids.erase( unique( all_trackids.begin(), all_trackids.end() ), all_trackids.end() );
 //	  std::cout<<" NO OF DIFFERENT TRACKIDS IN THIS EVENT IS: "<<all_trackids.size()<<std::endl;
 // std::cout<<"They are: ";
 for(unsigned int i=0;i<all_trackids.size();++i){
   // std::cout<<all_trackids[i]<<" ";
   
 }
 // std::cout<<std::endl;
 
 //now I have a vector(all_trackids) that only contains unique trackids
 //go to it and search for each trackid in every cluster
 
 
 
 sort(ids.begin(),ids.end() );
 ids.erase( unique( ids.begin(), ids.end() ), ids.end() );
 // std::cout<<"**NEW**  NO OF USED TRACKIDS IN THIS EVENT IS: "<<ids.size()<<std::endl;
 double no_of_clusters_per_track=noCluster/ids.size();
 //  std::cout<<"alright so i have no_of_clusters_per_track= "<<no_of_clusters_per_track<<std::endl;
 
 
 
 //Filling Histograms:
 
 if(no_cl_for_muon!=0){
   fCl_for_Muon->Fill(no_cl_for_muon);}
//  if(no_cl_for_electron!=0){
//    fCl_for_Electron->Fill(no_cl_for_electron);}
//  if(no_cl_for_positron!=0){
//    fCl_for_Positron->Fill(no_cl_for_positron);}
//  if(no_cl_for_pion_111!=0){
//    fCl_for_Pion_111->Fill(no_cl_for_pion_111);}
//  if(no_cl_for_pion_211!=0){
//    fCl_for_Pion_211->Fill(no_cl_for_pion_211);}
//  if(no_cl_for_pion_m211!=0){
//    fCl_for_Pion_m211->Fill(no_cl_for_pion_m211);}
//  if(no_cl_for_proton!=0){
//    fCl_for_Proton->Fill(no_cl_for_proton);}
 
 fno_of_clusters_per_track->Fill(no_of_clusters_per_track);
 
 no_cl_for_muon=0;
 no_cl_for_electron=0;
 no_cl_for_positron=0;
 no_cl_for_pion_111=0;
 no_cl_for_pion_211=0;
 no_cl_for_pion_m211=0;
 no_cl_for_proton=0;
 noCluster=0;
 //   std::cout<<"*************************************************"<<std::endl;
 //  std::cout<<"*************************************************"<<std::endl;
 // std::cout<<"                  NEW PLANE                      "<<std::endl;
 //  std::cout<<"*************************************************"<<std::endl;
 //  std::cout<<"*************************************************"<<std::endl;
 ids.clear();
 
    }//for each plane
    //std::cout<<"no_of_particles_in_cluster= "<<no_of_particles_in_cluster<<std::endl;
    //std::cout<<" no_of_clusters (with hits that are non-zero)= "<< no_of_clusters<<std::endl;
    double result=no_of_particles_in_cluster/no_of_clusters;
    // std::cout<<"FINALLY NO OF PARTICLES PER CLUSTER IS: "<<result<<std::endl;
    // std::cout<<"Sum of all hits in clusters= "<<total_no_hits_in_clusters<<std::endl;
    double no_noise_hits=hits.size()-total_no_hits_in_clusters;
    // std::cout<<"No of hits marked as noise= "<< no_noise_hits<<std::endl;
    double percent_noise=double(no_noise_hits/hits.size())*100;
    //  std::cout<<"%%%%%%%%% NOISE IS: "<< percent_noise<<std::endl;
    

    
    double no_trackid_per_cl_per_event  = sum_vec_trackid/no_of_clusters;
    // std::cout<<"sum_vec_trackid= "<<sum_vec_trackid<<" no_of_clusters: "<<no_of_clusters<<" no_trackid_per_cl_per_event= "<<no_trackid_per_cl_per_event<<std::endl;
    fNoParticles_trackid_per_event->Fill(no_trackid_per_cl_per_event);
    fNoParticles_pdg_per_event->Fill(result);
    fNoClustersInEvent->Fill(clusters.size());
    
    fPercentNoise->Fill(percent_noise);
    
    //????????here clean sum_vec_trackid

    /////////////////////////////////////////////////////////////////
      //-----------------------------------------------------------------
    
      
      
      
      
      //-----------------------------------------------------------------------
       for( unsigned int i = 0; i < mclist.size(); ++i ){
	 //edm::Ptr<const simb::MCTruth> mc(mctruthListHandle,i);
	 edm::Ptr<simb::MCTruth> mc(mclist[i]);
	 // simb::MCTruth mcp(*(mclist[i]));// We don't ever use this. EC, 6-Oct-2010.
	// std::cout<<"Number of particles is : "<<mc->NParticles()<<std::endl;
// 	std::cout<<" mc->NParticles()="<< mc->NParticles()<<std::endl;
//  TParticle part(mc->GetParticle(1));
//  std::cout<<"part.Gte<PDG()->PdgCode()= "<<part.GetPDG()->PdgCode()<<std::endl;

	for(int i = 0; i < mc->NParticles(); ++i){
	  simb::MCParticle part(mc->GetParticle(i));
	  std::cout<<"FROM MC TRUTH,the particle's pdg code is: "<<part.PdgCode()<<std::endl;
 std::cout<<"with energy= "<<part.E();
	  if(abs(part.PdgCode()) == 13){std::cout<<" I have a muon!!!"<<std::endl;
	    // std::cout<<"with energy= "<<part.Energy();
}
	  if(abs(part.PdgCode()) == 111){std::cout<<" I have a pi zero!!!"<<std::endl;}
	  //std::cout<<"here0"<<std::endl;
	}
	//std::cout<<"here1"<<std::endl;
      }
  //std::cout<<"here2"<<std::endl;
      //-----------------------------------------------------------------------
      
      
      
      
      
  }

  //  std::cout<<std::endl;
  //   std::cout<<" ********************************************************"<<std::endl;
  // std::cout<<" *************   HITS ONLY  *****************************"<<std::endl;
  //     std::cout<<" 88888888888888888888888888888888888888888888888888888888"<<std::endl;
  //......................................................................
  //......................................................................
  //......................................................................
  // ***        *****************       ***         *****************
     
  //      Now I can load just hits (not clusters) and see how many of what kind are lost due to dbscan marking them as noise:
     
     //Load hits and copy most of the code from above:
     
     // std::cout<<"Take care of "<<hits.size()<<" hits"<<std::endl;
  // std::cout<<"_en_11= "<<_en_11<<" _en_13= "<<_en_13<<std::endl;
  double hit_13=0,hit_11=0,hit_m_11=0,hit_111=0,hit_22=0,hit_211=0,hit_m211=0,hit_2212=0,hit_2112=0;
  double en_13=0,en_11=0,en_m11=0,en_111=0,en_22=0,en_211=0,en_m211=0,en_2212=0,en_2112=0;
     int no_hits=0;
     unsigned int plane_k=0;
     double total_eng_hits_p0=0;
     double total_eng_hits_p1=0;
     unsigned int w=0;
     _electrons=0;
     electrons=0;
     // if(hits.size()!=0){
     // for( int plane_k=0;plane_k<2;++plane_k){
     // std::cout<<"plane: "<<plane_k<<"  ";
       for(unsigned int j = 0; j < hits.size(); ++j) {
	 
	 no_hits++;
	 geom->ChannelToWire(hits[j]->Wire()->RawDigit()->Channel(), plane_k, w);

	 //std::cout<<"0:TOTAL ENERGY FROM HITS for P=0 = "<<total_eng_hits_p0<<std::endl;
	 //std::cout<<"0:TOTAL ENERGY FROM HITS for P=1 = "<<total_eng_hits_p1<<std::endl;


	 // std::cout<<"channel: "<<hits[j]->Wire()->RawDigit()->Channel()<<"  time= "<<(hits[j]->StartTime()+hits[j]->EndTime())/2.<<" w= "<<w<<" plane= "<<plane_k<<std::endl;
	   double XTime=hits[j]-> CrossingTime();
	   //if(XTime >1650){std::cout<<"possible fake hit line ***********"<<std::endl;}
	   const sim::SimDigit* simdigit = static_cast<const sim::SimDigit*>(hits[j]->Wire()->RawDigit().get());
	   int numberOfElectrons = simdigit->NumberOfElectrons();
	   //  std::cout<<"Hits only, numberOfElectrons= "<< numberOfElectrons<<std::endl;
	   if(numberOfElectrons==0){std::cout<<"  ZERO ELEC!!!"<<std::endl;}
	   //std::cout<<"# of elec: "<<numberOfElectrons<<"  ";
	   if(simdigit==0){std::cout<<"simdigit=0 !!!!!!!!!!"<<std::endl;}
	   //  std::cout<<"simdigit is: "<<simdigit<<std::endl;
	   for ( int i = 0; i != numberOfElectrons; ++i )
	     {
	       
	       _electrons = simdigit->GetElectrons(i);
	       double ArrivalTime=(_electrons->ArrivalT())/200;
	       double diff=XTime-ArrivalTime;
	       //	std::cout<<"e's ArrivalT = "<<ArrivalTime<<" diff= "<<diff<<std::endl;
	       diff_vec.push_back(diff);
	       if(plane_k==0)  { 
		 if((diff<22)&&(diff>13))
		   {
		     electrons = simdigit->GetElectrons(i);
		     // double _ArrivalTime=(electrons->ArrivalT())/200;
		     // double _diff=XTime-_ArrivalTime;
		     
		     
		     // std::cout<<"PLANE 0,diff= "<<_diff<<std::endl;
		   }
	       }
	       if(plane_k==1)  { 
		 if((diff<36)&&(diff>27))
		   {
		     electrons = simdigit->GetElectrons(i);
		     //double _ArrivalTime=(electrons->ArrivalT())/200;
		     //double _diff=XTime-_ArrivalTime;
		     
		     
		   }
	       }
	       
	       
	     }//numberofElectrons
	   
	   
	   //	std::cout<<"MIN of DIFF= "<<*min_element(diff_vec.begin(),diff_vec.end())<<" MAX of DIFF= "<<*max_element(diff_vec.begin(),diff_vec.end())<<std::endl;
	   diff_vec.clear();
	   
	   if(electrons!=0){
	     const sim::LArVoxelID voxelIDval = electrons->Voxel()->VoxelID(); 
	     const sim::LArVoxelID* voxelID = &voxelIDval;
	     if(voxelID==0){std::cout<<"voxelID =0!!!!!!!!! ************"<<std::endl;}
	     // sim::LArVoxelData& voxelData =const_cast<sim::LArVoxelList*> (voxelList)->at(*voxelID);
	     const sim::LArVoxelData& voxelData = voxelList->at(*voxelID);
	     
	     
	     int numberParticles = voxelData.NumberParticles();
	     //  std::cout<<"numberParticles= "<<numberParticles<<std::endl;
	     for ( int i = 0; i != numberParticles; ++i )
	       {
		 int trackID = voxelData.TrackID(i);
		 //	 std::cout<<"trackid= "<<trackID<<std::endl;
		 
		 for( unsigned int i=0; i<_particleList.size(); ++i )
		   {
		     edm::Ptr<sim::ParticleList> particleList(partListHandle,i);
		     // This barfs, EC, 6-Oct-2010.
		     const sim::Particle* particle = particleList->at( trackID );
		     int pdg = particle->PdgCode();
		     //  std::cout<<"pdg= "<<pdg<<std::endl;
		     
		     double energy3=voxelData.Energy(i);
		     //std::cout<<"plane= "<<plane_k<<std::endl;
		     if(plane_k==0){total_eng_hits_p0+=energy3;}
		     if(plane_k==1){total_eng_hits_p1+=energy3;}


		     if(pdg==13){hit_13++;
		       en_13+=energy3;}
		     if(pdg==11){hit_11++;
		       en_11+=energy3;
		       //  std::cout<<"in hits: en_11="<<en_11<<std::endl;
}
		     if(pdg==-11){hit_m_11++;
		       en_m11+=energy3;}
		     if(pdg==111){hit_111++;
		       en_111+=energy3;}
		     if(pdg==22){hit_22++;
		       en_22+=energy3;}
		     if(pdg==211){hit_211++;
		       en_211+=energy3;}
		     if(pdg==-211){hit_m211++;
		       en_m211+=energy3;}
		     if(pdg==2212){hit_2212++;
		       en_2212+=energy3;}
		     if(pdg==2112){hit_2112++;
		       en_2112+=energy3;}
		     if(pdg !=22 && pdg!=111 && pdg!=-11 && pdg !=11 && pdg!=13 && pdg!=211 && pdg!=-211 && pdg!=2212 && pdg!=2112){
		       
		       std::cout<<"SOMETHING ELSE!!! PDG= "<<pdg<<std::endl;
		     }
		     
		     //  std::cout<<"True PDG= "<<pdg<<std::endl;
		     
		   }
	       }
	     // std::cout<<"hit_13= "<<hit_13<<"  "<<"hit_11= "<<hit_11<<"  "<<"hit_m_11= "<<hit_m_11<<"  "<<"hit_111= "<<hit_111<<"  "<<"hit_22= "<<hit_22<<"  ";
	     
	     // int sum=hit_13+hit_11+hit_m_11+hit_111+hit_22;
	     // std::cout<<"sum= "<<sum<<" no_hits= "<<no_hits<<" DIFF= "<<sum-no_hits<<std::endl;
	     
	     
	     
	   }//non-zero elec
	   _electrons=0;
	   electrons=0;
	   // std::cout<<"PLANE_K= "<<plane_k<<std::endl;
       }//hits


       std::cout<<"TOTAL ENERGY FROM HITS for P=0 = "<<total_eng_hits_p0<<std::endl;
       std::cout<<"TOTAL ENERGY FROM HITS for P=1 = "<<total_eng_hits_p1<<std::endl;




	// std::cout<<"After hits,PLANE_K= "<<plane_k<<std::endl;
	 //  }//plane
       // }//non-zero hits
       // std::cout<<"no_hits= "<<no_hits<<std::endl;
      // 	std::cout<<"TOTAL # of muon hits= "<<hit_13<<std::endl;
      // 	std::cout<<"TOTAL # of electron hits= "<<hit_11<<std::endl;
      // 	std::cout<<"TOTAL # of positron hits= "<<hit_m_11<<std::endl;
      // 	std::cout<<"TOTAL # of 111 hits= "<<hit_111<<std::endl;
      // 	std::cout<<"TOTAL # of 211 hits= "<<hit_211<<std::endl;
      // 	std::cout<<"TOTAL # of -211 hits= "<<hit_m211<<std::endl;
      // 	std::cout<<"TOTAL # of 2212 hits= "<<hit_2212<<std::endl;
      // 	std::cout<<"IN CLUSTERS # of muon hits= "<<_hit_13<<std::endl;
      // 	std::cout<<"IN CLUSTERS # of electron hits= "<<_hit_11<<std::endl;
      // 	std::cout<<"IN CLUSTERS # of positron hits= "<<_hit_m_11<<std::endl;
      // 	std::cout<<"IN CLUSTERS # of 111 hits= "<<_hit_111<<std::endl;
      // 	std::cout<<"IN CLUSTERS # of 211 hits= "<<_hit_211<<std::endl;
      // 	std::cout<<"IN CLUSTERS # of -211 hits= "<<_hit_m211<<std::endl;
      // 	std::cout<<"IN CLUSTERS # of 2212 hits= "<<_hit_2212<<std::endl;
      
      // 	std::cout<<"WE MISSED % of muon hits= "<<100-((_hit_13/hit_13)*100)<<"%"<<std::endl;
      
       //	std::cout<<"WE MISSED % of electron hits= "<<100-((_hit_11/hit_11)*100)<<"%"<<std::endl;
      // 	std::cout<<"WE MISSED % of positron hits= "<<100-((_hit_m_11/hit_m_11)*100)<<"%"<<std::endl;
      // 	std::cout<<"WE MISSED % of 111 hits= "<<100-((_hit_111/hit_111)*100)<<"%"<<std::endl;
      // 	std::cout<<"WE MISSED % of 211 hits= "<<100-((_hit_211/hit_211)*100)<<"%"<<std::endl;
      // 	std::cout<<"WE MISSED % of -211 hits= "<<100-((_hit_m211/hit_m211)*100)<<"%"<<std::endl;
      // 	std::cout<<"WE MISSED % of 2212 hits= "<<100-((_hit_2212/hit_2212)*100)<<"%"<<std::endl;
      
      if(hit_13!=0){	fPercent_lost_muon_hits->Fill(100-((_hit_13/hit_13)*100));}
      if(hit_11!=0){	fPercent_lost_electron_hits->Fill(100-((_hit_11/hit_11)*100));}
      if(hit_m_11!=0){	fPercent_lost_positron_hits->Fill(100-((_hit_m_11/hit_m_11)*100));}
      if(hit_111!=0){	fPercent_lost_111_hits->Fill(100-((_hit_111/hit_111)*100));}
      
      if(hit_211!=0){	fPercent_lost_211_hits->Fill(100-((_hit_211/hit_211)*100));}
      
      if(hit_m211!=0){ fPercent_lost_m211_hits->Fill(100-((_hit_m211/hit_m211)*100));}
      if(hit_2212!=0){	fPercent_lost_2212_hits->Fill(100-((_hit_2212/hit_2212)*100));}
      if(hit_2112!=0){	fPercent_lost_2112_hits->Fill(100-((_hit_2112/hit_2112)*100));}

     //  std::cout<<"*** _en_11= "<<_en_11<<" en_11= "<<en_11<<std::endl;
// 	std::cout<<"WE MISSED % of muon energy= "<<100-((_en_13/en_13)*100)<<"%"<<std::endl;
      
//        	std::cout<<"WE MISSED % of electron energy= "<<100-((_en_11/en_11)*100)<<"%"<<std::endl;
//        	std::cout<<"WE MISSED % of positron energy= "<<100-((_en_m11/en_m11)*100)<<"%"<<std::endl;
//       	std::cout<<"WE MISSED % of 111 energy= "<<100-((_en_111/en_111)*100)<<"%"<<std::endl;
//       	std::cout<<"WE MISSED % of 211 energy= "<<100-((_en_211/en_211)*100)<<"%"<<std::endl;
//       	std::cout<<"WE MISSED % of -211 energy= "<<100-((_en_m211/en_m211)*100)<<"%"<<std::endl;
//       	std::cout<<"WE MISSED % of 2212 energy= "<<100-((_en_2212/en_2212)*100)<<"%"<<std::endl;
// 	std::cout<<"WE MISSED % of 2112 energy= "<<100-((_en_2112/en_2112)*100)<<"%"<<std::endl;
      // if(en_13==0){std::cout<<"NO MU IN THIS EVENT (en) $$$$$$$$$$$$$$$$$"<<std::endl;}
      //if(hit_13==0){std::cout<<"NO MU IN THIS EVENT (hit)$$$$$$$$$$$$$$$$$"<<std::endl;}
	if(en_13!=0){	fPercent_lost_muon_energy->Fill(100-((_en_13/en_13)*100));}
	if(en_11!=0){	fPercent_lost_electron_energy->Fill(100-((_en_11/en_11)*100));}
	if(en_m11!=0){	fPercent_lost_positron_energy->Fill(100-((_en_m11/en_m11)*100));
	  // std::cout<<"POSITRON E= "<<100-((_en_m11/en_m11)*100)<<std::endl;
	}
	if(en_111!=0){	fPercent_lost_111_energy->Fill(100-((_en_111/en_111)*100));}
	if(en_211!=0){	fPercent_lost_211_energy->Fill(100-((_en_211/en_211)*100));}
	if(en_m211!=0){ fPercent_lost_m211_energy->Fill(100-((_en_m211/en_m211)*100));}
	if(en_2212!=0){	fPercent_lost_2212_energy->Fill(100-((_en_2212/en_2212)*100));}
	if(en_2112!=0){	fPercent_lost_2112_energy->Fill(100-((_en_2112/en_2112)*100));}





	/////////////////////////////////////////////////////////////////////////////////
	  //-----------------------------------------------------------------------

	  //  FOR BRIAN:


	  //-------------------------------------------------------------------

	  std::cout<<"hello for Brian-->>>"<<std::endl;
  //------------------------------------------------------------  
//  std::vector<const recob::Wire*> wirelist;
//   try{
//     evt.Reco().Get(fInputFolder.c_str(),wirelist);
//   } 
//   catch(edm::Exception e){
//     std::cerr << "Error retrieving wire list, while looking for wires "
// 	      << "in hit::FFTHitFinder::Ana(),  "<< "directory : " 
// 	      << fInputFolder.c_str() << std::endl;
//     return jobc::kFailed;
//   }
// 	  _electrons=0;
// 	  electrons=0;
// 	  int  Total_Elec_p0=0;
// 	  int  Total_Elec_p1=0;
// 	  int wire;             
// 	  int pl=0;
// 	  unsigned int channel; 
// geo::Geometry * geom = geo::Geometry::Instance();
// 	  //loop through all wires:
	  
// 	  for(std::vector<const recob::Wire*>::iterator wireIter = wirelist.begin();
// 	      wireIter != wirelist.end();  wireIter++) {
	    
	    
	    
// 	    channel=(*wireIter)->RawDigit()->Channel();
	    
// 	    geom->ChannelToWire(channel,pl,wire);
	    
	    
// 	    _rawdigit2 = (*wireIter)->RawDigit();
// 	    sim::SimDigit* simdigit = dynamic_cast< sim::SimDigit*>(_rawdigit2);
// 	    int numberOfElectrons = simdigit->NumberOfElectrons();
	    
// 	    if(numberOfElectrons==0){std::cout<<"  ZERO ELEC!!!"<<std::endl;}
//     //std::cout<<"# of elec: "<<numberOfElectrons<<"  ";
// 	    if(simdigit==0){std::cout<<"simdigit=0 !!!!!!!!!!"<<std::endl;}
	   
// 	    if(pl==0)  { 
// 	      Total_Elec_p0 += numberOfElectrons;}
// 	    if(pl==1)  { 
// 	      Total_Elec_p1 += numberOfElectrons;}
	    
	    
	    
// 	  }//loop wires
	  
	  
// 	  std::cout<<"NO OF ELECTRONS, p0 = "<< Total_Elec_p0<<std::endl;
// 	  std::cout<<"NO OF ELECTRONS, p1 = "<< Total_Elec_p1<<std::endl;
	  
	  //now determine number of electrons for hits
	  
	  
	  //----------------------------------------------------------------------------------------
	  // Now, do the same thing for the found hits and compare the number of ionization electrons.  
	  
	  //----------------------------------------------------------------------------------------
	  
	  
	  _electrons=0;
	  electrons=0;
	  
	  double sum=0;
	  double sum0=0;
	  
	  
	  
	  double no_ele_p0=0;
	  double no_ele_p1=0;
	  unsigned int plane=0;
	  unsigned int w_=0;
	  double Tno_ele_p0=0;
	  double Tno_ele_p1=0;
	  
	  
	  for(unsigned int j = 0; j < hits.size(); ++j) {

	    
	    geom->ChannelToWire(hits[j]->Wire()->RawDigit()->Channel(), plane, w_);
	    
	    // std::cout<<"channel= "<<w_<<std::endl;
	    double XTime=hits[j]-> CrossingTime();
	    const sim::SimDigit* simdigit2 = static_cast<const sim::SimDigit*>(hits[j]->Wire()->RawDigit().get());
	    int numberOfElectrons = simdigit2->NumberOfElectrons();
	    
	    if(numberOfElectrons==0){std::cout<<"  ZERO ELEC!!!"<<std::endl;}
	    // std::cout<<"# of elec: "<<"for plane: "<<plane<<" is: "<<numberOfElectrons<<std::endl;
	    if(simdigit2==0){std::cout<<"simdigit=0 !!!!!!!!!!"<<std::endl;}
	    //  std::cout<<"simdigit is: "<<simdigit<<std::endl;
	    for ( int i = 0; i != numberOfElectrons; ++i )
	      {
		
		_electrons = simdigit2->GetElectrons(i);
		
		double ArrivalTime=(_electrons->ArrivalT())/200;
		
		double diff=XTime-ArrivalTime;
		//	std::cout<<"e's ArrivalT = "<<ArrivalTime<<" diff= "<<diff<<std::endl;
		
		if(plane==0)  { 
		  if((diff<22)&&(diff>13))
		    {
		      
		      electrons = simdigit2->GetElectrons(i);
		      
		      no_ele_p0= electrons-> NumElectrons();
		      // std::cout<<"p0,channel= "<<w_<<" e= "<<no_ele_p0<<" Hits "<<std::endl;
	Tno_ele_p0+= no_ele_p0;
	 sum0=Tno_ele_p0;

	 //std::cout<<"sum= "<<sum0<<std::endl;
		
		    }
		  //std::cout<<" no_ele_p0= "<< no_ele_p0<<std::endl;
		}
		
		
	
		if(plane==1)  { 
		  if((diff<36)&&(diff>27))
	    {
	      electrons = simdigit2->GetElectrons(i);
	      no_ele_p1= electrons-> NumElectrons();
	      //double _ArrivalTime=(electrons->ArrivalT())/200;
	      //double _diff=XTime-_ArrivalTime;
	      // std::cout<<"p1,channel= "<<w_<<" e= "<<no_ele_p1<<" Hits"<<std::endl;
	      	Tno_ele_p1+= no_ele_p1;
			 sum=Tno_ele_p1;

			 //std::cout<<"sum= "<<sum<<std::endl;
	    }
		  //std::cout<<" no_ele_p1= "<< no_ele_p1<<std::endl;
		}
		
	
		
	      }//numberofElectrons
	   //  std::cout<<"TOTAL for p0 is: "<< Tno_ele_p0<<std::endl;
// 	    std::cout<<"TOTAL for p1 is: "<< Tno_ele_p1<<std::endl;
	     sum=0;
	    sum0=0;
	  }//hits
	  
	  
	  std::cout<<"***TOTAL for p0 is: "<< Tno_ele_p0<<std::endl;
	  std::cout<<"***TOTAL for p1 is: "<< Tno_ele_p1<<std::endl;
	  

	  //-------------------first part FOR BRIAN done---------------------------

 //------------------------------------------------------------  
no_ele_p0=0;
	  no_ele_p1=0;
	  double Tno_ele_p0_w=0;
	  double Tno_ele_p1_w=0;
	   sum=0;
	    sum0=0;
	  std::cout<<"hi"<<std::endl;
	  _electrons=0;
	  electrons=0;
	  // int  Total_Elec_p0=0;
	  // int  Total_Elec_p1=0;
	  unsigned int wire=0;             
	  unsigned int pl=0;
	  unsigned int channel=0; 

	  //loop through all wires:
	  
	  for(edm::PtrVectorItr<recob::Wire> wireIter2 = wirelist.begin();
	      wireIter2 != wirelist.end();  wireIter2++) {
	    
	    
	    
	    channel=(*wireIter2)->RawDigit()->Channel();
	   
	    geom->ChannelToWire(channel,pl,wire);
	    // std::cout<<"channel: "<<wire<<std::endl;
	    
	    const sim::SimDigit* simdigit = static_cast<const sim::SimDigit*>((*wireIter2)->RawDigit().get());
	    int numberOfElectrons = simdigit->NumberOfElectrons();
	    
	    if(numberOfElectrons==0){std::cout<<"  ZERO ELEC!!!"<<std::endl;}
    //std::cout<<"# of elec: "<<numberOfElectrons<<"  ";
	    if(simdigit==0){std::cout<<"simdigit=0 !!!!!!!!!!"<<std::endl;}



 for ( int i = 0; i != numberOfElectrons; ++i )
	      {
		
		_electrons = simdigit->GetElectrons(i);
		
	
		if(pl==0)  { 
		 
		      
		      electrons = simdigit->GetElectrons(i);
		      
		      no_ele_p0= electrons-> NumElectrons();
		      // std::cout<<"p0,channel= "<<wire<<" e= "<<no_ele_p0<<" Wires"<<std::endl;
		    
		  //std::cout<<" no_ele_p0= "<< no_ele_p0<<std::endl;
	Tno_ele_p0_w+= no_ele_p0;
	 sum0=Tno_ele_p0;

	 //	std::cout<<"sum= "<<sum0<<std::endl;
		}
		
		
	
		//	std::cout<<"sum= "<<	Tno_ele_p0<<std::endl;



		if(pl==1)  { 
		 
	    
	      electrons = simdigit->GetElectrons(i);
	      no_ele_p1= electrons-> NumElectrons();
	      
	      //std::cout<<"p1,channel= "<<wire<<" e= "<<no_ele_p1<<" Wires"<<std::endl;
	    
		  //std::cout<<" no_ele_p1= "<< no_ele_p1<<std::endl;
	Tno_ele_p1_w+= no_ele_p1;
		 sum=Tno_ele_p1;

		 //std::cout<<"sum= "<<sum<<std::endl;
		
		}
		
	
	      }//numberofElectrons
 // std::cout<<"TOTAL for p0 is: "<< Tno_ele_p0_w<<std::endl;
 //    std::cout<<"TOTAL for p1 is: "<< Tno_ele_p1_w<<std::endl;

	    no_ele_p0=0;
	    no_ele_p1=0;





	    //--------------------------------------
	   //  if(pl==0)  { 
// 	      Total_Elec_p0 += numberOfElectrons;}
// 	    if(pl==1)  { 
// 	      Total_Elec_p1 += numberOfElectrons;}
// 	    //------------------------------------
	    
	    sum=0;
	    sum0=0;
	  }//loop wires
	  
	    std::cout<<"(wires)TOTAL for p0 is: "<< Tno_ele_p0_w<<std::endl;
	    std::cout<<"(wires)TOTAL for p1 is: "<< Tno_ele_p1_w<<std::endl;
	    fbrian_in->Fill(Tno_ele_p0_w, Tno_ele_p0 );
	    fbrian_coll->Fill(Tno_ele_p1_w, Tno_ele_p1 );

      
}

