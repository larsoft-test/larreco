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
  fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel")        ),
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

  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;

  //fNoParticles_pdg_per_event = tfs->make<TH1F>("fNoParticles_pdg_per_event","Average # of Particles per cluster for each event", 500,0 ,5);
  
}

void cluster::KingaClusterAna::analyze(const art::Event& evt)
{
  
  std::cout << "run    : " << evt.run() << std::endl;
  //std::cout << "subrun : " << evt.subRun() << std::endl; // Doesn't compile w. or w.o. id().
  std::cout << "event  : " << evt.id().event() << std::endl;
  //----------------------------------------------------------------

  /* This is basically a module for studying MC efficiency/purity. Kick out now if not MC. EC, 8-Oct-2010 */
  if (evt.isRealData()) 
    {
      std::cout<<"**** KingaClusterAna: Bailing. Don't call this module if you're not MC. "<<std::endl;
      exit (1);
    }

  art::Handle< std::vector<recob::Cluster>  > kingaListHandle;
  evt.getByLabel(fKingaModuleLabel,kingaListHandle);
  
  
  art::Handle< std::vector<raw::RawDigit>  > rdListHandle;
  evt.getByLabel(fDigitModuleLabel,rdListHandle);
 
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitsModuleLabel,hitListHandle);
  
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterFinderModuleLabel,clusterListHandle);
 

  
  //  std::cout<<"****simdigit.size()= "<<simdigit.size()<<std::endl;
  //   for(int i=0; i<simdigit.size();i++){
  //     std::cout<<"***numofele: "<<simdigit[i]->NumberOfElectrons()<<"  ";}
  
  // std::cout<<"simdigit->getelecetrons: "<<*simdigit[0]->GetElectrons(0)<<"  ";
  
  //----------------------------------------------------------------

  //----------------------------------------------------------------   

 
  
  art::PtrVector<raw::RawDigit> rawdigits;
  
  for (unsigned int ii = 0; ii <  rdListHandle->size(); ++ii)
    {
      art::Ptr<raw::RawDigit> rawdigit(rdListHandle,ii);
      rawdigits.push_back(rawdigit);
    }

 
  

  //get the sim::Particle collection from the art::Event and then use the Simulation/SimListUtils object to create a sim::ParticleList from the art::Event.  

 
  
 
 
  
  
  //---------------------------------------------------------------- 
 

  art::PtrVector<recob::Hit> hits;
  for (unsigned int ii = 0; ii <  hitListHandle->size(); ++ii)
    {
      art::Ptr<recob::Hit> hitHolder(hitListHandle,ii);
      hits.push_back(hitHolder);
    }


  art::PtrVector<recob::Cluster> clusters;
  for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
      clusters.push_back(clusterHolder);
    }


art::PtrVector<recob::Cluster> kingaclusters;
  for (unsigned int ii = 0; ii <  kingaListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> kingaHolder(kingaListHandle,ii);
      kingaclusters.push_back(kingaHolder);
    }
    
    
    
  std::cout<<"in Efficiency, kingaclusters.size()= "<<kingaclusters.size()<<std::endl;
  std::cout<<"in Efficiency, dbscanclusters.size()= "<<clusters.size()<<std::endl;
  //---------------------------------------------------------------
 //  art::PtrVector<recob::Wire> wirelist;
//   
//   for (unsigned int ii = 0; ii <  wireListHandle->size(); ++ii)
//     {
//       art::Ptr<recob::Wire> wireHolder(wireListHandle,ii);
//       
//       wirelist.push_back(wireHolder);
//       
//     }
//   
  
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
  art::Ptr<raw::RawDigit > _rawdigit;
  art::Ptr<raw::RawDigit > _rawdigit2;
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
  
  art::ServiceHandle<geo::Geometry> geom;  
  /*
    for(unsigned int i = 0; i < hits.size(); ++i) {
    std::cout<<"channel: "<<hits[i]->Wire()->RawDigit()->Channel()<<"  time= "<<(hits[i]->StartTime()+hits[i]->EndTime())/2.<<" X time= "<<hits[i]-> CrossingTime()<<std::endl;
    }
  */
  if(kingaclusters.size()!=0 && hits.size()!=0){
    for(unsigned int plane=0;plane<geom->Nplanes();++plane){
      geo::View_t view = geom->Plane(plane).View();
      //art::PtrVectorItr<recob::Cluster> clusterIter = clusters.begin();      
      for(unsigned int j=0; j<clusters.size();++j) 
	//while (clusterIter != clusters.end()) 
	{
	  //	std::cout<<"I AM ON PLANE #"<<plane<<std::endl;
	  if( kingaclusters[j]->View() == view){
	    art::PtrVector<recob::Hit> _hits; 
	    
	    //if (hits.size() <= 0) {clusterIter++; continue;}
	    //_hits = (*clusterIter)->Hits();
	    _hits=kingaclusters[j]->Hits();
	    
	    if(_hits.size()!=0){ //need this b/c of plane
	      
	      for(unsigned int i = 0; i < _hits.size(); ++i) {
		
		
		//std::cout<<"channel: "<<_hits[i]->Wire()->RawDigit()->Channel()<<"  time= "<<(_hits[i]->StartTime()+_hits[i]->EndTime())/2.<<" X time= "<<_hits[i]-> CrossingTime()<<std::endl;
		
		double XTime=_hits[i]->PeakTime();
		
		// grab the channel from this hit
		unsigned int channel = _hits[i]->Wire()->RawDigit()->Channel();
		
		// loop over the SimChannels to find this one
		art::Ptr<sim::SimChannel> sc;
		
		
	   


		
	   
	      }//for hits
	
	    
	  
	     
	  
	  
	   
	      
	  
	  
	  
	    }//non-zero hits
	  }//end if cluster is in correct view
	  //clusterIter++;
	
	}//for each cluster
      
      // std::cout<<"sum_vec_trackid= "<<sum_vec_trackid<<std::endl;
     
      
     
 
 
    

    }//for each plane
  
      
      
    //-----------------------------------------------------------------------
   
    //-----------------------------------------------------------------------
      
      
      
      
      
  }

 
}

