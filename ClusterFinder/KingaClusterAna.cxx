////////////////////////////////////////////////////////////////////////
//
// KingaCluster analyzer
//
// kinga.partyka@yale.edu
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

#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
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
#include "TTree.h"
#include "ClusterFinder/ClusterCheater.h"
#include "MCCheater/BackTracker.h"


 
//-------------------------------------------------
cluster::KingaClusterAna::KingaClusterAna(fhicl::ParameterSet const& pset) : 
  
  fLineMergerModuleLabel       (pset.get< std::string >("LineMergerModuleLabel")),
  fKingaModuleLabel         (pset.get< std::string >("KingaModuleLabel")        ),
  fEndPoint2DModuleLabel        (pset.get< std::string >("EndPoint2DModuleLabel")),
  fClusterCheaterModuleLabel        (pset.get< std::string >("ClusterCheaterModuleLabel")),
  fGenieGenModuleLabel        (pset.get< std::string >("GenieGenModuleLabel")),
  fLArGeantModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")),
  frun(0),
  fevent(0),
  ftime_vertex_true(0),
  fno_clusters_true(200),
  fno_clusters_reco(200),
  fno_clusters_linemerger(200),
  fno_primaries(200),
  fgenie_no_primaries(200)
  
{

}

//-------------------------------------------------
cluster::KingaClusterAna::~KingaClusterAna()
{
delete fwire_vertex_true;
delete fTTree_wire_vertex_reco;
delete fTTree_time_vertex_reco;
delete fclusters_planeNo_true;
delete fclusters_planeNo_reco;
delete flinemergerclusters_planeNo;
delete fStart_pt_w_true;
delete fStart_pt_t_true;
delete fStart_pt_w_reco;
delete fStart_pt_t_reco;
delete fStart_pt_w_linemerger;
delete fStart_pt_t_linemerger;
delete fcheated_cluster_size;
delete flinemerger_cluster_size;
delete fprimaries_pdg;
delete fEng;
delete fPx;
 delete fPy;
 delete fPz;
 delete fStartPointx;
 delete fStartPointy;
 delete fStartPointz;
 delete fEndPointx;
 delete fEndPointy;
 delete fEndPointz;
 delete fNumberDaughters;
 
 delete fgenie_primaries_pdg;
delete fgenie_Eng;
delete fgenie_Px;
 delete fgenie_Py;
 delete fgenie_Pz;
 delete fgenie_P;
 delete fgenie_status_code;
 delete fgenie_mass;
 delete fgenie_trackID;
 delete fgenie_ND;
 delete fgenie_mother;
}

void cluster::KingaClusterAna::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;
  
  fTree = tfs->make<TTree>("anatree","KingaCluster analysis");
 
 fwire_vertex_true= new double[2];
 fTTree_wire_vertex_reco= new double[2];
 fTTree_time_vertex_reco= new double[2];
 fStart_pt_w_true=new double[fno_clusters_true];
 fStart_pt_t_true=new double[fno_clusters_true];
 fStart_pt_w_reco=new double[fno_clusters_reco];
 fStart_pt_t_reco=new double[fno_clusters_reco];
 fStart_pt_w_linemerger=new double[fno_clusters_linemerger];
 fStart_pt_t_linemerger=new double[fno_clusters_linemerger];
 fclusters_planeNo_true=new int[fno_clusters_true];
 fclusters_planeNo_reco=new int[fno_clusters_reco];
 flinemergerclusters_planeNo=new int[fno_clusters_linemerger];
 fcheated_cluster_size=new int[fno_clusters_true];
 flinemerger_cluster_size=new int[fno_clusters_linemerger];
 
 fprimaries_pdg= new int[fno_primaries];
  fEng= new double[fno_primaries];
  fPx= new double[fno_primaries];
  fPy= new double[fno_primaries];
  fPz= new double[fno_primaries];
  fStartPointx= new double[fno_primaries];
  fStartPointy= new double[fno_primaries];
  fStartPointz= new double[fno_primaries];
  fEndPointx= new double[fno_primaries];
  fEndPointy= new double[fno_primaries];
  fEndPointz= new double[fno_primaries];
  fNumberDaughters= new int[fno_primaries];
 
 
fgenie_primaries_pdg= new double[fgenie_no_primaries];
 fgenie_Eng= new double[fgenie_no_primaries];
 fgenie_Px= new double[fgenie_no_primaries];
  fgenie_Py= new double[fgenie_no_primaries];
 fgenie_Pz= new double[fgenie_no_primaries];
 fgenie_P= new double[fgenie_no_primaries];
fgenie_status_code= new int[fgenie_no_primaries];
fgenie_mass= new double[fgenie_no_primaries];
fgenie_trackID= new int[fgenie_no_primaries];
fgenie_ND= new int[fgenie_no_primaries];
 fgenie_mother= new int[fgenie_no_primaries];
 
 
 
 fTree->Branch("run",&frun,"run/I");
 fTree->Branch("event",&fevent,"event/I");
 fTree->Branch("wire_vertex_true", fwire_vertex_true, "wire_vertex_true[2]/D");
 fTree->Branch("time_vertex_true",&ftime_vertex_true,"time_vertex_true/D");
 fTree->Branch("TTree_wire_vertex_reco", fTTree_wire_vertex_reco, "TTree_wire_vertex_reco[2]/D");
 fTree->Branch("TTree_time_vertex_reco", fTTree_time_vertex_reco, "TTree_time_vertex_reco[2]/D");
 
 
fTree->Branch("no_clusters_true",&fno_clusters_true,"no_clusters_true/I");
fTree->Branch("no_clusters_reco",&fno_clusters_reco,"no_clusters_reco/I");
fTree->Branch("no_clusters_linemerger",&fno_clusters_linemerger,"no_clusters_linemerger/I");

fTree->Branch("clusters_planeNo_true",fclusters_planeNo_true,"clusters_planeNo_true[no_clusters_true]/I");
fTree->Branch("clusters_planeNo_reco",fclusters_planeNo_reco,"clusters_planeNo_reco[no_clusters_reco]/I");
fTree->Branch("linemergerclusters_planeNo",flinemergerclusters_planeNo,"linemergerclusters_planeNo[no_clusters_linemerger]/I");

 fTree->Branch("Start_pt_w_true", fStart_pt_w_true, "Start_pt_w_true[no_clusters_true]/D");
 fTree->Branch("Start_pt_t_true", fStart_pt_t_true, "Start_pt_t_true[no_clusters_true]/D");
 fTree->Branch("Start_pt_w_reco", fStart_pt_w_reco, "Start_pt_w_reco[no_clusters_reco]/D");
 fTree->Branch("Start_pt_t_reco", fStart_pt_t_reco, "Start_pt_t_reco[no_clusters_reco]/D");
 fTree->Branch("Start_pt_w_linemerger", fStart_pt_w_linemerger, "Start_pt_w_linemerger[no_clusters_linemerger]/D");
 fTree->Branch("Start_pt_t_linemerger", fStart_pt_t_linemerger, "Start_pt_t_linemerger[no_clusters_linemerger]/D");
 
 
 
 fTree->Branch("cheated_cluster_size", fcheated_cluster_size, "cheated_cluster_size[no_clusters_true]/I"); // no of hits in each cluster
 fTree->Branch("linemerger_cluster_size", flinemerger_cluster_size, "linemerger_cluster_size[no_clusters_linemerger]/I");
 
  //......................................................
// from geant4:

  fTree->Branch("no_primaries",&fno_primaries,"no_primaries/I");
  fTree->Branch("primaries_pdg",fprimaries_pdg,"primaries_pdg[no_primaries]/I");
  fTree->Branch("Eng",fEng,"Eng[no_primaries]/D");
  fTree->Branch("Px",fPx,"Px[no_primaries]/D");
  fTree->Branch("Py",fPy,"Py[no_primaries]/D");
  fTree->Branch("Pz",fPz,"Pz[no_primaries]/D");
  fTree->Branch("StartPointx",fStartPointx,"StartPointx[no_primaries]/D");
  fTree->Branch("StartPointy",fStartPointy,"StartPointy[no_primaries]/D");
  fTree->Branch("StartPointz",fStartPointz,"StartPointz[no_primaries]/D");
  fTree->Branch("EndPointx",fEndPointx,"EndPointx[no_primaries]/D");
  fTree->Branch("EndPointy",fEndPointy,"EndPointy[no_primaries]/D");
  fTree->Branch("EndPointz",fEndPointz,"EndPointz[no_primaries]/D");
  fTree->Branch("NumberDaughters",fNumberDaughters,"NumberDaughters[no_primaries]/I");
  
  fTree->Branch("ccnc_truth",&ccnc_truth,"ccnc_truth/I");
  fTree->Branch("mode_truth",&mode_truth,"mode_truth/I");
  
  
 //..................................
 // now from genie:
  fTree->Branch("genie_no_primaries",&fgenie_no_primaries,"genie_no_primaries/I");
  fTree->Branch("genie_primaries_pdg",fgenie_primaries_pdg,"genie_primaries_pdg[genie_no_primaries]/D");
  fTree->Branch("genie_Eng",fgenie_Eng,"genie_Eng[genie_no_primaries]/D");
  fTree->Branch("genie_Px",fgenie_Px,"genie_Px[genie_no_primaries]/D");
  fTree->Branch("genie_Py",fgenie_Py,"genie_Py[genie_no_primaries]/D");
  fTree->Branch("genie_Pz",fgenie_Pz,"genie_Pz[genie_no_primaries]/D");
  fTree->Branch("genie_P",fgenie_P,"genie_P[genie_no_primaries]/D");
  fTree->Branch("genie_status_code",fgenie_status_code,"genie_status_code[genie_no_primaries]/D");
  fTree->Branch("genie_mass",fgenie_mass,"genie_mass[genie_no_primaries]/D");
  fTree->Branch("genie_trackID",fgenie_trackID,"genie_trackID[genie_no_primaries]/I");
  fTree->Branch("genie_ND",fgenie_ND,"genie_ND[genie_no_primaries]/I");
  fTree->Branch("genie_mother",fgenie_mother,"genie_mother[genie_no_primaries]/I");
  
 
 

 
 
 
 
//............................................... 
// DIFFERENCE BETWEEN KINGA CLUSTERS AND CHEATED(=TRUE) CLUSTERS
 fdiff_no_vertex_clusters_p0= tfs->make<TH1F>("fdiff_no_vertex_clusters_p0","Difference between No of Clusters around a vertex found by ClusterCheater vs KingaCluster, Induction Plane", 20,-10 ,10  );
  fdiff_no_vertex_clusters_p1= tfs->make<TH1F>("fdiff_no_vertex_clusters_p1","Difference between No of Clusters around a vertex found by ClusterCheater vs KingaCluster, Collection Plane", 20,-10 ,10  );
  
  // DIFFERENCE BETWEEN LINE MERGER CLUSTERS AND CHEATED(=TRUE) CLUSTERS
 fdiff_no_vertex_linemergerclusters_p0= tfs->make<TH1F>("fdiff_no_vertex_linemergerclusters_p0","Difference between No of Clusters around a vertex found by ClusterCheater vs LineMerger, Induction Plane", 20,-10 ,10  );
  fdiff_no_vertex_linemergerclusters_p1= tfs->make<TH1F>("fdiff_no_vertex_linemergerclusters_p1","Difference between No of Clusters around a vertex found by ClusterCheater vs LineMerger, Collection Plane", 20,-10 ,10  );
 
 
 
 //vertex truth vs reco:
 fdiff_time_vtx_p0= tfs->make<TH1F>("fdiff_time_vtx_p0","Difference between truth and reco vertex time position, Induction Plane", 20000,-1000 ,1000  );
 fdiff_wire_vtx_p0= tfs->make<TH1F>("fdiff_wire_vtx_p0","Difference between truth and reco vertex wire position, Induction Plane", 100,-50 ,50  );
 fdiff_wire_vtx_p1= tfs->make<TH1F>("fdiff_wire_vtx_p1","Difference between truth and reco vertex wire position, Collection Plane", 100,-50 ,50  );
 
 
 // No of Proton Tracks:
 fNoProtonTracks_p0_cheatedCl= tfs->make<TH1F>("fNoProtonTracks_p0_cheatedCl","No of proton tracks, Induction Plane", 15,0 ,15  );
 fNoProtonTracks_p1_cheatedCl= tfs->make<TH1F>("fNoProtonTracks_p1_cheatedCl","No of proton tracks, Collection Plane", 15,0 ,15  );
 fNoProtonTracks_p0_linemergerCl= tfs->make<TH1F>("fNoProtonTracks_p0_linemergerCl","No of proton tracks, Induction Plane", 15,0 ,15  );
 fNoProtonTracks_p1_linemergerCl= tfs->make<TH1F>("fNoProtonTracks_p1_linemergerCl","No of proton tracks, Collection Plane", 15,0 ,15  );
 fNoProtonTracks_p0_kingaCl= tfs->make<TH1F>("fNoProtonTracks_p0_kingaCl","No of proton tracks, Induction Plane", 15,0 ,15  );
  fNoProtonTracks_p1_kingaCl= tfs->make<TH1F>("fNoProtonTracks_p1_kingaCl","No of proton tracks, Induction Plane", 15,0 ,15  );
  
}

void cluster::KingaClusterAna::analyze(const art::Event& evt)
{
  std::cout<<"Hello, You are in KingaClusterAna::analyze"<<std::endl;
  std::cout << "run    : " <<evt.run()<<"event  : " << evt.id().event() << std::endl;
 frun= evt.run();
 fevent=evt.id().event();


 
  if (evt.isRealData()) 
    {
      std::cout<<"ATTENTION, THIS IS A DATA FILE, CANNOT DO COMPARISON WITH TRUTH !!!! "<<std::endl;
      return;
    }
//.....................................................................
 art::ServiceHandle<geo::Geometry> geom;
 art::ServiceHandle<util::LArProperties> larp;
 sim::ParticleList _particleList = sim::SimListUtils::GetParticleList(evt, fLArGeantModuleLabel);
 
 std::cout<<"geom->Nchannels()= "<<geom->Nchannels()<<std::endl;
    // get the sim::SimChannels
  std::vector<const sim::SimChannel*> sccol;
  evt.getView(fLArGeantModuleLabel, sccol);
  std::cout<<" ^^^^^^^^ sccol.size()= "<<sccol.size()<<std::endl;
  std::vector<const sim::SimChannel*> scs(geom->Nchannels(),0);
  for(size_t i = 0; i < sccol.size(); ++i) scs[sccol[i]->Channel()] = sccol[i];
    
    
    
    
//..................................................................



   
ftime_vertex.clear();
fwire_vertex.clear();
ftime_vertex_reco.clear();
fwire_vertex_reco.clear();

 //............MC TRUTH VERTEX:...........................................
 
 art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
 evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle);
 art::PtrVector<simb::MCTruth> mclist;
  for (unsigned int ii = 0; ii <  mctruthListHandle->size(); ++ii)
    {
      art::Ptr<simb::MCTruth> mctparticle(mctruthListHandle,ii);
      mclist.push_back(mctparticle);
    } 
    
    std::cout<<" mclist.size()= "<<mclist.size()<<std::endl;
    
    for( unsigned int i = 0; i < 1; ++i ){
    //for( unsigned int i = 0; i < mclist.size(); ++i ){

    art::Ptr<simb::MCTruth> mc(mclist[i]);

    simb::MCParticle neut(mc->GetParticle(i));

    // std::cout<<"vertex: "<<neut.Nu().Vx()<<" "<<neut.Nu().Vy()<<" "<<neut.Nu().Vz()<<std::endl;
    fMCvertex[0] =neut.Vx();
    fMCvertex[1] =neut.Vy();
    fMCvertex[2] =neut.Vz();
    
    ccnc_truth = mclist[0]->GetNeutrino().CCNC();
    mode_truth = mclist[0]->GetNeutrino().Mode();
    
    std::cout<<"MCvertex[0]= "<<fMCvertex[0]<<std::endl;
    std::cout<<"driftvelocity= "<<larp->DriftVelocity(larp->Efield(),larp->Temperature())<<std::endl;
    double presamplings=60.0;
    double drifttick=(fMCvertex[0]/larp->DriftVelocity(larp->Efield(),larp->Temperature()))*(1./.198)+presamplings;
    
    std::cout<<"%%%%%%%%%%%%%%%%%%   drifttick= "<<std::setprecision(10)<<drifttick<<std::endl;
    ftime_vertex.push_back(drifttick);
    ftime_vertex.push_back(drifttick);
    
    ftime_vertex_true=drifttick;

  }
  
  
  //......................
  for( unsigned int i = 0; i < 1; ++i ){
  //for( unsigned int i = 0; i < mclist.size(); ++i ){
    art::Ptr<simb::MCTruth> mc(mclist[i]);
    simb::MCParticle neut(mc->GetParticle(i));
    
  std::cout<<"List from Genie: "<<std::endl;
  
  fgenie_no_primaries=mc->NParticles();
  
     for(int j = 0; j < mc->NParticles(); ++j){
    simb::MCParticle part(mc->GetParticle(j));
    
    std::cout<<"pdg= "<<part.PdgCode()<<" ,Process="<<part.Process()<<" StatusCode= "<<part.StatusCode()<<" mass= "<<part.Mass()<<" p= "<<part.P()<<" E= "<<part.E()<<" trackID= "<<part.TrackId()<<" ND= "<<part.NumberDaughters()<<" Mother= "<<part.Mother()<<std::endl;
    
    
    
  fgenie_primaries_pdg[j]=part.PdgCode();
fgenie_Eng[j]=part.E();
fgenie_Px[j]=part.Px();
 fgenie_Py[j]=part.Py();
fgenie_Pz[j]=part.Pz();
fgenie_P[j]=part.Px();
 fgenie_status_code[j]=part.StatusCode();
 fgenie_mass[j]=part.Mass();
 fgenie_trackID[j]=part.TrackId();
 fgenie_ND[j]=part.NumberDaughters();
fgenie_mother[j]=part.Mother();
    
    
    
    }
    }
    //.............
  // now wire vertex:
   unsigned int channel2,plane2,wire2,tpc2; 
  for(size_t tpc = 0; tpc < geom->NTPC(); ++tpc){
   std::cout<<"No of planes = "<<geom->Nplanes(tpc)<<std::endl;
  for(unsigned int plane=0;plane<geom->Nplanes(tpc);plane++){
  if(plane==0){
	fMCvertex[0]=.3;//force time coordinate to be closer to induction plane 
	}
      else{
	fMCvertex[0]=-.3;//force time coordinate to be closer to collection plane
     }
     
  
     
  try{
  channel2 = geom->NearestChannel(fMCvertex,plane);
   
  }
  catch(cet::exception &e){
  mf::LogWarning("KingaClusterexc")<<e;
  
    
  //std::cout<<"fMCvertex[2]= "<<fMCvertex[2]<<" plane="<<plane<<" DetLength= "<<geom->DetLength()<<" geom->Nchannels()= "<<geom->Nchannels()<<std::endl;
  
  if(plane==0 && fMCvertex[2]<5) channel2=0;
  else if(plane==0 && fMCvertex[2]>geom->DetLength()-5) channel2=(geom->Nchannels())/2 -1;
  else if(plane==1 && fMCvertex[2]>geom->DetLength()-5) channel2=geom->Nchannels()-1;
  else if(plane==1 && fMCvertex[2]<5) channel2=(geom->Nchannels())/2 -1;

  
  }
      geom->ChannelToWire(channel2,tpc2,plane2,wire2);   
   std::cout<<"%%%%%%%%%%%%%%%%%%   WIRE VERTEX IS: "<<wire2<<std::endl;
   fwire_vertex.push_back(wire2);
   
   fwire_vertex_true[plane]=wire2;
   }
   }
  
   //............END OF TRUTH VERTEX.........................................
  
   //............RECO VERTEX:...........................................
   art::Handle< std::vector<recob::EndPoint2D> > endpointListHandle;
  evt.getByLabel(fEndPoint2DModuleLabel,endpointListHandle);
  art::PtrVector<recob::EndPoint2D> endpointlist;
 
    for (unsigned int i = 0; i < endpointListHandle->size(); ++i){
      art::Ptr<recob::EndPoint2D> endpointHolder(endpointListHandle,i);
      endpointlist.push_back(endpointHolder);
    }
    
  std::cout<<"SHOULD BE GETTING RECO VERTEX, endpointlist.size()= "<<endpointlist.size()<<std::endl;
  
  
  for (unsigned int j = 0; j<endpointlist.size();j++){

      std::cout<<"j="<<j<<" W_VERTEX_RECO"<<endpointlist[j]->WireNum()<<" T_VERTEX_RECO="<<endpointlist[j]->DriftTime()<<std::endl;
          
          fTTree_wire_vertex_reco[j]=endpointlist[j]->WireNum();
          fTTree_time_vertex_reco[j]=endpointlist[j]->DriftTime();
          
          ftime_vertex_reco.push_back(endpointlist[j]->DriftTime());
          fwire_vertex_reco.push_back(endpointlist[j]->WireNum());

        }
  
  
     //............END OF RECO VERTEX.........................................


    int proton_track_ind=0;
    int proton_track_coll=0;
    int proton_hit_ind=0;
    int proton_hit_coll=0;
    art::PtrVector<recob::Hit> hits;
    std::vector<int> vec_trackid;
    vec_trackid.clear();





  // get LineMergerClusters............
  std::cout<<"Trying to get line merger clusters***"<<std::endl;

  art::Handle< std::vector<recob::Cluster> > linemergerclusterListHandle;
  evt.getByLabel(fLineMergerModuleLabel,linemergerclusterListHandle);
  art::PtrVector<recob::Cluster> LineMergerClusIn;    
     
  flinemergerCl_p0=0;
    flinemergerCl_p1=0;
    flinemergerCl_near_vertex_p0=0;
    flinemergerCl_near_vertex_p1=0;
     
     
     for(unsigned int ii = 0; ii < linemergerclusterListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster(linemergerclusterListHandle, ii);
      LineMergerClusIn.push_back(cluster);
      
      
      //Fill TTree:
      flinemerger_cluster_size[ii]=cluster->Hits().size();
      
      if(cluster->View()==geo::kU){
      flinemergerclusters_planeNo[ii]=0;
      fStart_pt_w_linemerger[ii]=cluster->StartPos()[0];
      fStart_pt_t_linemerger[ii]=cluster->StartPos()[1];
      
      }
       else if(cluster->View()==geo::kV){
       flinemergerclusters_planeNo[ii]=1;
      fStart_pt_w_linemerger[ii]=cluster->StartPos()[0];
      fStart_pt_t_linemerger[ii]=cluster->StartPos()[1];
      
      }
      
      
      
      
      if(cluster->Hits().size()>2 && cluster->View()==geo::kU){
      
      flinemergerCl_p0++;
      
      //std::cout<<cluster->StartPos()[0]<<", "<<cluster->StartPos()[1]<<" --> "<<cluster->EndPos()[0]<<", "<<cluster->EndPos()[1]<<" No hits= "<<cluster->Hits().size()<<std::endl;
      
      
      if(fabs(cluster->StartPos()[0]-fwire_vertex_true[0])<6 && fabs(cluster->StartPos()[1]-ftime_vertex[0])<90 ){
      flinemergerCl_near_vertex_p0++;
      
      }
      
      
      }
      //..........
      if(cluster->Hits().size()>2 && cluster->View()==geo::kV){
      
      flinemergerCl_p1++;
      
       if(fabs(cluster->StartPos()[0]-fwire_vertex_true[1])<6 && fabs(cluster->StartPos()[1]-ftime_vertex[1])<90 ){
      flinemergerCl_near_vertex_p1++;
      
      }
      
      }
      
       //*********************************
     // find out what particle each cluster belongs to:
     
     
     hits=cluster->Hits();
     
    for(unsigned int h=0; h<hits.size(); h++){
    
    //std::cout<<"hits[h] channel= "<<hits[h]->Wire()->RawDigit()->Channel()<<std::endl;
   
    std::vector<cheat::TrackIDE> trackides = cheat::BackTracker::HitToTrackID(*scs[hits[h]->Channel()], hits[h]);
   
    std::vector<cheat::TrackIDE>::iterator idesitr = trackides.begin();
    //std::cout<<"trackides= "<<trackides.size();
    
    
   
    while( idesitr != trackides.end() ){
    
    vec_trackid.push_back((*idesitr).trackID);
    //std::cout<<"EngFraction= "<<(*idesitr).energyFrac<<std::endl;
    const sim::Particle* particle = _particleList.at( (*idesitr).trackID);
    int pdg = particle->PdgCode();
    //std::cout<<"pdg= "<<pdg<<std::endl;
    
    if(pdg==2212 || pdg==-2212){
       if(cluster->View()==geo::kU) proton_hit_ind++;
       else if(cluster->View()==geo::kV) proton_hit_coll++;
    
    }
    
    idesitr++;
    }//trackids
    
    vec_trackid.clear();
    
    }//hits
      
      if (cluster->View()==geo::kU && hits.size()>=2 && (proton_hit_ind/hits.size())>0.7){
      proton_track_ind++;
      std::cout<<" cluster #"<<ii<<" is a proton track! (IND)"<<std::endl;
      }
      else if (cluster->View()==geo::kV && hits.size()>=2 && (proton_hit_coll/hits.size())>0.7){
      proton_track_coll++;
      std::cout<<" cluster #"<<ii<<" is a proton track! (COLL)"<<std::endl;
      }
     proton_hit_ind=0;
     proton_hit_coll=0;
      
      
      
      
    }
    
    fno_clusters_linemerger=LineMergerClusIn.size();
    
std::cout<<"Total No of LINE MERGER clusters= "<<LineMergerClusIn.size()<<std::endl;
std::cout<<"for plane 0:"<<flinemergerCl_p0<<std::endl;
std::cout<<"for plane 1:"<<flinemergerCl_p1<<std::endl;
std::cout<<"Total No of LINE MERGER clusters ***NEAR THE VERTEX*** :"<<std::endl;
std::cout<<"for plane 0:"<<flinemergerCl_near_vertex_p0<<std::endl;
std::cout<<"for plane 1:"<<flinemergerCl_near_vertex_p1<<std::endl;
std::cout<<"no of proton tracks for plane 0: "<<proton_track_ind<<std::endl;
std::cout<<"no of proton tracks for plane 1: "<<proton_track_coll<<std::endl;

fNoProtonTracks_p0_linemergerCl->Fill(proton_track_ind);
fNoProtonTracks_p1_linemergerCl->Fill(proton_track_coll);
//............ End of analyzing linemerger clusters............................
  
  proton_track_ind=0;
  proton_track_coll=0;
  proton_hit_ind=0;
  proton_hit_coll=0;
  vec_trackid.clear();
  
  
  //........GET KINGACLUSTERS:
  
  art::Handle< std::vector<recob::Cluster>  > kingaListHandle;
  evt.getByLabel(fKingaModuleLabel,kingaListHandle);
   art::PtrVector<recob::Cluster> KingaClusIn;    
     
    fkingaCl_p0=0;
    fkingaCl_p1=0;
    fkingaCl_near_vertex_p0=0;
    fkingaCl_near_vertex_p1=0;
     
     
     for(unsigned int ii = 0; ii < kingaListHandle->size(); ++ii)
    {
      art::Ptr<recob::Cluster> cluster(kingaListHandle, ii);
      KingaClusIn.push_back(cluster);
      
      
      
      
      if(cluster->View()==geo::kU){
      fclusters_planeNo_reco[ii]=0;
      fStart_pt_w_reco[ii]=cluster->StartPos()[0];
      fStart_pt_t_reco[ii]=cluster->StartPos()[1];
      fkingaCl_p0++;
      //std::cout<<"p0, cluster# "<<ii<<" startPoint: "<<cluster->StartPos()[0]<<" , "<<cluster->StartPos()[1]<<std::endl;
      
       if(fabs(cluster->StartPos()[0]-fwire_vertex_true[0])<6 && fabs(cluster->StartPos()[1]-ftime_vertex[0])<90 ){
      fkingaCl_near_vertex_p0++;
      
        }
      
      }
       else if(cluster->View()==geo::kV){
       fclusters_planeNo_reco[ii]=1;
      fStart_pt_w_reco[ii]=cluster->StartPos()[0];
      fStart_pt_t_reco[ii]=cluster->StartPos()[1];
       fkingaCl_p1++;
      // std::cout<<"p1, cluster# "<<ii<<" startPoint: "<<cluster->StartPos()[0]<<" , "<<cluster->StartPos()[1]<<std::endl;
       
         if(fabs(cluster->StartPos()[0]-fwire_vertex_true[1])<6 && fabs(cluster->StartPos()[1]-ftime_vertex[1])<90 ){
         fkingaCl_near_vertex_p1++;
      
        }
      }
      
     
      
      //*********************************
     // find out what particle each cluster belongs to:
     
     
     hits=cluster->Hits();
     
    for(unsigned int h=0; h<hits.size(); h++){
    
    //std::cout<<"hits[h] channel= "<<hits[h]->Wire()->RawDigit()->Channel()<<std::endl;
   
    std::vector<cheat::TrackIDE> trackides = cheat::BackTracker::HitToTrackID(*scs[hits[h]->Channel()], hits[h]);
   
    std::vector<cheat::TrackIDE>::iterator idesitr = trackides.begin();
    //std::cout<<"trackides= "<<trackides.size();
    
    
   
    while( idesitr != trackides.end() ){
    
    vec_trackid.push_back((*idesitr).trackID);
    //std::cout<<"EngFraction= "<<(*idesitr).energyFrac<<std::endl;
    const sim::Particle* particle = _particleList.at( (*idesitr).trackID);
    int pdg = particle->PdgCode();
    //std::cout<<"pdg= "<<pdg<<std::endl;
    
    if(pdg==2212 || pdg==-2212){
       if(cluster->View()==geo::kU) proton_hit_ind++;
       else if(cluster->View()==geo::kV) proton_hit_coll++;
    
    }
    
    idesitr++;
    }//trackids
    
    vec_trackid.clear();
    
    }//hits
      
      if (cluster->View()==geo::kU && (proton_hit_ind/hits.size())>0.7){
      proton_track_ind++;
      std::cout<<" cluster #"<<ii<<" is a proton track! (IND)"<<std::endl;
      }
      else if (cluster->View()==geo::kV && (proton_hit_coll/hits.size())>0.7){
      proton_track_coll++;
      std::cout<<" cluster #"<<ii<<" is a proton track! (COLL)"<<std::endl;
      }
     proton_hit_ind=0;
     proton_hit_coll=0;
      
      
      
      
    }
    
    fno_clusters_reco=KingaClusIn.size();
    
std::cout<<"Total No of KINGA clusters= "<<KingaClusIn.size()<<std::endl;
std::cout<<"for plane 0:"<<fkingaCl_p0<<std::endl;
std::cout<<"for plane 1:"<<fkingaCl_p1<<std::endl;
std::cout<<"Total No of KINGA clusters ***NEAR THE VERTEX*** :"<<std::endl;
std::cout<<"for plane 0:"<<fkingaCl_near_vertex_p0<<std::endl;
std::cout<<"for plane 1:"<<fkingaCl_near_vertex_p1<<std::endl;
std::cout<<"no of proton tracks for plane 0: "<<proton_track_ind<<std::endl;
std::cout<<"no of proton tracks for plane 1: "<<proton_track_coll<<std::endl;


 fNoProtonTracks_p0_kingaCl->Fill(proton_track_ind);
 fNoProtonTracks_p1_kingaCl->Fill(proton_track_coll);
 
//............ End of analyzing KINGAclusters............................
  
   
  
   proton_track_ind=0;
   proton_track_coll=0;
   proton_hit_ind=0;
   proton_hit_coll=0;
   vec_trackid.clear();
 
  
  
  
  // get CheatedClusters............
std::cout<<"Trying to get cheated clusters***"<<std::endl;








    fcheatedCl_p0=0;
    fcheatedCl_p1=0;
    fcheatedCl_near_vertex_p0=0;
    fcheatedCl_near_vertex_p1=0;
    
    
    
 

 art::Handle< std::vector<recob::Cluster> > cheatedclusterListHandle;
  evt.getByLabel(fClusterCheaterModuleLabel,cheatedclusterListHandle);
  art::PtrVector<recob::Cluster> CheatedClusIn;
  
      for(unsigned int ii = 0; ii < cheatedclusterListHandle->size(); ++ii)
    {
    std::cout<<"working on cluster #"<<ii<<std::endl;
      art::Ptr<recob::Cluster> cluster(cheatedclusterListHandle, ii);
      CheatedClusIn.push_back(cluster);
      
      //Fill TTree:
      fcheated_cluster_size[ii]=cluster->Hits().size();
      
      if(cluster->View()==geo::kU){
      fclusters_planeNo_true[ii]=0;
      fStart_pt_w_true[ii]=cluster->StartPos()[0];
      fStart_pt_t_true[ii]=cluster->StartPos()[1];
      //std::cout<<"plane=0, cheated cl# "<<ii<<" size= "<<cluster->Hits().size()<<" startPos_w= "<<cluster->StartPos()[0]<<" startPos_t= "<<cluster->StartPos()[1]<<std::endl;
      
      }
       else if(cluster->View()==geo::kV){
      fclusters_planeNo_true[ii]=1;
      fStart_pt_w_true[ii]=cluster->StartPos()[0];
      fStart_pt_t_true[ii]=cluster->StartPos()[1];
     // std::cout<<"plane=1, cheated cl# "<<ii<<" size= "<<cluster->Hits().size()<<" startPos_w= "<<cluster->StartPos()[0]<<" startPos_t= "<<cluster->StartPos()[1]<<std::endl;
      
      }
      
      
      
     
      if(cluster->Hits().size()>2 && cluster->View()==geo::kU){
      
      fcheatedCl_p0++;
      
     // std::cout<<cluster->StartPos()[0]<<", "<<cluster->StartPos()[1]<<" --> "<<cluster->EndPos()[0]<<", "<<cluster->EndPos()[1]<<" No hits= "<<cluster->Hits().size()<<std::endl;
      
      
      if(fabs(cluster->StartPos()[0]-fwire_vertex_true[0])<6 && fabs(cluster->StartPos()[1]-ftime_vertex[0])<90 ){
      fcheatedCl_near_vertex_p0++;
      
      }
      
      
      }
      if(cluster->Hits().size()>2 && cluster->View()==geo::kV){
      
      fcheatedCl_p1++;
      
       if(fabs(cluster->StartPos()[0]-fwire_vertex_true[1])<6 && fabs(cluster->StartPos()[1]-ftime_vertex[1])<90 ){
      fcheatedCl_near_vertex_p1++;
      
      }
      
      }
      
     //*********************************
     // find out what particle each cluster belongs to:
     
     
     hits=cluster->Hits();
     
    for(unsigned int h=0; h<hits.size(); h++){
    
    //std::cout<<"hits[h] channel= "<<hits[h]->Wire()->RawDigit()->Channel()<<std::endl;
   
    std::vector<cheat::TrackIDE> trackides = cheat::BackTracker::HitToTrackID(*scs[hits[h]->Channel()], hits[h]);
   
    std::vector<cheat::TrackIDE>::iterator idesitr = trackides.begin();
    //std::cout<<"trackides= "<<trackides.size();
    
    
   
    while( idesitr != trackides.end() ){
    
    vec_trackid.push_back((*idesitr).trackID);
    //std::cout<<"EngFraction= "<<(*idesitr).energyFrac<<std::endl;
    const sim::Particle* particle = _particleList.at( (*idesitr).trackID);
    int pdg = particle->PdgCode();
    //std::cout<<"pdg= "<<pdg<<std::endl;
    
    if(pdg==2212 || pdg==-2212){
       if(cluster->View()==geo::kU) proton_hit_ind++;
       else if(cluster->View()==geo::kV) proton_hit_coll++;
    
    }
    
    idesitr++;
    }//trackids
    
    vec_trackid.clear();
    
    }//hits
      
      if (cluster->View()==geo::kU && hits.size()>=2 && (proton_hit_ind/hits.size())>0.7){
      proton_track_ind++;
      std::cout<<" cluster #"<<ii<<" is a proton track! (IND)"<<std::endl;
      }
      else if (cluster->View()==geo::kV && hits.size()>=2 && (proton_hit_coll/hits.size())>0.7){
      proton_track_coll++;
      std::cout<<" cluster #"<<ii<<" is a proton track! (COLL)"<<std::endl;
      }
     proton_hit_ind=0;
     proton_hit_coll=0;
      
    }//clusters
    
    fno_clusters_true=cheatedclusterListHandle->size();
    
std::cout<<"Total No of CHEATED clusters= "<<cheatedclusterListHandle->size()<<std::endl;
std::cout<<"for plane 0:"<<fcheatedCl_p0<<std::endl;
std::cout<<"for plane 1:"<<fcheatedCl_p1<<std::endl;
std::cout<<"Total No of CHEATED clusters ***NEAR THE VERTEX*** :"<<std::endl;
std::cout<<"for plane 0:"<<fcheatedCl_near_vertex_p0<<std::endl;
std::cout<<"for plane 1:"<<fcheatedCl_near_vertex_p1<<std::endl;
std::cout<<"no of proton tracks for plane 0: "<<proton_track_ind<<std::endl;
std::cout<<"no of proton tracks for plane 1: "<<proton_track_coll<<std::endl;

//............ End of analyzing cheated clusters............................

 
 
 
 fNoProtonTracks_p0_cheatedCl->Fill(proton_track_ind);
 fNoProtonTracks_p1_cheatedCl->Fill(proton_track_coll);
 
fdiff_no_vertex_clusters_p0->Fill(fkingaCl_near_vertex_p0-fcheatedCl_near_vertex_p0);
 fdiff_no_vertex_clusters_p1->Fill(fkingaCl_near_vertex_p1-fcheatedCl_near_vertex_p1);
 fdiff_no_vertex_linemergerclusters_p0->Fill(flinemergerCl_near_vertex_p0-fcheatedCl_near_vertex_p0);
 fdiff_no_vertex_linemergerclusters_p1->Fill(flinemergerCl_near_vertex_p1-fcheatedCl_near_vertex_p1);
 
 fdiff_time_vtx_p0->Fill(ftime_vertex[0]-ftime_vertex_reco[0]);
      
 fdiff_wire_vtx_p0->Fill(fwire_vertex[0]-fwire_vertex_reco[0]);
      
 fdiff_wire_vtx_p1->Fill(fwire_vertex[1]-fwire_vertex_reco[1]);
 
 
 
 
 
 
 //--------------------------------------------------------------//
 //        NOW I WILL GET INFO FROM GEANT4 TO FIND OUT HOW MANY 
 //        PARTICLES WE CAN REALLY SEE IN OUR DETECTOR
 //        this is needed if you want to confirm that kingaclusters 
 ///       can correctly count tracks:
 //--------------------------------------------------------------//
 
 
 
 
 art::Handle< std::vector<sim::Particle> > geant_list;
    evt.getByLabel (fLArGeantModuleLabel,geant_list);
 
  art::PtrVector<sim::Particle> geant_part;
   for (unsigned int ii = 0; ii <  geant_list->size(); ++ii)
    {
      art::Ptr<sim::Particle> p(geant_list,ii);
      geant_part.push_back(p);
    } 
 std::string pri ("primary");
 int primary=0;
 //determine the number of primary particles from geant:
  
  for( unsigned int i = 0; i < geant_part.size(); ++i ){
   
    if(geant_part[i]->Process()==pri){
    primary++;
    }
   
   }
  
  fno_primaries=primary;
  
 std::cout<<"Geant4 list: "<<std::endl;
 
 for( unsigned int i = 0; i < geant_part.size(); ++i ){
   
   if(geant_part[i]->Process()==pri){
   
   std::cout<<"StatusCode= "<<geant_part[i]->StatusCode()<<" Mother= "<<geant_part[i]->Mother()<<std::endl;
   
    // fprimaries_pdg.push_back(geant_part[i]->PdgCode());
//     std::cout<<"geant_part[i]->E()= "<<geant_part[i]->E()<<std::endl;
//     fEng.push_back(geant_part[i]->E());
//     
   
    fprimaries_pdg[i]=geant_part[i]->PdgCode();
    
    fEng[i]=geant_part[i]->E();
    fPx[i]=geant_part[i]->Px();
   
    fPy[i]=geant_part[i]->Py();
    fPz[i]=geant_part[i]->Pz();
    
   fStartPointx[i]=geant_part[i]->Vx();
   fStartPointy[i]=geant_part[i]->Vy();
   fStartPointz[i]=geant_part[i]->Vz();
   fEndPointx[i]=geant_part[i]->EndPoint()[0];
   fEndPointy[i]=geant_part[i]->EndPoint()[1];
   fEndPointz[i]=geant_part[i]->EndPoint()[2];
   
   fNumberDaughters[i]=geant_part[i]->NumberDaughters();
   
   
  std::cout<<"pdg= "<<geant_part[i]->PdgCode()<<" trackId= "<<geant_part[i]->TrackId()<<" mother= "<<geant_part[i]->Mother()<<" NumberDaughters()= "<<geant_part[i]->NumberDaughters()<<" process= "<<geant_part[i]->Process()<<std::endl;
     
     }
     
     }
 
 
 
 
 
 
 
  
  fTree->Fill();
  
}

