////////////////////////////////////////////////////////////////////////
//
// ClusterFinder base class
//
// kinga.partyka@yale.edu
//
// joshua.spitz@yale.edu
////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional>
#include "Geometry/geo.h"
#include "ClusterFinder.h"

#include "TDatabasePDG.h"
#include "TSystem.h"

#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/ParameterSet/interface/ParameterSet.h" 
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h" 
#include "DataFormats/Common/interface/Handle.h" 
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



//-------------------------------------------------
cluster::ClusterFinder::ClusterFinder(edm::ParameterSet const& pset) : 
  fMakeClusters(1),
  fHitModuleLabel    (pset.getParameter< std::string >("HitModuleLabel")),
  fClusterFinderModuleLabel    (pset.getParameter< std::string >("ClusterFinderModuleLabel"))
{

  produces< std::vector<recob::Cluster> > ();

}

//-------------------------------------------------
cluster::ClusterFinder::~ClusterFinder()
{

}

void cluster::ClusterFinder::ClusterFinder::beginJob(edm::EventSetup const&)
{
     // get access to the TFile service
     edm::Service<edm::TFileService> tfs;


}

//-------------------------------------------------
void cluster::ClusterFinder::produce(edm::Event& evt, edm::EventSetup const&)
{
  /* 
     Breaking new ground. Let's stick this in ClusterFinder.h as a class member.
     Will it work? Dunno. But otherwise I don't see how this gets used as a
     base class, as it is now by DBScanfinder.cxx
  */

  std::auto_ptr<std::vector<recob::Cluster> > fClusterVecCol(new std::vector<recob::Cluster>);
  GetData(evt);


  evt.put(fClusterVecCol);



  
}

//-------------------------------------------------
void cluster::ClusterFinder::GetData(edm::Event& evt)
{


  // get the data from the appropriate place
  if(fMakeClusters==1){
    // create wires, and copy over raw data
    Make(evt);
  }
  else{
    // start with existing cal hits
    Read(evt);
  }
  
}

//-------------------------------------------------
void cluster::ClusterFinder::Make(edm::Event& evt) 
{
  edm::Service<geo::Geometry> geom;  

  // get hits
  edm::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitModuleLabel,hitListHandle);


  edm::PtrVector<recob::Hit> hitlist;
  for (unsigned int ii = 0; ii <  hitListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Hit> hitHolder(hitListHandle,ii);
      hitlist.push_back(hitHolder);
    }

  ///loop over all hits in the event and look for clusters
  edm::PtrVector<recob::Hit> hitlistInd, hitlistCol;
  unsigned int p(0),w(0), c(0);
  for(edm::PtrVectorItr<recob::Hit> a = hitlist.begin(); a!= hitlist.end(); a++){
    c=(*a)->Wire()->RawDigit()->Channel();
    geom->ChannelToWire(c,p,w);
    if(p == 0) hitlistInd.push_back(*a);
    else hitlistCol.push_back(*a);
  }  

  /*
    Not sure the point of this. EC, 7-Oct-2010.

  if(hitlistInd.size()>1){FindCluster(hitlistInd);}
  if(hitlistCol.size()>1){FindCluster(hitlistCol);}
  */

}

//-------------------------------------------------
void cluster::ClusterFinder::Read(edm::Event& evt)
{
  // get existing cluster collection specified via "Input" parameter.
    
  /* This confuses me. What is the module that this Base class can
     be counted upon to want to read from? EC, 5-Oct-2010.

     try{
     evt.Reco().Get(fInputFolder.c_str(),fClusterInput);
     } 
     catch(edm::Exception e){
     std::cerr << "Error retrieving Cluster list in "<<fInputFolder.c_str()
     <<std::endl;
     
     }


  */

  edm::Handle< std::vector<recob::Cluster> > clusterListHandle;
  evt.getByLabel(fClusterFinderModuleLabel,clusterListHandle);
  
  edm::PtrVector<recob::Cluster> fClustersInput;
  for (unsigned int ii = 0; ii <  clusterListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Cluster> clusterHolder(clusterListHandle,ii);
      fClustersInput.push_back(clusterHolder);
    }

}


