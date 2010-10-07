////////////////////////////////////////////////////////////////////////
//
// HoughLineFinder class
//
// joshua.spitz@yale.edu
//
//  This algorithm is designed to find lines (Houghclusters) from clusters found by DBSCAN after deconvolution and hit finding.
//  The algorithm is based on: 
//  Queisser, A. "Computing the Hough Transform", C/C++ Users Journal 21, 12 (Dec. 2003).
//  Niblack, W. and Petkovic, D. On Improving the Accuracy of the Hough Transform", Machine Vision and Applications 3, 87 (1990)  
////////////////////////////////////////////////////////////////////////

#include "HoughLineFinderAna.h"
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
// ROOT includes
#include <TH1D.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TTree.h>
#include "TDatabasePDG.h"
#include "TSystem.h"

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>


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
#include "Filters/ChannelFilter.h"
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




cluster::HoughLineFinderAna::HoughLineFinderAna(edm::ParameterSet const& pset) : 
  fHoughModuleLabel    (pset.getParameter< string >("HoughModuleLabel")),
  fDigitModuleLabel    (pset.getParameter< string >("DigitModuleLabel")),
  fHitModuleLabel    (pset.getParameter< string >("HitModuleLabel")),
  m_run(0), m_event(0), m_plane(0), m_clusterid(0), m_wirespan(0), 
  m_sizeClusterZ(10000), m_sizeHitZ(10000)
{


}

//-------------------------------------------------
cluster::HoughLineFinderAna::~HoughLineFinderAna()
{
  delete [] m_hitidZ;
  delete [] m_mipZ;
  delete [] m_drifttimeZ;
  delete [] m_widthZ;
  delete [] m_upadcZ;
}

//-------------------------------------------------
void cluster::HoughLineFinderAna::beginJob(edm::EventSetup const&)
{

    // get access to the TFile service
     edm::Service<edm::TFileService> tfs;

     tree= tfs->make<TTree>("HoughTree","HoughTree");
     m_hitidZ = new Int_t[m_sizeHitZ];
     m_mipZ = new Float_t[m_sizeHitZ];
     m_drifttimeZ = new Float_t[m_sizeHitZ];
     m_widthZ = new Float_t[m_sizeHitZ];
     m_upadcZ = new Float_t[m_sizeHitZ];
     tree->Branch("run", &m_run, "run/I");
     tree->Branch("run_timestamp", &m_run_timestamp, "run_timestamp/D");
     tree->Branch("event", &m_event, "event/I");
     tree->Branch("plane",&m_plane,"plane/I");
     tree->Branch("clusterid",&m_clusterid,"clusterid/I");
     tree->Branch("wirespan",&m_wirespan,"wirespan/I");
     tree->Branch("numberHits",&m_sizeHitZ,"numberHits/I");
     tree->Branch("numberClusters",&m_sizeClusterZ,"numberClusters/I");
     tree->Branch("hitidZ",m_hitidZ,"hitidZ[numberHits]/I");
     tree->Branch("mipZ",m_mipZ,"mipZ[numberHits]/F");
     tree->Branch("drifttimeZ",m_drifttimeZ,"drifttitmeZ[numberHits]/F");
     tree->Branch("widthZ",m_widthZ,"widthZ[numberHits]/F");
}



cluster::HoughLineFinderAna::analyze(const edm::EventHandle& evt, edm::EventSetup const&)
{


  edm::Handle< std::vector<recob::Cluster> > hlfListHandle;
  evt.getByLabel(fHoughModuleLabel,hlfListHandle);
  edm::Handle< std::vector<raw::RawDigit> > rdListHandle;
  evt.getByLabel(fDigitModuleLabel,rdListHandle);
  edm::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitModuleLabel,hitListHandle);

  edm::PtrVector<recob::Cluster> clusters;   
  edm::PtrVector<raw::RawDigit> rawdigit;// unused, as yet. EC, 5-Oct-2010.
  edm::PtrVector<recob::Hit> hits;// unused, as yet. EC, 5-Oct-2010.
    
  for (unsigned int ii = 0; ii <  hlfListHandle->size(); ++ii)
    {
      edm::Ptr<const recob::Cluster> cluster(hlfListHandle,ii);
      clusters.push_back(cluster);
    }

  
    std::cout << "run    : " << evt.id().run() << std::endl;
    std::cout << "subrun : " << evt.id().subrun() << std::endl;
    std::cout << "event  : " << evt.id().event() << std::endl;
    m_run=evt.id().run();
    m_event=evt.id().event();
    m_run_timestamp=evt.time(); // no id(), EC, 5-Oct-2010.

    Int_t firstwire=0;
    Int_t lastwire=0;
    Int_t p;
    m_sizeClusterZ=0;
    m_sizeHitZ=0;

    edm::Service<geo::Geometry> geo;

    for(unsigned int plane=0;plane < geo->Nplanes();++plane)
 {
      m_plane=plane;
      m_sizeClusterZ=clusters.size();
  for(unsigned int j=0; j<clusters.size();++j) 
  {
       m_clusterid=clusters[j]->ID();
       edm::PtrVector<recob::Hit> _hits=clusters[j]->Hits(plane,-1);
       //std::vector<const recob::Hit*> _hits=clusters[j]->Hits(plane,-1);
	  
   if(_hits.size()!=0)
   {
	  geo->ChannelToWire(_hits[0]->Wire()->RawDigit()->Channel(), p, firstwire);
	  geo->ChannelToWire(_hits[_hits.size()-1]->Wire()->RawDigit()->Channel(), p, lastwire);
	  m_wirespan=lastwire-firstwire;
	  m_sizeHitZ=_hits.size();
	  for(unsigned int i = 0; i < _hits.size(); ++i) 
	  {
         m_hitidZ[i]=i;
         m_mipZ[i]=_hits[i]->MIPs();
         m_drifttimeZ[i]= _hits[i]->CrossingTime();
         m_widthZ[i]=_hits[i]->EndTime()-_hits[i]->StartTime();
         m_upadcZ[i]=_hits[i]->UpADC();
	  } 
    
   tree->Fill();  
   }
  }
 }


}

