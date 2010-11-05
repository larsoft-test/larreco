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
 
#include "Filters/ChannelFilter.h"
#include "SimulationBase/simbase.h"
#include "RawData/RawDigit.h"
#include "RecoBase/recobase.h"
#include "Geometry/geo.h"

cluster::HoughLineFinderAna::HoughLineFinderAna(edm::ParameterSet const& pset) : 
  fHoughModuleLabel    (pset.getParameter< std::string >("HoughModuleLabel")),
  fDigitModuleLabel    (pset.getParameter< std::string >("DigitModuleLabel")),
  fHitsModuleLabel      (pset.getParameter< std::string >("HitsModuleLabel")),
  fm_run(0), 
  fm_event(0), 
  fm_plane(0), 
  fm_clusterslope(0), 
  fm_clusterintercept(0),
  fm_clusterid(0), 
  fm_wirespan(0), 
  fm_sizeClusterZ(10000), 
  fm_sizeHitZ(10000)
{
}

//-------------------------------------------------
cluster::HoughLineFinderAna::~HoughLineFinderAna()
{

}

//-------------------------------------------------
void cluster::HoughLineFinderAna::beginJob(edm::EventSetup const&)
{


    // get access to the TFile service
     edm::Service<edm::TFileService> tfs;
     ftree= tfs->make<TTree>("HoughTree","HoughTree");
     fm_hitidZ = new Int_t[fm_sizeHitZ];
     fm_mipZ = new Float_t[fm_sizeHitZ];
     fm_drifttimeZ = new Float_t[fm_sizeHitZ];
     fm_widthZ = new Float_t[fm_sizeHitZ];
     fm_upadcZ = new Float_t[fm_sizeHitZ];
     fm_wireZ = new Int_t[fm_sizeHitZ];
     ftree->Branch("run", &fm_run, "run/I");
     ftree->Branch("run_timestamp", &fm_run_timestamp, "run_timestamp/l"); //l is for ULong64_t
     ftree->Branch("event", &fm_event, "event/I");
     ftree->Branch("plane",&fm_plane,"plane/I");
     ftree->Branch("clusterid",&fm_clusterid,"clusterid/I");
     ftree->Branch("clusterslope",&fm_clusterslope,"clusterslope/F");
     ftree->Branch("clusterintercept",&fm_clusterintercept,"clusterintecept/F");
     ftree->Branch("wirespan",&fm_wirespan,"wirespan/I");
     ftree->Branch("numberHits",&fm_sizeHitZ,"numberHits/I");
     ftree->Branch("numberClusters",&fm_sizeClusterZ,"numberClusters/I");
     ftree->Branch("hitidZ",fm_hitidZ,"hitidZ[numberHits]/I");
     ftree->Branch("wireZ",fm_wireZ,"wireZ[numberHits]/I");
     ftree->Branch("mipZ",fm_mipZ,"mipZ[numberHits]/F");
     ftree->Branch("drifttimeZ",fm_drifttimeZ,"drifttitmeZ[numberHits]/F");
     ftree->Branch("widthZ",fm_widthZ,"widthZ[numberHits]/F");
}


void cluster::HoughLineFinderAna::analyze(const edm::Event& evt, edm::EventSetup const&)
{

  edm::Handle< std::vector<recob::Cluster> > hlfListHandle;
  evt.getByLabel(fHoughModuleLabel,hlfListHandle);
  edm::Handle< std::vector<recob::Hit> > hitListHandle;
  evt.getByLabel(fHitsModuleLabel,hitListHandle);

  edm::PtrVector<recob::Cluster> clusters;   
//   edm::PtrVector<recob::Hit> hits;// unused, as yet. EC, 5-Oct-2010.
    
  for (unsigned int ii = 0; ii <  hlfListHandle->size(); ++ii)
    {
      edm::Ptr<recob::Cluster> cluster(hlfListHandle,ii);
      clusters.push_back(cluster);
    }
  
    std::cout << "run    : " << evt.id().run() << std::endl;
    //std::cout << "subrun : " << evt.subRun() << std::endl;
    std::cout << "event  : " << evt.id().event() << std::endl;
    fm_run=evt.id().run();
    fm_event=evt.id().event();
    fm_run_timestamp=evt.time().value(); // won't cast, EC, 7-Oct-2010.

    unsigned int firstwire=0;
    unsigned int lastwire=0;
    unsigned int p;
    fm_sizeClusterZ=0;
    fm_sizeHitZ=0;
    unsigned int wire=0;

    edm::Service<geo::Geometry> geo;

    for(unsigned int plane=0;plane < geo->Nplanes();++plane)
 {
      fm_plane=plane;
      fm_sizeClusterZ=clusters.size();
  for(unsigned int j=0; j<clusters.size();++j) 
  {
       fm_clusterid=clusters[j]->ID();
       edm::PtrVector<recob::Hit> _hits=clusters[j]->Hits(plane);
	   fm_clusterslope=(Float_t)clusters[j]->Slope();
       fm_clusterintercept=(Float_t)clusters[j]->Intercept();
   if(_hits.size()!=0)
   {
	  geo->ChannelToWire(_hits[0]->Wire()->RawDigit()->Channel(), p, firstwire);
	  geo->ChannelToWire(_hits[_hits.size()-1]->Wire()->RawDigit()->Channel(), p, lastwire);
	  fm_wirespan=lastwire-firstwire;
	  fm_sizeHitZ=_hits.size();
	  
	  for(unsigned int i = 0; i < _hits.size(); ++i) 
	  {	     
	  
	     geo->ChannelToWire(_hits[i]->Wire()->RawDigit()->Channel(), p, wire);
         fm_hitidZ[i]=i;         
         fm_wireZ[i]=wire;
         fm_mipZ[i]=(Float_t)_hits[i]->MIPs();
         fm_drifttimeZ[i]= (Float_t)_hits[i]->CrossingTime();
         fm_widthZ[i]=(Float_t)_hits[i]->EndTime()-_hits[i]->StartTime();
         fm_upadcZ[i]=(Float_t)_hits[i]->UpADC();
	  } 
    
   ftree->Fill();  
   }
  }
 }


}

