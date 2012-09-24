////////////////////////////////////////////////////////////////////////
//
// DBSCANfinder.cxx
//
// kinga.partyka@yale.edu
//
//  This algorithm finds clusters of hits, they can be of arbitrary shape.You need to specify 2(3) parameters: 
// epsilon, epsilon2 and MinPoints as explained in the corresponding xml file.In my comments a 'point' reference 
// appears quite often. A 'point' is basically a simple hit which only contains wire and time information. This 
// algorithm is based on DBSCAN(Density Based Spatial Clustering of Applications with Noise): M. Ester, H.-P. Kriegel, 
// J. Sander, and X. Xu, A density-based algorithm for discovering clusters in large spatial databases with noise, 
// Second International Conference on Knowledge Discovery and Data Mining, pp. 226-231, AAAI Press. 1996.
// ( Some of this code is from "Antonio Gulli's coding playground")  
////////////////////////////////////////////////////////////////////////


//Framework includes:
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "Utilities/AssociationUtil.h"
#include "ClusterFinder/DBcluster.h"
#include "Filters/ChannelFilter.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TGeoManager.h"
#include "TH1.h"

//-------------------------------------------------
cluster::DBcluster::DBcluster(fhicl::ParameterSet const& pset)
  : fDBScan(pset.get< fhicl::ParameterSet >("DBScanAlg")) 
{  
  this->reconfigure(pset);
  produces< std::vector<recob::Cluster> >();  
  produces< art::Assns<recob::Cluster, recob::Hit> >();
}

//-------------------------------------------------
cluster::DBcluster::~DBcluster()
{
}

//-------------------------------------------------
void cluster::DBcluster::reconfigure(fhicl::ParameterSet const& p)
{
  fhitsModuleLabel = p.get< std::string >("HitsModuleLabel");
  fDBScan.reconfigure(p.get< fhicl::ParameterSet >("DBScanAlg"));
}

//-------------------------------------------------
void cluster::DBcluster::beginJob(){
  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;

  fhitwidth= tfs->make<TH1F>(" fhitwidth","width of hits in cm", 50000,0 ,5  );
  fhitwidth_ind_test= tfs->make<TH1F>("fhitwidth_ind_test","width of hits in cm", 50000,0 ,5  );
  fhitwidth_coll_test= tfs->make<TH1F>("fhitwidth_coll_test","width of hits in cm", 50000,0 ,5  );
    
}

//-----------------------------------------------------------------
void cluster::DBcluster::produce(art::Event& evt)
{
   
  //get a collection of clusters   
  std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>);
  std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);

  art::ServiceHandle<geo::Geometry> geom;

  art::Handle< std::vector<recob::Hit> > hitcol;
  evt.getByLabel(fhitsModuleLabel,hitcol);
  
  // loop over all hits in the event and look for clusters (for each plane)
  art::PtrVector<recob::Hit> allhits;

  // get the ChannelFilter
  filter::ChannelFilter chanFilt;
      
  unsigned int p(0),w(0), t(0), c(0), channel(0);
  for(unsigned int cstat = 0; cstat < geom->Ncryostats(); ++cstat){
    for(unsigned int tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc){
      for(unsigned int plane = 0; plane < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++plane){
	geo::SigType_t sigType = geom->Cryostat(cstat).TPC(tpc).Plane(plane).SignalType();
	for(size_t i = 0; i< hitcol->size(); ++i){
      
	  art::Ptr<recob::Hit> hit(hitcol, i);
      
	  channel=hit->Wire()->RawDigit()->Channel();
	  geom->ChannelToWire(channel, c, t, p, w);
    
	  if(p == plane && t == tpc && c == cstat) allhits.push_back(hit);
	
	}  
      
	fDBScan.InitScan(allhits, chanFilt.SetOfBadChannels());

	//----------------------------------------------------------------
	for(unsigned int j = 0; j < fDBScan.fps.size(); ++j){
	  
	  if(allhits.size() != fDBScan.fps.size()) break;
	  
	  fhitwidth->Fill(fDBScan.fps[j][2]);
	  
	  if(sigType == geo::kInduction)  fhitwidth_ind_test->Fill(fDBScan.fps[j][2]);
	  if(sigType == geo::kCollection) fhitwidth_coll_test->Fill(fDBScan.fps[j][2]);
	}
 
	//*******************************************************************
	fDBScan.run_cluster();
	
	for(size_t i = 0; i < fDBScan.fclusters.size(); ++i){
	  art::PtrVector<recob::Hit> clusterHits;
	  double totalQ = 0.;

	  for(size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j){	  
	    if(fDBScan.fpointId_to_clusterId[j]==i){
	      clusterHits.push_back(allhits[j]);
	      totalQ += clusterHits.back()->Charge();
	    }
	  }
         
	  ////////
	  if (clusterHits.size()>0){
	    
	    /// \todo: need to define start and end positions for this cluster and slopes for dTdW, dQdW
	    unsigned int sw = 0;
	    unsigned int ew = 0;
	    geom->ChannelToWire(clusterHits[0]->Wire()->RawDigit()->Channel(), c, t, p, sw);
	    geom->ChannelToWire(clusterHits[clusterHits.size()-1]->Wire()->RawDigit()->Channel(), c, t, p, ew);
	 
	    recob::Cluster cluster(sw*1., 0.,
				   clusterHits[0]->PeakTime(), clusterHits[0]->SigmaPeakTime(),
				   ew*1., 0.,
				   clusterHits[clusterHits.size()-1]->PeakTime(), clusterHits[clusterHits.size()-1]->SigmaPeakTime(),
				   -999., 0., 
				   -999., 0.,
				   totalQ,
				   geom->Cryostat(c).TPC(t).Plane(p).View(),
				   ccol->size());
	    
	    ccol->push_back(cluster);

	    // associate the hits to this cluster
	    util::CreateAssn(*this, evt, *(ccol.get()), clusterHits, *(assn.get()));
	    
	    clusterHits.clear();
	    
	  }//end if clusterHits has at least one hit
   
	}//end loop over fclusters
	
	allhits.clear();
      } // end loop over planes
    } // end loop over tpcs
  } // end loop over cryostats

  mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
  mf::LogVerbatim("Summary") << "DBcluster Summary:";
  for(unsigned int i = 0; i<ccol->size(); ++i) mf::LogVerbatim("Summary") << ccol->at(i) ;

  evt.put(std::move(ccol));
  evt.put(std::move(assn));

  return;
}
