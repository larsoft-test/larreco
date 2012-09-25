////////////////////////////////////////////////////////////////////////
//
// fuzzyCluster.cxx
//
// Ben Carls, bcarls@fnal.gov
//
// This code looks for clusters using a fuzzy c-means algorithm. The
// clusters are then examined by the HoughClusAlg to identify Hough lines
// which can then be split off into their own clusters. See the webpage below
// for more information on the fuzzy clustering algorithm.
//
// http://home.dei.polimi.it/matteucc/Clustering/tutorial_html/cmeans.html
//
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
#include "CLHEP/Random/JamesRandom.h"

#include "Filters/ChannelFilter.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/SeedCreator.h"
#include "ClusterFinder/fuzzyCluster.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TGeoManager.h"
#include "TH1.h"

//-------------------------------------------------
//cluster::fuzzyCluster::fuzzyCluster(fhicl::ParameterSet const& pset) :
   //ffuzzyCluster(pset.get< fhicl::ParameterSet >("fuzzyClusterAlg")), 
   //fHCAlg(pset.get< fhicl::ParameterSet >("HoughClusAlg")) 
cluster::fuzzyCluster::fuzzyCluster(fhicl::ParameterSet const& pset) :
   ffuzzyCluster(pset.get< fhicl::ParameterSet >("fuzzyClusterAlg")) 
{  
  this->reconfigure(pset);
  produces< std::vector<recob::Cluster> >();  
  produces< art::Assns<recob::Cluster, recob::Hit> >();

  // Create random number engine needed for PPHT
  createEngine(SeedCreator::CreateRandomNumberSeed(),"HepJamesRandom");
}

//-------------------------------------------------
cluster::fuzzyCluster::~fuzzyCluster()
{
}

//-------------------------------------------------
void cluster::fuzzyCluster::reconfigure(fhicl::ParameterSet const& p)
{
  fhitsModuleLabel  = p.get< std::string >("HitsModuleLabel");
  ffuzzyCluster.reconfigure(p.get< fhicl::ParameterSet >("fuzzyClusterAlg"));
}

//-------------------------------------------------
void cluster::fuzzyCluster::beginJob(){
  // get access to the TFile service
  art::ServiceHandle<art::TFileService> tfs;

  fhitwidth= tfs->make<TH1F>(" fhitwidth","width of hits in cm", 50000,0 ,5  );
  fhitwidth_ind_test= tfs->make<TH1F>("fhitwidth_ind_test","width of hits in cm", 50000,0 ,5  );
  fhitwidth_coll_test= tfs->make<TH1F>("fhitwidth_coll_test","width of hits in cm", 50000,0 ,5  );
    
}

//-----------------------------------------------------------------
void cluster::fuzzyCluster::produce(art::Event& evt)
{
   
  //get a collection of clusters   
  std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>);
  std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);

  art::ServiceHandle<geo::Geometry> geom;

  art::Handle< std::vector<recob::Hit> > hitcol;
  evt.getByLabel(fhitsModuleLabel,hitcol);
  
  // loop over all hits in the event and look for clusters (for each plane)
  std::vector<art::Ptr<recob::Hit> > allhits;

  // Set event number as the random number seed needed for PPHT
  //std::cout << "Event number check: " << evt.event() << std::endl;
  //art::ServiceHandle<art::RandomNumberGenerator> rng;
  //CLHEP::HepRandomEngine &engine = rng->getEngine();
  //engine.setSeed((long int)evt.event(),0);

  //srand((unsigned)evt.event()); 
  
  
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
      
        //Begin clustering with fuzzy
        
        //Attempt to get number of clusters
        //std::cout << "Number of clusters: " << fHCAlg.Transform(allhits) << std::endl; 

	ffuzzyCluster.InitFuzzy(allhits, chanFilt.SetOfBadChannels());

	//----------------------------------------------------------------
	for(unsigned int j = 0; j < ffuzzyCluster.fps.size(); ++j){
	  
	  if(allhits.size() != ffuzzyCluster.fps.size()) break;
	  
	  fhitwidth->Fill(ffuzzyCluster.fps[j][2]);
	  
	  if(sigType == geo::kInduction)  fhitwidth_ind_test->Fill(ffuzzyCluster.fps[j][2]);
	  if(sigType == geo::kCollection) fhitwidth_coll_test->Fill(ffuzzyCluster.fps[j][2]);
	}
 
	//*******************************************************************
	ffuzzyCluster.run_fuzzy_cluster(allhits);

        //End clustering with fuzzy


	for(size_t i = 0; i < ffuzzyCluster.fclusters.size(); ++i){
          std::vector<art::Ptr<recob::Hit> > clusterHits;
	  double totalQ = 0.;
	  
	  for(size_t j = 0; j < ffuzzyCluster.fpointId_to_clusterId.size(); ++j){
            //std::cout << "ffuzzyCluster.fpointId_to_clusterId[j]: " << ffuzzyCluster.fpointId_to_clusterId[j] << " i: " << i << std::endl;
	    if(ffuzzyCluster.fpointId_to_clusterId[j]==i){ 
	      clusterHits.push_back(allhits[j]);
	      totalQ += clusterHits.back()->Charge();
            }
          } 


	  ////////
	  if (clusterHits.size()>0){
	    //std::cout << "i: " << i << std::endl; 
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
  mf::LogVerbatim("Summary") << "fuzzyCluster Summary:";
  for(unsigned int i = 0; i<ccol->size(); ++i) mf::LogVerbatim("Summary") << ccol->at(i) ;

  evt.put(std::move(ccol));
  evt.put(std::move(assn));

  return;
}
