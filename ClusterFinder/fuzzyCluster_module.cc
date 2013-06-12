////////////////////////////////////////////////////////////////////////
// $Id: DBSCANfinderAna.cxx,v 1.36 2010/09/15  bpage Exp $
//
// \file fuzzyCluster_module.cc
//
// \author kinga.partyka@yale.edu
//
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "TGeoManager.h"
#include "TH1.h"

//Framework includes:
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
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

// LArSoft includes
#include "RawData/raw.h"
#include "RawData/RawDigit.h"
#include "RecoAlg/fuzzyClusterAlg.h"
#include "RecoAlg/CornerFinderAlg.h"
#include "Filters/ChannelFilter.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/SeedCreator.h"


class TH1F;

namespace cluster{
   
  //--------------------------------------------------------------- 
  class fuzzyCluster : public art::EDProducer
  {
  public:
    explicit fuzzyCluster(fhicl::ParameterSet const& pset); 
    ~fuzzyCluster();
    void produce(art::Event& evt);
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);
    
  private:
       
    TH1F *fhitwidth;
    TH1F *fhitwidth_ind_test;  
    TH1F *fhitwidth_coll_test;  
  
    std::string fhitsModuleLabel;
    std::string fCalDataModuleLabel;
   
    fuzzyClusterAlg ffuzzyCluster; ///< object that implements the fuzzy cluster algorithm
    CornerFinderAlg fcornerfinder; ///< object that implements the corner finder algorithm
  };

}

namespace cluster{


  //-------------------------------------------------
  fuzzyCluster::fuzzyCluster(fhicl::ParameterSet const& pset) :
    ffuzzyCluster(pset.get< fhicl::ParameterSet >("fuzzyClusterAlg")),
    fcornerfinder(pset.get< fhicl::ParameterSet >("CornerFinderAlg"))
  {  
    this->reconfigure(pset);
    produces< std::vector<recob::Cluster> >();  
    produces< art::Assns<recob::Cluster, recob::Hit> >();
  
    // Create random number engine needed for PPHT
    createEngine(SeedCreator::CreateRandomNumberSeed(),"HepJamesRandom");
  }
  
  //-------------------------------------------------
  fuzzyCluster::~fuzzyCluster()
  {
  }
  
  //-------------------------------------------------
  void fuzzyCluster::reconfigure(fhicl::ParameterSet const& p)
  {
    fhitsModuleLabel  = p.get< std::string >("HitsModuleLabel");
    fCalDataModuleLabel  = p.get< std::string >("CalDataModuleLabel");
    ffuzzyCluster.reconfigure(p.get< fhicl::ParameterSet >("fuzzyClusterAlg"));
  }
  
  //-------------------------------------------------
  void fuzzyCluster::beginJob(){
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;
  
    fhitwidth= tfs->make<TH1F>(" fhitwidth","width of hits in cm", 50000,0 ,5  );
    fhitwidth_ind_test= tfs->make<TH1F>("fhitwidth_ind_test","width of hits in cm", 50000,0 ,5  );
    fhitwidth_coll_test= tfs->make<TH1F>("fhitwidth_coll_test","width of hits in cm", 50000,0 ,5  );
      
  }
  
  //-----------------------------------------------------------------
  void fuzzyCluster::produce(art::Event& evt)
  {
     
    //get a collection of clusters   
    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>);
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);
  
    art::ServiceHandle<geo::Geometry> geom;
  
    art::Handle< std::vector<recob::Hit> > hitcol;
    evt.getByLabel(fhitsModuleLabel,hitcol);

    art::Handle< std::vector<recob::Wire> > wirecol;
    evt.getByLabel(fCalDataModuleLabel,wirecol);
   
    // loop over all hits in the event and look for clusters (for each plane)
    std::vector<art::Ptr<recob::Hit> > allhits;

    // loop over all end points in the event to help look for clusters (for each plane)
    //std::vector<art::Ptr<recob::EndPoint2D> > allends;
    std::vector<recob::EndPoint2D> allends;

    // Set event number as the random number seed needed for PPHT
    //std::cout << "Event number check: " << evt.event() << std::endl;
    //art::ServiceHandle<art::RandomNumberGenerator> rng;
    //CLHEP::HepRandomEngine &engine = rng->getEngine();
    //engine.setSeed((long int)evt.event(),0);
  
    std::vector<recob::EndPoint2D> endcol; 

    // Pass information into CornerFinder
    //fcornerfinder.TakeInRaw(evt); //this makes the histograms
    //fcornerfinder.get_feature_points_LineIntegralScore(endcol); //this fills our vector with feature points

    // get the ChannelFilter
    filter::ChannelFilter chanFilt;
        
    for(unsigned int cstat = 0; cstat < geom->Ncryostats(); ++cstat){
      for(unsigned int tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc){
        for(unsigned int plane = 0; plane < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++plane){
	  geo::SigType_t sigType = geom->Cryostat(cstat).TPC(tpc).Plane(plane).SignalType();
	  for(size_t i = 0; i< hitcol->size(); ++i){
	    art::Ptr<recob::Hit> hit(hitcol, i);
	    if(hit->WireID().Plane    == plane && 
	       hit->WireID().TPC      == tpc   && 
	       hit->WireID().Cryostat == cstat) allhits.push_back(hit);  
	  }  


	  //for(size_t i = 0; i< endcol.size(); ++i){
	    //if(endcol[i].WireID().Plane    == plane && 
	       //endcol[i].WireID().TPC      == tpc   && 
	       //endcol[i].WireID().Cryostat == cstat){allends.push_back(endcol[i]);  
	      //std::cout << "wire: " << endcol[i].WireID().Wire << " time: " << endcol[i].DriftTime() << std::endl;
	    //}
	  //}  



          //Begin clustering with fuzzy
          
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
	      if(ffuzzyCluster.fpointId_to_clusterId[j]==i){ 
		clusterHits.push_back(allhits[j]);
		totalQ += clusterHits.back()->Charge();
	      }
	    } 
  
  
	    ////////
	    if (clusterHits.size()>0){
	      /// \todo: need to define start and end positions for this cluster and slopes for dTdW, dQdW
	      unsigned int sw = clusterHits[0]->WireID().Wire;
	      unsigned int ew = clusterHits[clusterHits.size()-1]->WireID().Wire;
  	 
	      recob::Cluster cluster(sw*1., 0.,
				     clusterHits[0]->PeakTime(), clusterHits[0]->SigmaPeakTime(),
				     ew*1., 0.,
				     clusterHits[clusterHits.size()-1]->PeakTime(), clusterHits[clusterHits.size()-1]->SigmaPeakTime(),
				     -999., 0., 
				     -999., 0.,
				     totalQ,
				     clusterHits[0]->View(),
				     ccol->size());
  	    
	      ccol->push_back(cluster);
  
	      // associate the hits to this cluster
	      util::CreateAssn(*this, evt, *ccol, clusterHits, *assn);
  	    
	      clusterHits.clear();
  	    
	    }//end if clusterHits has at least one hit
     
	  }//end loop over fclusters
  	
	  allhits.clear();
        } // end loop over planes
      } // end loop over tpcs
    } // end loop over cryostats
  
    mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    mf::LogVerbatim("Summary") << "fuzzyCluster Summary:";
    for(size_t i = 0; i<ccol->size(); ++i) mf::LogVerbatim("Summary") << ccol->at(i) ;
  
    evt.put(std::move(ccol));
    evt.put(std::move(assn));
  
    return;
  } // end produce
  
} // end namespace

namespace cluster{

  DEFINE_ART_MODULE(fuzzyCluster);
  
} 







