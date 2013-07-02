////////////////////////////////////////////////////////////////////////
// Class:       ClusterCrawler
// Module Type: producer
// File:        ClusterCrawler_module.cc
//
// Generated at Fri Jun  7 09:44:09 2013 by Bruce Baller using artmod 
// from cetpkgsupport v1_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include <memory>

//LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "Utilities/AssociationUtil.h"
#include "Filters/ChannelFilter.h"
#include "RecoAlg/ClusterCrawlerAlg.h"


namespace cluster {
  class ClusterCrawler;
}

class cluster::ClusterCrawler : public art::EDProducer {

  public:
    explicit ClusterCrawler(fhicl::ParameterSet const & pset);
    virtual ~ClusterCrawler();

    void reconfigure(fhicl::ParameterSet const & pset) override;
    void produce(art::Event & evt) override;
    void beginJob();

  private:
    std::string fhitsModuleLabel;
    ClusterCrawlerAlg fCCAlg; // define ClusterCrawlerAlg object
    
};


namespace cluster {

  ClusterCrawler::ClusterCrawler(fhicl::ParameterSet const& pset)
    : fCCAlg(pset.get< fhicl::ParameterSet >("ClusterCrawlerAlg")) 
  {  
    this->reconfigure(pset);
    produces< std::vector<recob::Cluster> >();  
    produces< art::Assns<recob::Cluster, recob::Hit> >();
  }

  ClusterCrawler::~ClusterCrawler()
  {
  }

  void ClusterCrawler::reconfigure(fhicl::ParameterSet const & pset)
  {
    fhitsModuleLabel = pset.get< std::string >("HitsModuleLabel");
    fCCAlg.reconfigure(pset.get< fhicl::ParameterSet >("ClusterCrawlerAlg"));
  }
  
  void ClusterCrawler::beginJob(){
  }
  
  struct SortByWire {
    bool operator() (recob::Hit const& h1, recob::Hit const& h2) const 
    { 
      return h1.Wire()->Channel() < h2.Wire()->Channel();
    }
  };
  
  void ClusterCrawler::produce(art::Event & evt)
  {
    art::ServiceHandle<geo::Geometry> geo;
  
    art::Handle< std::vector<recob::Hit> > hitcol;
    evt.getByLabel(fhitsModuleLabel,hitcol);

    std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);
    
    // loop over all hits in the event and look for clusters in each plane
    art::PtrVector<recob::Hit> plnhits;
  
    // get the ChannelFilter
    filter::ChannelFilter chanFilt;
        
    for(unsigned int cstat = 0; cstat < geo->Ncryostats(); ++cstat){
      for(unsigned int tpc = 0; tpc < geo->Cryostat(cstat).NTPC(); ++tpc){
        for(unsigned int plane = 0; plane < geo->Cryostat(cstat).TPC(tpc).Nplanes(); ++plane){
          // load the hits in this plane
//          geo::SigType_t sigType = geo->Cryostat(cstat).TPC(tpc).Plane(plane).SignalType();
// debugging
//          if(sigType != geo::kCollection) continue;
          for(size_t i = 0; i< hitcol->size(); ++i){
            art::Ptr<recob::Hit> hit(hitcol, i);
            if(hit->WireID().Plane    == plane && 
               hit->WireID().TPC      == tpc   && 
               hit->WireID().Cryostat == cstat) plnhits.push_back(hit);
          }  // i
          plnhits.sort(cluster::SortByWire());
          // vector of temporary clusters returned from RunCrawler
          std::vector< ClusterStore > tcl; 
          // make clusters and fill the tcl vector
          fCCAlg.RunCrawler(plnhits, plane, tcl);
          for(std::vector<ClusterStore>::const_iterator it = tcl.begin();
            it != tcl.end(); ++it) {
            ClusterStore clstr = *it;
            // ignore deleted clusters
            if(clstr.ID < 0) continue;
            std::vector<int>::const_iterator iht = clstr.tclhits.begin();
            int hits = *iht;
            int startwire = plnhits[hits]->WireID().Wire;
            double starttime = plnhits[hits]->PeakTime();
            iht = clstr.tclhits.end()-1;
            int hite = *iht;
            int endwire = plnhits[hite]->WireID().Wire;
            double endtime = plnhits[hite]->PeakTime();
            art::PtrVector<recob::Hit> clusterHits;
            double totalQ = 0.;
            // find the average charge at the low and high wire ends
            // this will be used to determine the start and end of the cluster
            double LoWirChg = 0.;
            double HiWirChg = 0.;
            // determine how many hits to use to sum the charge
            int ncnt = clstr.tclhits.size() / 2;
            if(ncnt > 3) ncnt = 3;
            for(std::vector<int>::const_iterator itt = clstr.tclhits.begin();
                itt != clstr.tclhits.end(); ++itt) {
              int hit = *itt;
              int wire = plnhits[hit]->WireID().Wire;
              // Note that the endwire < startwire since tracking is from DS to US
              if(wire < endwire + ncnt) {
                LoWirChg += plnhits[hit]->Charge();
              }
              if(wire > startwire - ncnt - 1) {
                HiWirChg += plnhits[hit]->Charge();
              }
              totalQ += plnhits[hit]->Charge();
              clusterHits.push_back(plnhits[hit]);
            } // hit iterator
            // assume the track is going in the direction of increasing wire number
            bool GoingPos = true;
            // unless the charge is much higher at the low wire end
                  if(LoWirChg > 2 * HiWirChg) GoingPos = false;
            // swap the start and end
            double slope;
            if(GoingPos) {
              int itmp = startwire;
              startwire = endwire;
              endwire = itmp;
              double tmp = starttime;
              starttime = endtime;
              endtime = tmp;
              slope = clstr.slpstart;
            } else {
              slope = clstr.slpend;
            }
            std::vector<recob::Cluster> newclus;
            std::vector< art::PtrVector<recob::Hit> > clusterhits;
            art::ServiceHandle<geo::Geometry> geom;

            recob::Cluster cluster(startwire, 0.,
                                  starttime, 0.,
                                  endwire, 0.,
                                  endtime, 0.,
                                  slope, 0.,
                                  -999.,0.,
                                  totalQ,
                                  geom->View((*plnhits.begin())->Channel()),
                                  clstr.ID);
            ccol->push_back(cluster);
    
            // associate the hits to this cluster
            util::CreateAssn(*this, evt, *(ccol.get()), clusterHits, *(assn.get()));
            clusterHits.clear();
          } // cluster iterator
          plnhits.clear();
        } // plane
      } // tpc
    } // cstat
      
    evt.put(std::move(ccol));
    evt.put(std::move(assn));
    
    return;

  } // produce
} // namespace


namespace cluster{

  DEFINE_ART_MODULE(ClusterCrawler);
  
} 

