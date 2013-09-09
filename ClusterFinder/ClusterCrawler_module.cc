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
#include "RecoBase/EndPoint2D.h"
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
    produces< std::vector<recob::EndPoint2D> >();
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
    
    std::vector<recob::EndPoint2D> vtxOut;
    std::unique_ptr<std::vector<recob::EndPoint2D> > vcol(new std::vector<recob::EndPoint2D>(vtxOut));
    
    // loop over all hits in the event and look for clusters in each plane
    art::PtrVector<recob::Hit> plnhits;

    for(unsigned int cstat = 0; cstat < geo->Ncryostats(); ++cstat){
      for(unsigned int tpc = 0; tpc < geo->Cryostat(cstat).NTPC(); ++tpc){
        for(unsigned int plane = 0; plane < geo->Cryostat(cstat).TPC(tpc).Nplanes(); ++plane){
          // load the hits in this plane
          plnhits.clear();
          for(size_t i = 0; i< hitcol->size(); ++i){
            art::Ptr<recob::Hit> hit(hitcol, i);
            if(hit->WireID().Plane    == plane && 
               hit->WireID().TPC      == tpc   && 
               hit->WireID().Cryostat == cstat) plnhits.push_back(hit);
          }  // i
          if(plnhits.size() < 4 || plnhits.size() > 32767) {
            plnhits.clear();
            continue;
          }
          plnhits.sort(cluster::SortByWire());
          // convert to const
          fCCAlg.RunCrawler(plnhits, plane);
          for(size_t it = 0; it < fCCAlg.tcl.size(); it++) {
            ClusterCrawlerAlg::ClusterStore clstr = fCCAlg.tcl[it];
            // ignore deleted clusters
            if(clstr.ID < 0) continue;
            art::PtrVector<recob::Hit> clusterHits;
            double totalQ = 0.;
            for(auto const &itt : clstr.tclhits) {
              if(itt < 0 || itt > plnhits.size() - 1) {
                std::cout<<"Bad itt "<<itt<<std::endl;
                continue;
              }
              totalQ += plnhits[itt]->Charge();
              clusterHits.push_back(plnhits[itt]);
            } // hit iterator
            recob::Cluster cluster((double)clstr.BeginWir, 0.,
                                   (double)clstr.BeginTim, 0.,
                                   (double)clstr.EndWir, 0.,
                                   (double)clstr.EndTim, 0.,
                                   (double)clstr.EndSlp, (double)clstr.EndSlpErr,
                                    -999.,0.,
                                   totalQ,
                                   plnhits[0]->View(),
                                   (int)clstr.ID);
            ccol->push_back(cluster);
            // associate the hits to this cluster
            util::CreateAssn(*this, evt, *ccol, clusterHits, *assn);
            clusterHits.clear();
          } // cluster iterator
          fCCAlg.tcl.clear();
          for(unsigned int iv = 0; iv < fCCAlg.vtx.size(); iv++) {
            ClusterCrawlerAlg::VtxStore Vtx = fCCAlg.vtx[iv];
            double drtime = Vtx.Time;
            unsigned int wire = Vtx.Wire;
            uint32_t chan = geo->PlaneWireToChannel(plane, wire, tpc, cstat);
            std::vector<geo::WireID> wIDvec = geo->ChannelToWire(chan);
            geo::WireID wID = wIDvec[0];
            double strenth = Vtx.Wght;
            recob::EndPoint2D myvtx(drtime, wID, strenth, (int)iv, plnhits[0]->View(), 0.);
            vcol->push_back(myvtx);
          } // iv
          fCCAlg.vtx.clear();
        } // plane
      } // tpc
    } // cstat
    
    evt.put(std::move(ccol));
    evt.put(std::move(assn));
    evt.put(std::move(vcol));

  } // produce
} // namespace


namespace cluster{

  DEFINE_ART_MODULE(ClusterCrawler)
  
} 

