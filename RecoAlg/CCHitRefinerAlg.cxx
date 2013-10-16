//////////////////////////////////////////////////////////////////////
///
/// CCHitRefiner class
///
/// Bruce Baller, baller@fnal.gov
///
/// Refines hits found by CCHitFinder, using information from
/// ClusterCrawlerAlg
///
////////////////////////////////////////////////////////////////////////


extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include <stdint.h>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h"   
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/EDProducer.h" 


// LArSoft Includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Hit.h"

// ROOT Includes 
#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "TVirtualFitter.h"

#include "RecoAlg/CCHitRefinerAlg.h"

namespace cluster{

//------------------------------------------------------------------------------
  CCHitRefinerAlg::CCHitRefinerAlg(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  void CCHitRefinerAlg::reconfigure(fhicl::ParameterSet const& pset)
  {
    fRefineHits         = pset.get<             bool   >("RefineHits");
    fBEChgRat           = pset.get<             float  >("BEChgRat");
  }

//------------------------------------------------------------------------------
  CCHitRefinerAlg::~CCHitRefinerAlg()
  {
  }
  
  void CCHitRefinerAlg::RunCCHitRefiner(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      CCHitFinderAlg::HitCuts& hitcuts,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl,
      std::vector<ClusterCrawlerAlg::VtxStore>& vtx, 
      ClusterCrawlerAlg& fCCAlg) 
  {
    // try to refine hits near vertices. Hits on clusters are assumed to be
    // in reverse wire order, ala ClusterCrawler, i.e. Begin = DS = large wire
    // number and End = US = small wire number.
    // This alg also defines the Begin and End of clusters. This code
    // may be called independently of the request to refine hits

    // swap Begin/End w/o refining hits?
    if(!fRefineHits && fBEChgRat > 0.) {
      cl2SetBeginEnd(allhits, tcl);
      return;
    }
  
    std::cout<<"CCHitRefiner "<<std::endl;
    
    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if(vtx[iv].Wire != 71) continue;
      plane = vtx[iv].CTP - vtx[iv].CTP / 10;
      // get the hit range --> should only do this once for each plane...
      fCCAlg.GetHitRange(allhits, vtx[iv].CTP, WireHitRange, fFirstWire, fLastWire);
//  std::cout<<"first wire "<<fFirstWire<<" "<<fLastWire<<std::endl;
      // list of clusters that are US (DS) of the vertex
      std::vector<unsigned short> clBeg;
      std::vector<unsigned short> clEnd;
      for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
        if(tcl[icl].ID < 0) continue;
        // clusters that end DS of the vtx
        if(tcl[icl].BeginVtx == iv) clBeg.push_back(icl);
        // clusters that end US of the vtx
        if(tcl[icl].EndVtx == iv) clEnd.push_back(icl);
      }
      std::cout<<"Vtx "<<iv<<" P:W:T "<<plane<<":"<<vtx[iv].Wire<<":"<<(int)vtx[iv].Time;
      std::cout<<" Beg clusters ";
      for(unsigned int ii = 0; ii <clBeg.size(); ++ii) {
        std::cout<<clBeg[ii]<<" ";
      }
      std::cout<<" End clusters ";
      for(unsigned int ii = 0; ii <clEnd.size(); ++ii) {
        std::cout<<clEnd[ii]<<" ";
      }
      std::cout<<std::endl;
      // clusters are all DS of the vertex
      if(clBeg.size() == 0) {
        ChkTopo1Vtx(allhits, hitcuts, tcl, vtx, iv, clEnd);
      } // clUSVtx.size
    } // iv

    if(fBEChgRat > 0.) cl2SetBeginEnd(allhits, tcl);

  } //RunCCHitFinder

/////////////////////////////////////////
  void CCHitRefinerAlg::ChkTopo1Vtx(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      CCHitFinderAlg::HitCuts& hitcuts,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl, 
      std::vector<ClusterCrawlerAlg::VtxStore>& vtx, unsigned int iv,
      std::vector<unsigned short>& clusters)
  {
    // check hits near a topo1 vertex = all clusters are DS of the vertex
    // find the US wire where there are no hits in this region. This is
    // done in case the vertex was not reconstructed properly
    unsigned short usWire = 9999;
    for(unsigned short wire = vtx[iv].Wire - 1; wire > vtx[iv].Wire - 6; --wire) {
      unsigned short index = wire - fFirstWire;
      // dead wire?
      if(WireHitRange[index].first == -1) continue;
      // no hits on the wire?
      if(WireHitRange[index].first == -2) {
        usWire = wire + 1;
        break;
      }
      unsigned short firsthit = WireHitRange[index].first;
      unsigned short lasthit = WireHitRange[index].second;
      bool WireSigOK = false;
      for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
        float tdiff = fabs(vtx[iv].Time - allhits[khit].Time);
        if (tdiff < 2.5 * allhits[khit].RMS) {
          // found a signal. Skip checking on this wire
          WireSigOK = true;
          break;
        } // tdiff test
      } // khit
      if(!WireSigOK) {
        usWire = wire + 1;
        break;
      }
    } // wire
    std::cout<<"usWire "<<usWire<<std::endl;
    if(usWire == 9999) return;
    // find the DS wire where the clusters are well separated
//    unsigned short dsWire = 0;
    std::vector<unsigned short> clWireOK;
    for(unsigned short ii = 0; ii < clusters.size(); ++ii) {
      unsigned short icl = clusters[ii];
  std::cout<<"Cluster "<<icl<<" "<<tcl[icl].tclhits.size()<<std::endl;
      unsigned short ntest = 0;
      for(unsigned short jj = tcl[icl].tclhits.size() - 1; jj > 0; --jj) {
        unsigned short hit = tcl[icl].tclhits[jj];
  std::cout<<jj<<"  hit "<<hit<<" "<<allhits[hit].numHits<<std::endl;
        if(ntest > 10) break;
        ++ntest;
      }
    }
  } // ChkTopo1Vtx

/*
/////////////////////////////////////////
    void CCHitRefinerAlg::GetRAT(std::vector<CCHitFinderAlg::CCHit>& allhits,
      unsigned short inwire, unsigned short intime, std::vector<float>& rat,
      unsigned short& firsttick)
    {
      // gets a wire signal above threshold on wire inwire in the region
      // of intime. Returns the signal vector and the time tick of the
      // first entry of the vector
      
      // assume failure
      firsttick = -1;
      rat.clear();
      
      // find a hit on the wire to get the wire object
      short iht = WireHitRange[inwire - fFirstWire].first;
      if(iht < 0) return;
      art::Ptr< recob::Wire> wire = allhits[iht].Wire();
      std::vector<float> signal( wire->Signal() );
      if(intime > (short)signal.size()) return;
      
      // look for the tick where the signal falls below a threshold (= 1)
      for(short it = intime; it > 0; it--) {
        if(signal[it] < 1) {
          firsttick = it;
          break;
        }
      }
      firsttick++;
      // fill the vector
      for(short it = firsttick; it < (short)signal.size(); it++) {
        rat.push_back(signal[it]);
        if(signal[it] < 1) break;
      }
      return;
    }

*/
////////////////////////////////////////////
  void CCHitRefinerAlg::FitNG(unsigned short nGaus, unsigned short tstart, 
      unsigned short tend, std::vector<float>& signal, 
      CCHitFinderAlg::HitCuts& hitcuts)
  {
    // fit the Region Above Threshold (rat) to nGaus Gaussians. The parameters, bounds
    // and fixed parameters should be defined by the calling routine

    fitChiDOF = 9999.;

    // define the fit string to pass to TF1
    std::stringstream numConv;
    std::string eqn = "gaus";
    if(nGaus > 1) eqn = "gaus(0)";
    for(unsigned short ii = 3; ii < nGaus*3; ii+=3){
      eqn.append(" + gaus(");
      numConv.str("");
      numConv << ii;
      eqn.append(numConv.str());
      eqn.append(")");
    }

    // load the data
    unsigned short npt = tend - tstart + 1;
    float *ticks = new float[npt];
    float *rat = new float[npt];
    unsigned short ii = 0;
    for(unsigned short tick = tstart; tick <= tend; ++tick) {
      ticks[ii] = (float)ii;
      rat[ii] = signal[tick];
      ++ii;
    }

    TGraph *fitn = new TGraph(npt, ticks, rat);
    TF1 *Gn = new TF1("gn",eqn.c_str());
    
    unsigned short npar = 3 * nGaus;
    
    if(par.size() != npar) return;
    if(parlo.size() != npar) return;
    if(parhi.size() != npar) return;

    // load the parameters
    for(unsigned short ipar = 0; ipar < npar; ++ipar) {
      Gn->SetParameter(ipar, par[ipar]);
      Gn->SetParLimits(ipar, parlo[ipar], parhi[ipar]);
    }

    // W = set weights to 1, N = no drawing or storing, Q = quiet
    // B = use parameters specified above
    fitn->Fit(Gn,"WNQB");
    
    // load the fit into a temp vector
    std::vector<double> partmp;
    std::vector<double> partmperr;

    for(unsigned short ipar = 0; ipar < npar; ++ipar) {
      partmp.push_back(Gn->GetParameter(ipar));
      partmperr.push_back(Gn->GetParError(ipar));
    }
    float ndof = npt - npar - 1;
    float chinorm = hitcuts.ChiNorms[plane];
    fitChiDOF = Gn->GetChisquare() / ( ndof * chinorm);
    
    delete ticks;
    delete rat;
    delete fitn;
    delete Gn;
  }
  
/////////////////////////////////////////
    void CCHitRefinerAlg::cl2SetBeginEnd(
        std::vector<CCHitFinderAlg::CCHit>& allhits,
        std::vector<ClusterCrawlerAlg::ClusterStore>& tcl)
    {
      // This routine prepares the clusters in tcl for stuffing into
      // recob::cluster. The start and end wire, time and slope are
      // defined based on the ratio of start and end charge. Tracks are
      // assumed to be going "downstream" => increasing wire number =>
      // from the end to the beginning of the tcl hit vector, unless the
      // charge ratio is significantly different
      
      for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
        // ignore deleted clusters
        if(tcl[ii].ID < 0) continue;
        bool GoingDS = true;
        if(tcl[ii].EndChg > fBEChgRat * tcl[ii].BeginChg) GoingDS = false;
        // need to swap the beginning/end?
        if(GoingDS) {
          // slope
          float tmp = tcl[ii].BeginSlp;
          tcl[ii].BeginSlp = tcl[ii].EndSlp;
          tcl[ii].EndSlp = tmp;
          // wire
          unsigned short itmp = tcl[ii].BeginWir;
          tcl[ii].BeginWir = tcl[ii].EndWir;
          tcl[ii].BeginWir = itmp;
          // time
          tmp = tcl[ii].BeginTim;
          tcl[ii].BeginTim = tcl[ii].EndTim;
          tcl[ii].EndTim = tmp;
          // charge
          tmp = tcl[ii].BeginChg;
          tcl[ii].BeginChg = tcl[ii].EndChg;
          tcl[ii].EndChg = tmp;
          // vertex
          short jtmp = tcl[ii].BeginVtx;
          tcl[ii].BeginVtx = tcl[ii].EndVtx;
          tcl[ii].EndVtx = jtmp;
        }
      }
    } // cl2SetBeginEnd



} // namespace cluster

