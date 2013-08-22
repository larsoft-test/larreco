/////////////////////////////////////////////////////////////////////
///
/// ClusterCrawlerAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithms for crawling along a string of hits to make line clusters
///
////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <stdint.h>
#include <iostream>
#include <iomanip>

#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "RecoAlg/ClusterCrawlerAlg.h"
#include "Filters/ChannelFilter.h"
#include "RawData/RawDigit.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Cluster.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

#include "TLinearFitter.h"
#include "TF1.h"

namespace cluster {
//------------------------------------------------------------------------------
  ClusterCrawlerAlg::ClusterCrawlerAlg(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

//------------------------------------------------------------------------------
  ClusterCrawlerAlg::~ClusterCrawlerAlg()
  {
    delete xwir;
    delete ytim;
    delete ytimerr;
    delete ychg;
    delete ychgerr;
  }

  void ClusterCrawlerAlg::reconfigure(fhicl::ParameterSet const& pset)
  { 
    fNumPass            = pset.get<             short  >("NumPass");
    fMaxHitsFit         = pset.get< std::vector<short> >("MaxHitsFit");
    fMinHits            = pset.get< std::vector<short> >("MinHits");
    fNHitsAve           = pset.get< std::vector<short> >("NHitsAve");
    fChgCut             = pset.get< std::vector<float> >("ChgCut");
    fWidCut             = pset.get< std::vector<float> >("WidCut");
    fChiCut             = pset.get< std::vector<float> >("ChiCut");
    fMaxWirSkip         = pset.get< std::vector<short> >("MaxWirSkip");
    fMinWirAfterSkip    = pset.get< std::vector<short> >("MinWirAfterSkip");
    fKinkChiRat         = pset.get< std::vector<float> >("KinkChiRat");
    fKinkAngCut         = pset.get< std::vector<float> >("KinkAngCut");
    fDoMerge            = pset.get< std::vector<bool>  >("DoMerge");
    fTimeDelta          = pset.get< std::vector<float> >("TimeDelta");
    fTimeDeltaLA        = pset.get< std::vector<float> >("TimeDeltaLA");

    fHitErrFac          = pset.get<             float  >("HitErrFac");
    fHitWidFac          = pset.get<             float  >("HitWidFac");
    fBEChgRat           = pset.get<             float  >("BEChgRat");
    fFudgeBigHits       = pset.get<             float  >("FudgeBigHits");
    fPairAngCut         = pset.get<             float  >("PairAngCut");
    fCurlyMergeAngCut   = pset.get<             float  >("CurlyMergeAngCut");
    fDoVertex           = pset.get<             bool   >("DoVertex");
    fDebugWire          = pset.get<             short  >("DebugWire");
    fDebugHit           = pset.get<             short  >("DebugHit");

    // use the ROOT linear fitter. Define 1st order polynomial w 2 params
    TLinearFitter *lf =  new TLinearFitter(2);
    lf->SetFormula("pol1");
    // determine the length of the TLinearFitter arrays
    short maxhits = 2;
    for(short ii = 0; ii < fNumPass; ii++) {
      if(fMaxHitsFit[ii] > maxhits) maxhits = fMaxHitsFit[ii];
    }
    // fit arrays for time and charge
    xwir = new Double_t[maxhits];
    ytim = new Double_t[maxhits];
    ytimerr = new Double_t[maxhits];
    ychg = new Double_t[maxhits];
    ychgerr = new Double_t[maxhits];

  }

  void ClusterCrawlerAlg::beginJob()
  {
  }

  // used for sorting hits on wires
  bool SortByLowHit(short i, short j) {return ((i > j));}
  // used for sorting clusters by length
  typedef std::pair<unsigned int, unsigned int> mypair;
  bool SortByLen(const mypair& L, const mypair& R) {return (L.first > R.first);}


  void ClusterCrawlerAlg::RunCrawler(art::PtrVector<recob::Hit>& plnhits,int plane,
     std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx)
  {
    // The hit collection is assumed to be sorted in increasing wire order

    if(plnhits.size() <= 2) return;

    // get the scale factor to convert dTick/dWire to dX/dU. This is used
    // to make the kink and merging cuts
    uint32_t channel = plnhits[0]->Wire()->RawDigit()->Channel();
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;
    float wirePitch = geom->WirePitch(geom->View(channel));
    float tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
    tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
    ScaleF = tickToDist / wirePitch;

    art::PtrVector<recob::Hit>::const_iterator it = plnhits.begin();
    fFirstWire = (*it)->WireID().Wire;
    it = plnhits.end()-1;
    fLastWire = (*it)->WireID().Wire;
    
    WireHitRange.clear();
    hiterr2.clear();
    hitwid.clear();
    NClusters = 0;
    
    prt = (fDebugWire > 0);
    
    // initialize WireHitRange to "no hits on wire" condition
    for(short wire = fFirstWire; wire <= fLastWire; wire++) {
      WireHitRange.push_back(std::make_pair(-2, -2));
    }

    // find dead wires in this region
    filter::ChannelFilter cf;
    for(short wire = fFirstWire+1; wire < fLastWire; wire++) {
      unsigned int pln = plane;
      unsigned int wir = wire;
      uint32_t chan = geom->PlaneWireToChannel(pln,wir);
      // remember to offset references to WireHitRange by the FirstWire
      if(cf.BadChannel(chan)) WireHitRange[wire - fFirstWire] = std::make_pair(-1, -1);
    }
    
    short lastwire = -1;
    short thishit = 0;
    short lastfirsthit = 0;
    if(prt) std::cout<<" W:H Time Wid Err Chg  Chi Mult"<<std::endl;
    for(art::PtrVector<recob::Hit>::const_iterator hitIter = plnhits.begin();
          hitIter != plnhits.end(); ++hitIter) {
      short thiswire = (*hitIter)->WireID().Wire;
      hitwid[thishit] = (*hitIter)->EndTime()-(*hitIter)->PeakTime();
      // calculate the hit position uncertainty
      float arg = fHitErrFac * hitwid[thishit];
      // modify the hit width
      hitwid[thishit] = fHitWidFac * hitwid[thishit];
      hiterr2[thishit] = arg * arg;
  if(prt) {
  std::cout<<thiswire<<":"<<thishit<<" "<<(int)(*hitIter)->PeakTime();
  std::cout<<" "<<std::setprecision(2)<<hitwid[thishit];
  std::cout<<" "<<std::setprecision(2)<<arg;
  std::cout<<" "<<(int)(*hitIter)->Charge();
  std::cout<<" "<<std::setprecision(2)<<(*hitIter)->GoodnessOfFit();
  std::cout<<" "<<(*hitIter)->Multiplicity();
  std::cout<<std::endl;
  }
      if(thiswire > lastwire) {
        WireHitRange[lastwire - fFirstWire] = std::make_pair(lastfirsthit,thishit);
        lastwire = thiswire;
        lastfirsthit = thishit;
      } else if(thiswire < lastwire) {
        mf::LogError("ClusterCrawler")<<"ERROR: Hits not sorted!!";
        return;
      }
      thishit++;
    }
    WireHitRange[fLastWire - fFirstWire] = std::make_pair(lastfirsthit, plnhits.size());
    
    prt = false;
    unsigned int nHitsUsed = 0;
    bool AllDone = false;
    for(short thispass = 0; thispass < fNumPass; thispass++) {
      pass = thispass;
//      std::cout<<"******** ClusterCrawler pass "<<pass<<std::endl;
      // look for a starting cluster that spans a block of wires
      short span = 3;
      if(fMinHits[pass] < span) span = fMinHits[pass];
      for(short iwire = fLastWire; iwire > fFirstWire; iwire--) {
        short ifirsthit = WireHitRange[iwire - fFirstWire].first;
        // skip bad wires or no hits on the wire
        if(ifirsthit < 0) continue;
        short ilasthit = WireHitRange[iwire - fFirstWire].second;
        for(short ihit = ifirsthit; ihit < ilasthit; ihit++) {
          bool ClusterAdded = false;
          // skip used hits
          if(hiterr2[ihit] < 0) continue;
          // skip multiple hits
          if(plnhits[ihit]->Multiplicity() > 1) continue;
          for(short jwire = iwire - span + 1; jwire < iwire; jwire++) {
            short jfirsthit = WireHitRange[jwire - fFirstWire].first;
            if(jfirsthit < 0) continue;
            short jlasthit = WireHitRange[jwire - fFirstWire].second;
            for(short jhit = jfirsthit; jhit < jlasthit; jhit++) {
              if(hiterr2[jhit] < 0) continue;
              // skip multiple hits
              if(plnhits[jhit]->Multiplicity() > 1) continue;
              // start a cluster with these two hits
              fcl2hits.clear();
              fAveWid = -1.;
              fAveChg = -1.;
              clStopCode = -1;
              clProcCode = pass;
              clAssn = -1;
              
              fcl2hits.push_back(ihit);
              fcl2hits.push_back(jhit);
              // define the fit origin. Use the upstream hit
              wire0 = jwire;
              cl2Fit(plnhits);
              // kill it if something bad happened in the fitter
              if(clChisq > 5) {
                fcl2hits.clear();
                continue;
              }
              // now look for hits to add on the intervening wires
              bool SigOK = false;
              bool HitOK = false;
              bool clok = true;
              for(short kwire = jwire+1; kwire < iwire; kwire++) {
                cl2AddHit(plnhits, kwire, false, HitOK, SigOK);
                // no hit added and no nearby hit either
                if(!HitOK && !SigOK) {
                  clok = false;
                  break;
                }
              }
              // kill it?
              if((short)fcl2hits.size() < span || !clok) {
                fcl2hits.clear();
                continue;
              }
              // sort them by decreasing wire number
              // assume that this is the same as sorting by decreasing 
              // hit number. This only needs to be done on the starting cluster
              // hits will be added in the proper order by cl2Follow
              std::sort(fcl2hits.begin(), fcl2hits.end(), SortByLowHit);
              // re-fit
              // define the hit origin
              std::vector<short>::reverse_iterator ii = fcl2hits.rbegin();
              short jj = *ii;
              wire0 = plnhits[jj]->WireID().Wire;
              cl2Fit(plnhits);
              if(clChisq > 5) {
                fcl2hits.clear();
                continue;
              }
              // save the cluster begin info
              clBeginSlp = (float)clpar[1];
              clBeginSlpErr = (float)clparerr[1];
              // follow a trail of hits upstream
              cl2FollowUS(plnhits);
              if((short)fcl2hits.size() < fMinHits[pass]) {
                // is it long enough for the next pass?
                if(pass < fNumPass-1 && (short)fcl2hits.size() >= fMinHits[pass+1]) {
                  // try to follow further using next pass cuts
                  pass++;
                  cl2Fit(plnhits);
                  cl2FollowUS(plnhits);
                  pass--;
                  clEndSlp = (float)clpar[1]; // save the slope at the end
                  clEndSlpErr = (float)clparerr[1];
                  clEndChg = fAveChg;
                  // set a special code for later cleaning
                  clProcCode = pass + 2001;
                  cl2TmpStore(plnhits, tcl);
                  ClusterAdded = true;
                  nHitsUsed += fcl2hits.size();
                  AllDone = (nHitsUsed == plnhits.size());
                  break;
                } else {
                  // kill it
                  fcl2hits.clear();
                }
              } else {
                clEndSlp = (float)clpar[1]; // save the slope at the end
                clEndSlpErr = (float)clparerr[1];
                clEndChg = fAveChg;
                cl2TmpStore(plnhits, tcl); // store the cluster
                ClusterAdded = true;
                nHitsUsed += fcl2hits.size();
                AllDone = (nHitsUsed == plnhits.size());
                break;
              }
            } // jhit
            if(ClusterAdded || AllDone) break;
          } // jwire
          if(AllDone) break;
        } // ihit
        if(AllDone) break;
      } // iwire
      if(AllDone) break;
      if(pass > 0 && fDebugWire != 0 ) {
        std::cout<<"Clustering done in plane = "<<plane<<std::endl;
        cl2Print(plnhits, tcl);
      }
      // try to merge clusters 
      if(fDoMerge[pass]) cl2ChkMerge(plnhits, tcl);
      if(fDoMerge[pass] && pass > 0 && fDebugWire < 0) {
        std::cout<<"After merging in plane = "<<plane<<std::endl;
        cl2Print(plnhits, tcl);
      }
    } // pass
    
    // prepare close pair clusters for 3D matching
    cl2ChkPair(plnhits, tcl);
    
    // try to merge short curly clusters
    if(fCurlyMergeAngCut > 0.) cl2CurlyMerge(plnhits, tcl);
    short ncl = 0;
    for(unsigned int ii = 0; ii < tcl.size(); ii++) {
      if(tcl[ii].ID > 0) ncl++;
    }
    
    if(fDoVertex) cl2DoVertex(plnhits, tcl, vtx);
    
    // re-define the beginning and end of the cluster using the average charge
    // ratio if it is significant
    if(fBEChgRat > 0.) cl2SetBeginEnd(plnhits, tcl);
    
//    std::cout<<"nverts "<<vtx.size()<<std::endl;

    WireHitRange.clear();
    hiterr2.clear();
    hitwid.clear();
    
    return;
  } // RunCrawler

/*
/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2Clean(art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl)
    {
      // attaches hits that were missed during clustering to the ends of
      // isolated clusters. These hits may have been missed if the cluster
      // was found using cuts on one pass but would be included if the 
      // cluster was formed on a later pass
      
      for(unsigned int ii = 0; ii < tcl.size(); ii++) {
        // look at the End first
        if(tcl[ii].ID < 0) continue;
        // get the pass that
        short cpass = tcl[ii].ProcCode - 10 * (tcl[ii].ProcCode);
        // check for isolation
        bool isIsolated = true;
        for(unsigned int jj = 0; jj < tcl.size(); jj++) {
          short dw = abs(tcl[ii].EndWir - tcl[jj].EndWir);
          short dt = abs(tcl[ii].EndTim - tcl[jj].EndTim);
          if(dw < 4 && dt < 10) {
            isIsolated = false;
            break;
          }
        } // jj cluster loop
        if(isIsolated) {
        }
      } // ii cluster loop
      
      return;
    }
*/
/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2DoVertex(art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx)
    {
      // try to make 2D vertices
      
      if(tcl.size() < 2) return;
      
      // initialize the begin and end vertex IDs
      for(unsigned int ii = 0; ii < tcl.size(); ii++) {
        tcl[ii].BeginVtx = -99;
        tcl[ii].EndVtx = -99;
      }
      vtx.clear();

      art::ServiceHandle<geo::Geometry> geom;
      art::ServiceHandle<util::DetectorProperties> detprop;
      
      unsigned int plane = plnhits[0]->WireID().Plane;
      float nwires = geom->Nwires(plane);
      float maxtime = detprop->NumberTimeSamples();

  if(fDebugHit < 0) {
    std::cout<<"DoVertex plane "<<plane<<std::endl;
    cl2Print(plnhits, tcl);
  }
      for (unsigned int it1 = 0; it1 < tcl.size() - 1; it1++) {
        // ignore abandoned clusters
        if(tcl[it1].ID < 0) continue;
        // ignore already attached clusters
        if(tcl[it1].BeginVtx >= 0 && tcl[it1].EndVtx >= 0) continue;
        // ignore short clusters
        if(tcl[it1].tclhits.size() < 10) continue;
        float es1 = tcl[it1].EndSlp;
        short ew1 = tcl[it1].EndWir;
        float et1 = tcl[it1].EndTim;
        float bs1 = tcl[it1].BeginSlp;
        short bw1 = tcl[it1].BeginWir;
        float bt1 = tcl[it1].BeginTim;
        for (unsigned int it2 = it1 + 1; it2 < tcl.size(); it2++) {
          // ignore abandoned clusters
          if(tcl[it2].ID < 0) continue;
          // ignore already attached clusters
          if(tcl[it2].BeginVtx >= 0 && tcl[it2].EndVtx >= 0) continue;
          // try to attach cluster to existing vertices at either end
          cl2ClsVertex(plnhits, tcl, vtx, it2);
          // ignore short clusters
          if(tcl[it2].tclhits.size() < 10) continue;
          float es2 = tcl[it2].EndSlp;
          short ew2 = tcl[it2].EndWir;
          float et2 = tcl[it2].EndTim;
          float bs2 = tcl[it2].BeginSlp;
          short bw2 = tcl[it2].BeginWir;
          float bt2 = tcl[it2].BeginTim;
  if(fDebugHit < 0) std::cout<<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID<<std::endl;
  // topo 1: check for vtx US of the ends of both clusters
          if(tcl[it1].EndVtx < 0 && tcl[it2].EndVtx < 0) {
            float dsl = es2 - es1;
            if(fabs(dsl) < 0.001) dsl = 0.001;
            // find vertex wire and vertex time in float precision (fvw, fvt)
            float fvw = 0.5 + (et1 - ew1 * es1 - et2 + ew2 * es2) / dsl;
            if(fvw > 0. && fvw < nwires) {
              // vertex wire in the detector
              short vw = fvw;
              // require vtx in the range of wires with hits AND
              // vtx US of both clusters AND
              // vtx not too far US of both clusters
              if(vw >= fFirstWire && 
                 vw <= ew1 && vw <= ew2 &&
                 vw  > ew1 - 10 && vw  > ew2 - 10) {
                float fvt = et1 + (vw - ew1) * es1;
  if(fDebugHit < 0) std::cout<<"topo1 vtx wire "<<vw<<" time "<<(int)fvt<<std::endl;
                if(fvt > 0. && fvt < maxtime) {
                  // vertex wire US of cluster ends and time in the detector
                  // Check this against existing vertices and update
                  cl2ChkVertex(plnhits, tcl, vtx, vw, fvt, it1, it2, 1);
                } // fvt in detector
              } // vw topo 1 check
            } // fvw in detector
          } // topo 1
  // topo 2: check for vtx US of it1 and DS of it2
          if(tcl[it1].EndVtx < 0 && tcl[it2].BeginVtx < 0) {
            float dsl = bs2 - es1;
            if(fabs(dsl) < 0.001) dsl = 0.001;
            float fvw = 0.5 + (et1 - ew1 * es1 - bt2 + bw2 * bs2) / dsl;
            if(fvw > 0 && fvw < nwires) {
              // vertex wire in the detector
              short vw = fvw;
              if(vw <= ew1 && vw >= bw2) {
                float fvt = et1 + (vw - ew1) * es1;
  if(fDebugHit < 0) std::cout<<"topo2 vtx wire "<<vw<<" time "<<(int)fvt<<std::endl;
                if(fvt > 0. && fvt < maxtime) {
                  cl2ChkVertex(plnhits, tcl, vtx, vw, fvt, it1, it2, 2);
                } // fvt in detector
              } // vw topo 2 check
            } // fvw in detector
          } // topo 2
  // topo 3: check for vtx DS of it1 and US of it2
          if(tcl[it1].BeginVtx < 0 && tcl[it2].EndVtx < 0) {
            float dsl = bs1 - es2;
            if(fabs(dsl) < 0.001) dsl = 0.001;
            float fvw = 0.5 + (et2 - ew2 * es2 - bt1 + bw1 * bs1) / dsl;
            if(fvw > 0 && fvw < nwires) {
              short vw = fvw;
              if(vw <= ew2 && vw >= bw1) {
                float fvt = et2 + (vw - ew2) * es2;
  if(fDebugHit < 0) std::cout<<"topo3 vtx wire "<<vw<<" time "<<(int)fvt<<std::endl;
                if(fvt > 0. && fvt < maxtime) {
                  cl2ChkVertex(plnhits, tcl, vtx, vw, fvt, it1, it2, 3);
                } // fvt in detector
              } // vw topo 3 check
            } // fvw in detector
          } // topo 3
  // topo 4: check for vtx DS of it1 and DS of it2
          if(tcl[it1].BeginVtx < 0 && tcl[it2].BeginVtx < 0) {
            float dsl = bs2 - bs1;
            if(fabs(dsl) < 0.001) dsl = 0.001;
            // find vertex wire and vertex time in float precision (fvw, fvt)
            // convert to integer if within the detector (vw, vt)
            float fvw = 0.5 + (bt1 - bw1 * bs1 - bt2 + bw2 * bs2) / dsl;
            if(fvw > 0 && fvw < nwires) {
              // vertex wire in the detector
              short vw = fvw;
              // require vtx in the range of wires with hits AND
              // vtx DS of both clusters AND
              // vtx not too far DS of both clusters
              if(vw >= bw1 && 
                 vw >= bw2 && vw <= fLastWire &&
                 vw <  bw2 + 10 && vw <  bw1 + 10) {
                float fvt = bt1 + (vw - bw1) * bs1;
  if(fDebugHit < 0) std::cout<<"topo4 vtx wire "<<vw<<" time "<<(int)fvt<<std::endl;
                if(fvt > 0. && fvt < maxtime) {
                  // vertex wire US of cluster ends and time in the detector
                  // Check this against existing vertices and update
                  cl2ChkVertex(plnhits, tcl, vtx, vw, fvt, it1, it2, 4);
                } // fvt in detector
              } // vw topo 4 check
            } // fvw in detector
          } // it2
        } // topo 4
      } // it1
      
      // set the vertex weights
      for(unsigned int it = 0; it < tcl.size(); it++) {
        if(tcl[it].ID < 0) continue;
        // cluster weight = number of hits, saturated at 10
        float cw = tcl[it].tclhits.size();
        if(cw > 10) cw = 10;
        if(tcl[it].BeginVtx >=0) vtx[tcl[it].BeginVtx].Wght += cw;
        if(tcl[it].EndVtx >=0) vtx[tcl[it].EndVtx].Wght += cw;
      }
      

  if(fDebugHit < 0) {
    std::cout<<"Vertices "<<vtx.size()<<std::endl;
    for(unsigned int iv = 0; iv < vtx.size(); iv++) {
      std::cout<<"vtx "<<iv<<" wire "<<vtx[iv].Wire<<" time "<<(int)vtx[iv].Time<<" wght "<<(int)vtx[iv].Wght;
      std::cout<<" topo "<<vtx[iv].Topo<<std::endl;
    }
    cl2Print(plnhits, tcl);
  }
      return;
    }

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ClsVertex(art::PtrVector<recob::Hit>& plnhits, 
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
        unsigned int it)
    {
      // try to attach cluster it2 to an existing vertex      
      if(vtx.size() == 0) return;
      
      for(unsigned int iv = 0; iv < vtx.size(); iv++) {
        // determine which end to match - begin or end
        if(vtx[iv].Wire <= tcl[it].EndWir) {
          // project cluster to US vertex
          float tdiff = tcl[it].EndTim + (vtx[iv].Wire - tcl[it].EndWir) * tcl[it].EndSlp - vtx[iv].Time;
          if(fabs(tdiff) < 10) {
            bool SigOK = false;
            cl2ChkSignal(plnhits, vtx[iv].Wire, vtx[iv].Time, tcl[it].EndWir, tcl[it].EndTim, SigOK);
            if(SigOK) {
              // good match
              tcl[it].EndVtx = iv;
              return;
            } // SigOK
          } // tdiff
        } else if(vtx[iv].Wire >= tcl[it].BeginWir) {
          // project cluster to DS vertex
          float tdiff = tcl[it].BeginTim + (vtx[iv].Wire - tcl[it].BeginWir) * tcl[it].BeginSlp - vtx[iv].Time;
          if(fabs(tdiff) < 10) {
            bool SigOK = false;
            cl2ChkSignal(plnhits, vtx[iv].Wire, vtx[iv].Time, tcl[it].BeginWir, tcl[it].BeginTim, SigOK);
            if(SigOK) {
              // good match
              tcl[it].BeginVtx = iv;
              return;
            } // SigOK
          } // tdiff
        } // vtx[iv].Wire
      } // iv
      return;
    }



/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ChkVertex(art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
        short vw, float fvt, unsigned int it1, unsigned int it2, short topo)
      {
        // Checks the vertex (vw, vt) against the existing set of vertices.
        // If there a match, clusters it1 and/or it2 are associated with it
        // if there is signal between the existing vertex and the start of
        // the cluster. The topo flag indicates the type of vertex that was
        // found: 1 = US of it1 && US of it2, 2 = US of it1 && DS of it2,
        // 3 = DS of it1 && US of it2, 4 = DS of it1 and DS of it2.
        // didit is set true if a cluster is attached to a (new or old) vertex
        

  if(fDebugHit < 0) {
    std::cout<<"ChkVertex "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" topo "<<topo;
    std::cout<<" vw "<<vw<<" vt "<<(int)fvt<<std::endl;
  }
        // check vertex and clusters for proximity to existing vertices
        bool SigOK = false;
        for(unsigned int iv = 0; iv < vtx.size(); iv++) {
          if( abs( vw - vtx[iv].Wire) < 2 &&
             fabs(fvt - vtx[iv].Time) < 10) {
            // got a match. Check the appropriate cluster end and attach
            if( (topo == 1 || topo == 2) && tcl[it1].EndVtx < 0) {
              cl2ChkSignal(plnhits, vw, fvt, tcl[it1].EndWir, tcl[it1].EndTim, SigOK);
              if(SigOK) tcl[it1].EndVtx = iv;
  if(fDebugHit < 0)  std::cout<<"12 Attach cl "<<tcl[it1].ID<<" to vtx "<<iv<<std::endl;
            } else if( (topo == 3 || topo == 4) && tcl[it1].BeginVtx < 0) {
              cl2ChkSignal(plnhits, vw, fvt, tcl[it1].BeginWir, tcl[it1].BeginTim, SigOK);
              if(SigOK) tcl[it1].BeginVtx = iv;
  if(fDebugHit < 0)  std::cout<<"34 Attach cl "<<tcl[it1].ID<<" to vtx "<<iv<<std::endl;
            } // cluster it1
            if( (topo == 1 || topo == 3) && tcl[it2].EndVtx < 0) {
              cl2ChkSignal(plnhits, vw, fvt, tcl[it2].EndWir, tcl[it2].EndTim, SigOK);
              if(SigOK) tcl[it2].EndVtx = iv;
  if(fDebugHit < 0) std::cout<<"13 Attach cl "<<tcl[it2].ID<<" to vtx "<<iv<<std::endl;
            } else if( (topo == 2 || topo == 4) && tcl[it2].BeginVtx < 0) {
              cl2ChkSignal(plnhits, vw, fvt, tcl[it2].BeginWir, tcl[it2].BeginTim, SigOK);
              if(SigOK) tcl[it2].BeginVtx = iv;
  if(fDebugHit < 0) std::cout<<"24 Attach cl "<<tcl[it2].ID<<" to vtx "<<iv<<std::endl;
            } // cluster it2
            return;
          } // matched vertex
        } // iv

        // no match to existing vertices. Ensure that there is a wire signal between
        // the vertex and the appropriate ends of the clusters
        bool Sig1OK = false;
        bool Sig2OK = false;
        if(topo == 1 || topo == 2) {
          cl2ChkSignal(plnhits, vw, fvt, tcl[it1].EndWir, tcl[it1].EndTim, Sig1OK);
        } else {
          cl2ChkSignal(plnhits, vw, fvt, tcl[it1].BeginWir, tcl[it1].BeginTim, Sig1OK);
        }
        if(topo == 1 || topo == 3) {
          cl2ChkSignal(plnhits, vw, fvt, tcl[it2].EndWir, tcl[it2].EndTim, Sig2OK);
        } else {
          cl2ChkSignal(plnhits, vw, fvt, tcl[it2].BeginWir, tcl[it2].BeginTim, Sig2OK);
        }
  if(fDebugHit < 0) {
    std::cout<<"Chk new "<<tcl[it1].ID<<" OK "<<Sig1OK<<" "<<tcl[it2].ID<<" OK "<<Sig2OK;
    std::cout<<" Vtx at "<<vw<<" "<<(int)fvt<<std::endl;
  }
        // both clusters must have an OK signal
        if(Sig1OK && Sig2OK) {
          VtxStore newvx;
          newvx.Wire = vw;
          newvx.Time = fvt;
          newvx.Topo = topo;
          newvx.Wght = 0;
          vtx.push_back(newvx);
          short iv = vtx.size() - 1;
          if(topo == 1 || topo == 2) {
            tcl[it1].EndVtx = iv;
          } else {
            tcl[it1].BeginVtx = iv;
          }
          if(topo == 1 || topo == 3) {
            tcl[it2].EndVtx = iv;
          } else {
            tcl[it2].BeginVtx = iv;
          }
  if(fDebugHit < 0) std::cout<<"New vtx "<<iv<<" topo "<<topo<<" cls "<<tcl[it1].ID<<" "<<tcl[it2].ID<<std::endl;
        }

        return;
      }

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ChkSignal(art::PtrVector<recob::Hit>& plnhits,
      short wire1, float time1, short wire2, float time2, bool& SigOK)
    {
      // returns SigOK true if there is a signal on the line between
      // (wire1, time1) and (wire2, time2).
      SigOK = false;
      // get the begin and end right
      short wireb = wire1;
      float timeb = time1;
      short wiree = wire2;
      float timee = time2;
      // swap them?
      if(wiree > wireb) {
        wireb = wire2;
        timeb = time2;
        wiree = wire1;
        timee = time1;
      }
      wire0 = wiree;
      // checking a single wire?
      if(wireb == wiree) {
        clpar[1] = 0.;
      } else {
        clpar[1] = (timeb - timee) / (wireb - wiree);
      }
      for(short wire = wiree; wire < wireb + 1; wire++) {
        // assume there is no signal on this wire
        bool WireSigOK = false;
        float prtime = timee + (wire - wire0) * clpar[1];
        // skip dead wires
        short firsthit = WireHitRange[wire - fFirstWire].first;
        if(firsthit < -1) continue;
        short lasthit = WireHitRange[wire - fFirstWire].second;
//  std::cout<<"ChkSignal "<<wiree<<" "<<wire<<" "<<wireb<<" "<<(int)prtime;
//  std::cout<<" first "<<firsthit<<" last "<<lasthit;
        for(short khit = firsthit; khit < lasthit; khit++) {
          float tdiff = fabs(prtime - plnhits[khit]->PeakTime());
          if (tdiff < 2 * hitwid[khit]) {
            // found a signal. Skip checking on this wire
            WireSigOK = true;
            break;
          }
        } // khit
//        std::cout<<" SigOK "<<WireSigOK<<std::endl;
        if(!WireSigOK) return;
      } // wire
      SigOK = true;
      return;
    }

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2CurlyMerge(art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl)
    {
      // try to merge short curly clusters. The algorithm uses the fact that
      // upstream clusters are added to the tcl struct after downstream 
      // clusters.
      
      if(tcl.size() < 2) return;
      
      short tclsize = tcl.size();
      for(short it1 = 0; it1 < tclsize - 1; it1++) {
        if(tcl[it1].ID < 0) continue;
        // ignore long clusters
        if(tcl[it1].tclhits.size() > 20) continue;
        // ignore close pair clusters
        if(tcl[it1].ProcCode == 1000) continue;
        // ignore longish straight tracks
        short ew1 = tcl[it1].EndWir;
        short et1 = tcl[it1].EndTim;
        // we will look for a 2nd cluster on the next upstream wire = pw1
        short pw1 = ew1 - 1;
        // ignore dead wires for now
        if(WireHitRange[pw1 - fFirstWire].first == -1) continue;
        // project the time to this wire using the two end hits of cluster 1
        std::vector<short>::reverse_iterator i4 = tcl[it1].tclhits.rbegin() + 1;
        short ih4 = *i4;
        float dt34 = plnhits[ih4]->PeakTime() - et1;
        float dw34 = plnhits[ih4]->WireID().Wire - ew1;
        // slope dT/dW
        float sl34 = dt34 / dw34;
        float th34 = atan(ScaleF * sl34);
        // projected time to the next wire (dW = 1)
        float pt1 = et1 - sl34;
        // use the fit array to hold the new(?) cluster
        fcl2hits.clear();
        fcl2hits = tcl[it1].tclhits;
        clBeginSlp = tcl[it1].BeginSlp;
        bool didit = false;
  if(prt) std::cout<<"Curly cl1 "<<tcl[it1].ID<<" pt1 "<<pt1<<std::endl;
        for(short it2 = it1 + 1; it2 < tclsize; it2++) {
          if(tcl[it2].ID < 0) continue;
          if(tcl[it2].tclhits.size() > 20) continue;
          // ignore close pair clusters
          if(tcl[it2].ProcCode == 1000) continue;
          short bw2 = tcl[it2].BeginWir;
          // cluster 2 begins on wire US of the end of cluster 1?
          if(bw2 != pw1) continue;
          short bt2 = tcl[it2].BeginTim;
  if(prt) std::cout<<" >> cl2 "<<tcl[it2].ID<<" dt "<<fabs(bt2 - pt1)<<std::endl;
          // make a rough cut on the time difference
          if(fabs(bt2 - pt1) > 30) continue;
          didit = true;
          // find the begin angle of this cluster
          std::vector<short>::iterator i2 = tcl[it2].tclhits.begin() + 1;
          short ih2 = *i2;
          float dt12 = bt2 - plnhits[ih2]->PeakTime();
          float dw12 = bw2 - plnhits[ih2]->WireID().Wire;
          float sl12 = dt12 / dw12;
          float th12 = atan(ScaleF * sl12);
          float dth = fabs(th34 - th12);
  if(prt) std::cout<<" >> dth "<<dth<<" cut "<<fCurlyMergeAngCut<<std::endl;
          if(dth < fCurlyMergeAngCut) {
            // add it to fcl2hits
            for(std::vector<short>::iterator iht2 = tcl[it2].tclhits.begin();
                iht2 != tcl[it2].tclhits.end(); ++iht2) {
              fcl2hits.push_back(*iht2);
            }
            // declare cluster 2 obsolete
            tcl[it2].ID = -tcl[it2].ID;
            clStopCode = tcl[it2].StopCode;
            // update the parameters of end of the cluster
            ew1 = tcl[it2].EndWir;
            et1 = tcl[it2].EndTim;
            // we will look for a 2nd cluster on the next upstream wire = pw1
            pw1 = ew1 - 1;
            // ignore dead wires for now
            if(WireHitRange[pw1 - fFirstWire].first == - 1) continue;
            // project the time to this wire using the two end hits of cluster 1
            i4 = tcl[it2].tclhits.rbegin() + 1;
            ih4 = *i4;
            dt34 = plnhits[ih4]->PeakTime() - et1;
            dw34 = plnhits[ih4]->WireID().Wire - ew1;
            // slope dT/dW
            clEndSlp = dt34 / dw34;
            th34 = atan(ScaleF * clEndSlp);
            // projected time to the next wire (dW = 1)
            pt1 = et1 - sl34;
          }
        } // it2
        if(didit) {
  if(prt) std::cout<<"New cluster "<<std::endl;
          // store the new cluster
          clProcCode = 10000;
          cl2TmpStore(plnhits, tcl);
        }
      } // it1
  prt = false;
      return;
    } // cl2CurlyMerge
  
/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2SetBeginEnd(art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl)
    {
      // This routine prepares the clusters in tcl for stuffing into
      // recob::cluster. The start and end wire, time and slope are
      // defined based on the ratio of start and end charge. Tracks are
      // assumed to be going "downstream" => increasing wire number =>
      // from the end to the beginning of the tcl hit vector, unless the
      // charge ratio is significantly different
      
      for(unsigned int ii = 0; ii < tcl.size(); ii++) {
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
          short itmp = tcl[ii].BeginWir;
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
          // no obvious need to swap the hit order...
        }
      }
      return;
    } // cl2SetBeginEnd

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ChkPair(art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl)
    {
      // break clusters that are a part of a pair of clusters that
      // are less than PairAngCut
      
      if(tcl.size() < 2) return;
      // The size of the ClusterStore vector will increase after merging
      // is done.
      short tclsize = tcl.size();

      for (short it1 = 0; it1 < tclsize - 1; it1++) {
        // ignore abandoned clusters
        if(tcl[it1].ID < 0) continue;
        // ignore short clusters
        if(tcl[it1].tclhits.size() < 10) continue;
        float bs1 = tcl[it1].BeginSlp;
        // convert slope to angle
        float bth1 = atan(ScaleF * bs1);
        short bw1 = tcl[it1].BeginWir;
        float bt1 = tcl[it1].BeginTim;
        short ew1 = tcl[it1].EndWir;
        for (short it2 = it1 + 1; it2 < tclsize; it2++) {
          // ignore abandoned clusters
          if(tcl[it2].ID < 0) continue;
          if(tcl[it1].ID < 0) continue;
          // ignore short clusters
          if(tcl[it2].tclhits.size() < 10) continue;
          float bs2 = tcl[it2].BeginSlp;
          // convert slope to angle
          float bth2 = atan(ScaleF * bs2);
          short bw2 = tcl[it2].BeginWir;
          float bt2 = tcl[it2].BeginTim;
          short ew2 = tcl[it2].EndWir;
          if(fabs(bth1 - bth2) < fPairAngCut) {
            // check for begin angle difference
            if(ew2 > ew1 && ew2 < bw1) {
              // cluster 2 end is in the wire bounds of cluster 1
              // find the vertex position using the begin slope of both
              // clusters. 
              float dsl = bs2 - bs1;
              if(fabs(dsl) < 0.001) continue;
//  std::cout<<"Pair1 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<std::endl;
              short vw = (int)(0.5 + (bt1 - bw1 * bs1 - bt2 + bw2 * bs2) / dsl);
//  std::cout<<"Vtx wire "<<vw<<std::endl;
              // compare the vertex wire with the end of cluster 1
              if(abs(vw - ew1) < 5) cl2DoSplit(plnhits, tcl, it1, it2);
            } else if(ew1 > ew2 && ew1 < bw2) {
              // cluster 1 end is in the wire bounds of cluster 2
              // find the vertex position using the begin slope of both
              // clusters. 
              float dsl = bs1 - bs2;
              if(fabs(dsl) < 0.001) continue;
//  std::cout<<"Pair2 "<<tcl[it2].ID<<" "<<tcl[it1].ID<<std::endl;
              short vw = (int)(0.5 + (bt2 - bw2 * bs2 - bt1 + bw1 * bs1) / dsl);
//  std::cout<<"Vtx wire "<<vw<<std::endl;
              // compare the vertex wire with the end of cluster 2
              if(abs(vw - ew2) < 5) cl2DoSplit(plnhits, tcl, it2, it1);
            } // begin/end wire test
          } // PairAngCut test
        } // it2
      } // it1

      return;
    }

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2DoSplit(art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl, unsigned int it1, unsigned int it2)
    {
      // split cluster 1 at the gap in the middle of it. If no gap is found
      // split the cluster at the start of cluster 2

//  std::cout<<"DoSplit cl1 "<<tcl[it1].ID<<" cl2 "<<tcl[it2].ID<<std::endl;

      short splitwire = -1;
      short lastwire = -1;
      for(std::vector<short>::iterator iht1 = tcl[it1].tclhits.begin();
          iht1 != tcl[it1].tclhits.end(); ++iht1) {
        short hit = *iht1;
        short wire = plnhits[hit]->WireID().Wire;
        if(lastwire < 0) {
          lastwire = wire;
        } else {
          // require a significant gap
          if(wire < lastwire - 3) {
            splitwire = lastwire;
            break;
          }
          lastwire = wire;
        } // lastwire test
      } // cluster 1 hits
//  std::cout<<"Split at wire "<<splitwire<<std::endl;
      // split at the end of cluster 2 if no gap in cluster 1 was found
      if(splitwire < 0) splitwire = tcl[it2].EndWir;

      // split up cluster 1 - the begin end first
      clBeginSlp = tcl[it1].BeginSlp;
      clBeginSlpErr = tcl[it1].BeginSlpErr;
      clBeginWir = tcl[it1].BeginWir;
      clBeginTim = tcl[it1].BeginTim;
      clBeginChg = tcl[it1].BeginChg;
      clStopCode = 5;
      clProcCode = tcl[it1].ProcCode + 1000;
      fcl2hits.clear();
      for(std::vector<short>::iterator iht1 = tcl[it1].tclhits.begin();
          iht1 != tcl[it1].tclhits.end(); ++iht1) {
        short hit = *iht1;
        fcl2hits.push_back(hit);
        short wire = plnhits[hit]->WireID().Wire;
        if(wire == splitwire) {
          clEndTim = plnhits[hit]->PeakTime();
          break;
        }
      }
      // re-fit the end of the cluster
      wire0 = splitwire;
      cl2Fit(plnhits);
      clEndSlp = (float)clpar[1];
      clEndSlpErr = (float)clparerr[1];
      clEndWir = splitwire;
      clEndChg = fAveChg;
      cl2TmpStore(plnhits, tcl);
      // index of the first DownStream cluster
      unsigned int DScl1 = tcl.size() - 1;
      
      // look for an associated cluster with cluster1
      short jj = -1;
      for(unsigned int ii = 0; ii < tcl.size(); ii++) {
        if(tcl[ii].Assn == (int)it1) {
          jj = ii;
          break;
        }
      }
      unsigned int it3 = 9999;
//  std::cout<<"Assoc cluster? "<<jj<<std::endl;
      
      if(jj >= 0) {
        it3 = jj;
        // cluster 1 has an associated cluster embedded within it.
        // attach the hits on cluster 1 to the embedded cluster
        clStopCode = tcl[it3].StopCode;
        clProcCode = tcl[it3].ProcCode + 1000;
        fcl2hits = tcl[it3].tclhits;
        short wire = 0;
        for(std::vector<short>::iterator iht1 = tcl[it1].tclhits.begin();
          iht1 != tcl[it1].tclhits.end(); ++iht1) {
          short hit = *iht1;
          wire = plnhits[hit]->WireID().Wire;
          if(wire < splitwire) fcl2hits.push_back(hit);
        }
        wire0 = wire;
        cl2Fit(plnhits);
        clEndSlp = (float)clpar[1]; // save the slope at the end
        clEndSlpErr = (float)clparerr[1];
        clEndChg = fAveChg;
        cl2TmpStore(plnhits, tcl);
      } else {
        // create a new cluster using the hits on cluster 1 that are
        // upstream of splitwire. 
        clStopCode = tcl[it1].StopCode;
        clProcCode = tcl[it1].ProcCode + 1000;
        fcl2hits.clear();
        bool didfit = false;
        short wire = 0;
        for(std::vector<short>::iterator iht1 = tcl[it1].tclhits.begin();
          iht1 != tcl[it1].tclhits.end(); ++iht1) {
          short hit = *iht1;
          wire = plnhits[hit]->WireID().Wire;
          if(wire < splitwire) {
            fcl2hits.push_back(hit);
            // re-fit the begin end when there are enough hits
            if(fcl2hits.size() == 3) {
//  std::cout<<"refit begin "<<tcl[it1].ID<<std::endl;
              didfit = true;
              wire0 = wire;
              cl2Fit(plnhits);
              clBeginSlp = (float)clpar[1];
              clBeginSlpErr = (float)clparerr[1];
            } // fcl2hits test
          } // wire < splitwire
        } // iht1 iterator
        if(didfit) {
          // did a fit at the beginning of the cluster
          if(fcl2hits.size() > 4) {
            // long cluster: re-fit the end
            wire0 = wire;
            cl2Fit(plnhits);
            clEndSlp = (float)clpar[1];
            clEndSlpErr = (float)clparerr[1];
          } else {
            // short cluster: set the end params to the begin params
            clEndSlp = clBeginSlp;
            clEndSlpErr = clBeginSlpErr;
          }
          cl2TmpStore(plnhits, tcl);
        } else {
          // didn't do the fit for some reason.
          std::cout<<"Coding error in cl2DoSplit! "<<std::endl;
        } // didfit test
      } // tcl[it1].Assn > 0 else if test
      
      // At this point we have 3 clusters. The most upstream cluster
      // is the last one added to tcl
      unsigned int UScl = tcl.size() - 1;
      // The first downstream cluster was defined above
      // the second downstream cluster is it2
      unsigned int DScl2 = it2;

//  std::cout<<"UScl "<<UScl<<" DScl1 "<<DScl1<<" DScl2 "<<DScl2<<std::endl;

      // mark clusters it1 and it3 as obsolete and set the processor code to
      // show that this routine clobbered them
      tcl[it1].ID = -tcl[it1].ID;
      tcl[it1].ProcCode += 1000;
      if(it3 < tcl.size()) {
        tcl[it3].ID = -tcl[it3].ID;
        tcl[it3].ProcCode += 1000;
      }
      // make the new association. This will turn into a vertex later (maybe)
      tcl[DScl1].Assn = UScl;
      tcl[DScl2].Assn = UScl;
      
      return;
    }
/*
/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2SortByLength(std::vector<ClusterStore>& tcl,
        std::vector<unsigned int>& sortindex)
    {
      // sorts the temporary cluster vector by decreasing number of hits,
      // while ignoring abandoned clusters. Returns index vector with the
      // sort order
      
      // form a vector of pairs of the number of hits and the index
      std::vector< std::pair<unsigned int, unsigned int> > index;
      for(unsigned int ii = 0; ii < tcl.size(); ii++) {
        if(tcl[ii].ID > 0) index.push_back(std::make_pair(tcl[ii].tclhits.size(),ii));
      }
      std::sort(index.begin(), index.end(), SortByLen);
      sortindex.clear();
      for(unsigned int ii = 0; ii < tcl.size(); ii++) {
        sortindex.push_back(index[ii].second);
//    std::cout<<"sortd "<<index[ii].first<<" "<<index[ii].second<<std::endl;
      }
      return; 
    }
*/
/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ChkMerge(art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl)
    {
      // Try to merge clusters. Clusters that have been subsumed in other
      // clusters, i.e. no longer valid, have ID < 0
      
      if(tcl.size() < 2) return;
      // The size of the ClusterStore vector will increase while merging
      // is done.
      
      unsigned int tclsize = tcl.size();
            
      for(unsigned int it1 = 0; it1 < tclsize - 1; it1++) {
        // ignore already merged clusters
        if(tcl[it1].ID < 0) continue;
        float bs1 = tcl[it1].BeginSlp;
        // convert slope to angle
        float arg = ScaleF * bs1;
        float bth1 = atan(arg);
        // error on the angle
        float bth1e = ScaleF * tcl[it1].BeginSlpErr / (1 + arg * arg);
        // more compact notation for begin/end, wire/time/chg/slp/theta, 1/2
        short bw1 = tcl[it1].BeginWir;
        float bt1 = tcl[it1].BeginTim;
        float bc1 = tcl[it1].BeginChg;
        float es1 = tcl[it1].EndSlp;
        // convert slope to angle
        arg = ScaleF * es1;
        float eth1 = atan(arg);
        float eth1e =  ScaleF * tcl[it1].EndSlpErr / (1 + arg * arg);
        short ew1 = tcl[it1].EndWir;
        float et1 = tcl[it1].EndTim;
        float ec1 = tcl[it1].EndChg;
        short pass1 = tcl[it1].ProcCode - 10 * (tcl[it1].ProcCode / 10);
        for(unsigned int it2 = it1 + 1; it2 < tclsize; it2++) {
          // ignore already merged clusters
          if(tcl[it1].ID < 0) continue;
          if(tcl[it2].ID < 0) continue;
          float bs2 = tcl[it2].BeginSlp;
          float bs2e = tcl[it2].BeginSlpErr;
          // convert slope to angle
          arg = ScaleF * bs2;
          float bth2 = atan(arg);
          // error on the angle
          float bth2e = ScaleF * tcl[it2].BeginSlpErr / (1 + arg * arg);
          short bw2 = tcl[it2].BeginWir;
          float bt2 = tcl[it2].BeginTim;
          float bc2 = tcl[it2].BeginChg;
          float es2 = tcl[it2].EndSlp;
          arg = ScaleF * es2;
          float eth2 = atan(arg);
          float eth2e = ScaleF * tcl[it2].EndSlpErr / (1 + arg * arg);
          short ew2 = tcl[it2].EndWir;
          float et2 = tcl[it2].EndTim;
          float ec2 = tcl[it2].EndChg;
          short pass2 = tcl[it2].ProcCode - 10 * (tcl[it2].ProcCode / 10);
          // use the more promiscuous pass for cuts
          float angcut = fKinkAngCut[pass1];
          if(fKinkAngCut[pass2] > angcut) angcut = fKinkAngCut[pass2];
          short skipcut = fMaxWirSkip[pass1];
          if(fMaxWirSkip[pass2] > skipcut) skipcut = fMaxWirSkip[pass2];
          float chgcut = fChgCut[pass1];
          if(fChgCut[pass2] > chgcut) chgcut = fChgCut[pass2];
          float timecut = fTimeDelta[pass];
          if(fTimeDelta[pass2] > timecut) timecut = fTimeDelta[pass2];
          
          float dth = fabs(bth2 - eth1);
//  if(fDebugWire < 0 && bw2 < ew1 && (ew1 - bw2)  < skipcut) {
  if(fDebugWire < 0 && bw2 < ew1 ) {
    std::cout<<"Chk1 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" dW "<<(ew1 - bw2);
    std::cout<<" skipcut "<<skipcut;
    std::cout<<" dth "<<dth<<" therrs "<<bth2e<<" "<<eth1e<<" angcut "<<angcut<<std::endl;
  }
          if(bw2 < ew1 && (ew1 - bw2)  < skipcut && dth < angcut) {
            // look for US and DS broken clusters w similar angle.
            // US cluster 2 merge with DS cluster 1?
            // This is the most likely occurrence given the order in which
            // clusters are created so put it first.
            float chgrat = 2 * fabs((bc2 - ec1) / (bc2 + ec1));
            // project bw2,bt2 to ew1
            float dtim = fabs(bt2 + (ew1-bw2)*bs2 - et1);
            // error on the projection
            float dtime = fabs(bt2 + (ew1-bw2)*bs2e - et1);
  if(fDebugWire < 0) {
    std::cout<<" dtim "<<dtim<<" err "<<dtime<<" timecut "<<(int)timecut;
    std::cout<<" chgrat "<<chgrat<<" chgcut "<<chgcut;
    std::cout<<std::endl;
  }
            if(chgrat < chgcut && dtim < timecut) {
              cl2DoMerge(plnhits, tcl, it2, it1, 10);
              tclsize = tcl.size();
              break;
            }
          } // US cluster 2 merge with DS cluster 1?
          
          dth = fabs(bth1 - eth2);
  if(fDebugWire < 0 && bw1 < ew2 && (ew2 - bw1)  < skipcut) {
    std::cout<<"Chk2 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" dW "<<(bw1 - ew2);
    std::cout<<" skipcut "<<skipcut;
    std::cout<<" dth "<<dth<<" therrs "<<eth2e<<" "<<bth1e<<" angcut "<<angcut<<std::endl;
  }
          if( bw1 < ew2 && (ew2 - bw1)  < skipcut && dth < angcut ) {
            // look for US and DS broken clusters w similar angle
            // US cluster 1 merge with DS cluster 2?
            float chgrat = 2 * fabs((bc1 - ec2) / (bc1 + ec2));
            // project sw1,st1 to ew2
            float dtim = fabs(bt1 + (ew2-bw1)*bs1 - et2);
  if(fDebugWire < 0) {
    std::cout<<" dtim "<<dtim<<" err "<<dtim<<" timecut "<<(int)timecut;
    std::cout<<" chgrat "<<chgrat<<" chgcut "<<chgcut;
    std::cout<<std::endl;
  }
            if(chgrat < chgcut && dtim < timecut) {
              cl2DoMerge(plnhits, tcl, it1, it2, 10);
              tclsize = tcl.size();
              break;
            }
          } // US cluster 1 merge with DS cluster 2
          
          if(bw2 < bw1 && ew2 > ew1) {
            // look for small cl2 within the wire boundary of cl1
            // with similar times and slopes for both clusters
            float dth = fabs(eth2 - eth1);
            float dtim = fabs(et1 +(ew2 - ew1 - 1)*es1 - et2);
            // count the number of wires with no hits on cluster 1
            short nmiss1 = bw1 - ew1 + 1 - tcl[it1].tclhits.size();
            // compare with the number of hits in cluster 2
            short nin2 = tcl[it2].tclhits.size();
  if(fDebugWire < 0) {
    std::cout<<"cl2: "<<tcl[it2].ID<<" within cl1 "<<tcl[it1].ID;
    std::cout<<" ? dth "<<dth<<" dtim "<<dtim<<" nmissed "<<nmiss1<<std::endl;
  }
            // make rough cuts before calling cl2Merge12
            // this may not work well for long wandering clusters
            bool didit = false;
            if(dth < 1 && dtim < timecut && nmiss1 >= nin2) 
                cl2ChkMerge12(plnhits, tcl, it1, it2, didit);
            if(didit) {
              tclsize = tcl.size();
              break;
            } //didit
          } // small cl2 within the wire boundary of cl1
          
          if(bw1 < bw2 && ew1 > ew2) {
            // look for small cl1 within the wire boundary of cl2
            // with similar times and slopes for both clusters
            float dth = fabs(eth2 - eth1);
            float dtim = fabs(et2 +(ew1 - ew2 - 1)*es2 - et1);
            // count the number of wires with no hits on cluster 2
            short nmiss2 = bw2 - ew2 + 1 - tcl[it2].tclhits.size();
            // compare with the number of hits in cluster 1
            short nin1 = tcl[it1].tclhits.size();
  if(fDebugWire < 0) {
    std::cout<<"cl1: "<<tcl[it1].ID<<" within cl2 "<<tcl[it2].ID;
    std::cout<<" ? dth "<<dth<<" dtim "<<dtim<<" nmissed "<<nmiss2<<std::endl;
  }
            // make rough cuts before calling cl2Merge12
            // this may not work well for long wandering clusters
            bool didit = false;
            if(dth < 1 && dtim < timecut && nmiss2 >= nin1) 
                cl2ChkMerge12(plnhits, tcl, it2, it1, didit);
            if(didit) {
              tclsize = tcl.size();
              break;
            } // didit
          } // small cl1 within the wire boundary of cl2
          
          if(tcl[it1].ID < 0) break;
        } // cluster 2
        if(tcl[it1].ID < 0) continue;
      } // cluster 1
      
      return;
    }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2ChkMerge12(art::PtrVector<recob::Hit>& plnhits, 
     std::vector<ClusterStore>& tcl, unsigned int it1, unsigned int it2, bool& didit)
  {
    // Calling routine has done a rough check that cluster it2 is a candidate
    // for merging with cluster it1. The wire range spanned by it2 lies
    // within the wire range of it1 and the clusters are reasonably close
    // together in time.
    // Flag cluster 2 if angle and time cuts are met but the charge cut is not.
    // This is an indicator that cluster 2 is a delta ray or e+e- pair separating
    
    // assume that no merging was done
    didit = false;

  if(fDebugWire < 0) std::cout<<"cl2ChkMerge12 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<std::endl;
    
    ClusterStore& cl1 = tcl[it1];
    std::map<short, short> HitNum;
    // fill a map of the hit time and charge on each wire for cluster 1
    for(std::vector<short>::iterator iht = cl1.tclhits.begin();
          iht != cl1.tclhits.end(); ++iht) {
      short hit = *iht;
      unsigned int wire = plnhits[hit]->WireID().Wire;
      HitNum[wire] = hit;
    }
    short ew2 = tcl[it2].EndWir;
    short et2 = tcl[it2].EndTim;
    // look for the closest wire with a hit on cluster 1
    short cwon1 = -1;
    // count the number of missing hits
    short nmiss = 0;
    for(short wire = ew2 - 1; wire > tcl[it1].EndWir; wire--) {
      if(HitNum[wire] > 0) {
        cwon1 = wire;
        break;
      }
      nmiss++;
    } // wire
  if(fDebugWire < 0) std::cout<<"chk next US wire "<<cwon1<<" missed "<<nmiss<<std::endl;
    if(cwon1 < 0) return;
    if(nmiss > fMaxWirSkip[pass]) return;
    
    // decode the pass for both clusters and select the matching cuts
    short pass1 = tcl[it1].ProcCode - 10 * (tcl[it1].ProcCode / 10);
    short pass2 = tcl[it2].ProcCode - 10 * (tcl[it2].ProcCode / 10);
    short cpass = pass1;
    // use the tighter cuts
    if(pass2 < pass1) cpass = pass2;
    
    // ******************* Check End of Cluster 2 matching with cluster 1
    // fit this section of cluster 1 with 4 hits starting at the hit on the
    // closest wire and moving US
    short hiton1 = HitNum[cwon1];
    cl2FitMid(plnhits, tcl, it1, hiton1, 4);
    // fit parameters are now in clpar. Charge is in fAveChg
    // check for angle consistency
    float dth = fabs(atan(ScaleF * clpar[1]) - atan(ScaleF * tcl[it2].EndSlp));
  if(fDebugWire < 0) std::cout<<"US dtheta "<<dth<<" cut "<<fKinkAngCut[cpass]<<std::endl;
    if(dth > fKinkAngCut[cpass]) return;
    // project the fit to the end wire of cluster 2
    float dtim = fabs(clpar[0] + (ew2 - wire0) * clpar[1] - et2);
  if(fDebugWire < 0) std::cout<<"US dtime "<<dtim<<std::endl;
    if(dtim > fTimeDelta[cpass]) return;  
    // make a charge ratio cut
    float chgrat = 2 * fabs(fAveChg - tcl[it2].EndChg) / (fAveChg + tcl[it2].EndChg);
  if(fDebugWire < 0)  std::cout<<"US chgrat "<<chgrat<<" cut "<<fChgCut[pass]<<std::endl;
    
    // ensure that there is a signal on any missing wires at the US end of 1
    short wir1 = plnhits[hiton1]->WireID().Wire;
    short tim1 = plnhits[hiton1]->PeakTime();
    bool SigOK = false;
    cl2ChkSignal(plnhits, wir1, tim1, ew2, et2, SigOK);
  if(fDebugWire < 0)  std::cout<<"US SigOK? "<<SigOK<<std::endl;
    if( !SigOK ) return;


    // ******************* Check Begin of Cluster 2 matching with cluster 1
    // next make a downstream charge cut. Start by finding the closest 
    // hit on cluster 1 near the beginning of cluster 2
    short bw2 = tcl[it2].BeginWir;
    short bt2 = tcl[it2].BeginTim;
    nmiss = 0;
    for(short wire = bw2 + 1; wire < tcl[it1].BeginWir; wire++) {
      if(HitNum[wire] > 0) {
        cwon1 = HitNum[wire];
        break;
      }
      nmiss++;
    }
    if(cwon1 < 0) return;
    if(nmiss > fMaxWirSkip[pass]) return;
    // fit this section of cluster 1 with 4 hits starting at the hit on the
    // closest wire and moving DS
    hiton1 = HitNum[cwon1];
    cl2FitMid(plnhits, tcl, it1, hiton1, -4);
    // check for angle consistency
    dth = fabs(atan(ScaleF * clpar[1]) - atan(ScaleF * tcl[it2].BeginSlp));
  if(fDebugWire < 0) std::cout<<"DS dtheta "<<dth<<" cut "<<fKinkAngCut[cpass]<<std::endl;
    if(dth > fKinkAngCut[cpass]) return;
    // project the fit to the begin wire of cluster 2
    dtim = fabs(clpar[0] + (bw2 - wire0) * clpar[1] - bt2);
  if(fDebugWire < 0) std::cout<<"DS dtime "<<dtim<<std::endl;
    if(dtim > fTimeDelta[cpass]) return;  
    // make a charge ratio cut
    chgrat = 2 * fabs(fAveChg - tcl[it2].BeginChg) / (fAveChg + tcl[it2].BeginChg);
  if(fDebugWire < 0)  std::cout<<"DS chgrat "<<chgrat<<" cut "<<fChgCut[pass]<<std::endl;
    
    // ensure that there is a signal on any missing wires at the US end of 1
    wir1 = plnhits[hiton1]->WireID().Wire;
    tim1 = plnhits[hiton1]->PeakTime();
    SigOK = false;
    cl2ChkSignal(plnhits, wir1, tim1, bw2, bt2, SigOK);
  if(fDebugWire < 0)  std::cout<<"DS SigOK? "<<SigOK<<std::endl;
    if( !SigOK ) return;

  if(fDebugWire < 0)  std::cout<<"Merge em"<<std::endl;
    // success. Merge them
    cl2DoMerge(plnhits, tcl, it1, it2, 100);
    didit = true;
    return;
  }


/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2DoMerge(art::PtrVector<recob::Hit>& plnhits, 
     std::vector<ClusterStore>& tcl, unsigned int it1, unsigned int it2,
     short ProcCode)
  {
    // Merge clusters. Cluster 1 has precedence for assignment of hits
    ClusterStore& cl1 = tcl[it1];
    ClusterStore& cl2 = tcl[it2];
    // ensure that there is only one hit/wire on both clusters
    std::map<short, short> wirehit;
    short hiwire = -32700;
    short lowire = 32700;
    for(std::vector<short>::iterator iht = cl1.tclhits.begin();
          iht != cl1.tclhits.end(); ++iht) {
      short hit = *iht;
      short wire = plnhits[hit]->WireID().Wire;
      if(wire < lowire) lowire = wire;
      if(wire > hiwire) hiwire = wire;
      wirehit[wire] = hit;
    }
    for(std::vector<short>::const_iterator iht = cl2.tclhits.begin();
          iht != cl2.tclhits.end(); ++iht) {
      short hit = *iht;
      short wire = plnhits[hit]->WireID().Wire;
      if(wire < lowire) lowire = wire;
      if(wire > hiwire) hiwire = wire;
      if(wirehit[wire] == 0) {
        wirehit[wire] = hit;
      } else {
        // a hit from cluster cl1 is on this wire. Free it up for later use
        short freehit = wirehit[wire];
        hiterr2[freehit] = fabs(hiterr2[freehit]);
        wirehit[wire] = hit;
      }
    }
    // make a new cluster
    fcl2hits.clear();
    for(short wire = hiwire; wire >= lowire; wire--) {
      if(wirehit[wire] > 0) fcl2hits.push_back(wirehit[wire]);
    }
    // re-fit the end of the cluster
    std::vector<short>::iterator iend = fcl2hits.end() - 1;
    short jend = *iend;
    wire0 = plnhits[jend]->WireID().Wire;
    cl2Fit(plnhits);
    clEndSlp = (float)clpar[1];
    clEndSlpErr = (float)clparerr[1];
    clEndChg = fAveChg;
    // re-fit the beginning of the cluster, using the first 4 hits
    std::vector<short>::iterator ibeg3 = fcl2hits.begin() + 3;
    short jbeg3 = *ibeg3;
    wire0 =  plnhits[jbeg3]->WireID().Wire;
    cl2Fit(plnhits);
    clBeginSlp = (float)clpar[1];
    clBeginSlpErr = (float)clparerr[1];
    clBeginChg = fAveChg;
    clStopCode = cl1.StopCode;
    clAssn = -1;
    // append it to the tcl vector
    cl2TmpStore(plnhits, tcl);
    unsigned int itnew = tcl.size()-1;
  if(fDebugWire < 0) std::cout<<"DoMerge "<<cl1.ID<<" "<<cl2.ID<<" -> "<<tcl[itnew].ID<<std::endl;
    // stuff the processor code with the higher pass cluster
    short pass1 = tcl[it1].ProcCode - 10 * (tcl[it1].ProcCode / 10);
    short pass2 = tcl[it2].ProcCode - 10 * (tcl[it2].ProcCode / 10);
    if(pass1 > pass2) {
      tcl[itnew].ProcCode = ProcCode + pass1;
    } else {
      tcl[itnew].ProcCode = ProcCode + pass2;
    }
    // mark cl1 and cl2 obsolete
    tcl[it1].ID = -tcl[it1].ID;
    tcl[it2].ID = -tcl[it2].ID;
    // move any associations to clusters 1 and 2 to the new cluster
    for(unsigned int ii = 0; ii < itnew; ii++) {
      if(tcl[ii].ID > 0 && tcl[ii].Assn >= 0) {
        if(tcl[ii].Assn == (short)it1 || tcl[ii].Assn == (short)it2) tcl[ii].Assn = itnew;
      }
    }
    return;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2Print(art::PtrVector<recob::Hit>& plnhits, 
     std::vector<ClusterStore>& tcl)
  {
    // prints clusters to the screen for code development
    std::cout<<"  ID nht Stop  Proc   beg_W:H   begT  bTheta Therr begChg  end_W:H  endT  eTheta Therr endChg";
    std::cout<<" bVx eVx";
    std::cout<<std::endl;
    for(unsigned int ii = 0; ii < tcl.size(); ii++) {
      std::vector<short>::const_iterator ihtb = tcl[ii].tclhits.begin();
      short hitb = *ihtb;
      std::vector<short>::const_iterator ihte = tcl[ii].tclhits.end()-1;
      short hite = *ihte;
      std::cout<<std::right<<std::setw(4)<<tcl[ii].ID;
      std::cout<<std::right<<std::setw(5)<<tcl[ii].tclhits.size();
      std::cout<<std::right<<std::setw(4)<<tcl[ii].StopCode;
      std::cout<<std::right<<std::setw(6)<<tcl[ii].ProcCode;
      std::cout<<std::right<<std::setw(6)<<tcl[ii].BeginWir<<":"<<hitb;
      if(hitb < 10) {
        std::cout<<"   ";
      } else if(hitb < 100) {
        std::cout<<"  ";
      } else if(hitb < 1000) std::cout<<" ";
      std::cout<<std::right<<std::setw(6)<<(short)tcl[ii].BeginTim;
//      std::cout<<std::right<<std::setw(7)<<std::setprecision(3)<<tcl[ii].BeginSlp;
      float arg = ScaleF * tcl[ii].BeginSlp;
      float theta = atan(arg);
      std::cout<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<theta;
      float thetaerr =  ScaleF * tcl[ii].BeginSlpErr / (1 + arg * arg);
      std::cout<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<thetaerr;
      std::cout<<std::right<<std::setw(5)<<(short)tcl[ii].BeginChg;
      std::cout<<std::right<<std::setw(6)<<tcl[ii].EndWir<<":"<<hite;
      if(hite < 10) {
        std::cout<<"   ";
      } else if(hite < 100) {
        std::cout<<"  ";
      } else if(hite < 1000) std::cout<<" ";
      std::cout<<std::right<<std::setw(6)<<(short)tcl[ii].EndTim;
//      std::cout<<std::right<<std::setw(7)<<std::setprecision(3)<<tcl[ii].EndSlp;
      arg = ScaleF * tcl[ii].EndSlp;
      theta = atan(arg);
      std::cout<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<theta;
      thetaerr =  ScaleF * tcl[ii].EndSlpErr / (1 + arg * arg);
      std::cout<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<thetaerr;
      theta = atan(tcl[ii].EndSlp);
      std::cout<<std::right<<std::setw(5)<<(short)tcl[ii].EndChg;
      std::cout<<std::right<<std::setw(5)<<tcl[ii].BeginVtx;
      std::cout<<std::right<<std::setw(5)<<tcl[ii].EndVtx;
      std::cout<<std::endl;
    }
    // print out lots of stuff
/*
    for(unsigned int ii = 0; ii < tcl.size(); ii++) {
      std::vector<short>::const_iterator ihts = tcl[ii].tclhits.begin();
      if(tcl[ii].ID > -1) {
        std::cout<<tcl[ii].ID<<" W:H ";
        for(std::vector<short>::const_iterator iht = tcl[ii].tclhits.begin();
              iht != tcl[ii].tclhits.end(); ++iht) {
          int hit = *iht;
          int wir =  plnhits[hit]->WireID().Wire;
          std::cout<<wir<<":"<<hit<<" ";
        }
      std::cout<<std::endl;
      }
    }
*/
  }


/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2TmpStore(art::PtrVector<recob::Hit>& plnhits, 
     std::vector<ClusterStore>& tcl)
  {

    if(fcl2hits.size() == 0) return;
    
    NClusters++;

    // flag all the hits as used
    for(std::vector<short>::const_iterator it = fcl2hits.begin(); it != fcl2hits.end(); ++it) {
      short hit = *it;
      hiterr2[hit] = -fabs(hiterr2[hit]);
    }

    // ensure that the cluster begin/end info is correct
    std::vector<short>::const_iterator ibg = fcl2hits.begin();
    short hitb = *ibg;
    std::vector<short>::const_iterator iend = fcl2hits.end() - 1;
    short hite = *iend;

    // store the cluster in the temporary ClusterStore struct
    ClusterStore clstr;
    
    clstr.ID = NClusters;
    clstr.BeginSlp    = clBeginSlp;
    clstr.BeginSlpErr = clBeginSlpErr;
    clstr.BeginWir    = plnhits[hitb]->WireID().Wire;
    clstr.BeginTim    = plnhits[hitb]->PeakTime();
    clstr.BeginChg    = plnhits[hitb]->Charge();
    clstr.EndSlp      = clEndSlp;
    clstr.EndSlpErr   = clEndSlpErr;
    clstr.EndWir      = plnhits[hite]->WireID().Wire;
    clstr.EndTim      = plnhits[hite]->PeakTime();
    clstr.EndChg      = plnhits[hite]->Charge();
    clstr.StopCode    = clStopCode;
    clstr.ProcCode    = clProcCode;
    clstr.Assn        = clAssn;
    clstr.BeginVtx    = -99;
    clstr.EndVtx      = -99;
    clstr.tclhits     = fcl2hits;
    tcl.push_back(clstr);
    return;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2FollowUS(art::PtrVector<recob::Hit>& plnhits)
  {
    // follow the cluster upstream

  std::vector<short>::const_iterator itt = fcl2hits.begin();
  int fhit = *itt;
  int fwir = plnhits[fhit]->WireID().Wire;
  prt = (fwir == fDebugWire && fhit == fDebugHit);

  if(prt) {
    std::cout<<"cl2FollowUS PASS "<<pass<<" Hits: ";
    for(std::vector<short>::const_iterator it = fcl2hits.begin(); it != fcl2hits.end(); ++it) {
      std::cout<<*it<<" ";
    }
    std::cout<<std::endl;
  }
    // SigOK = true if there is a ADC signal near the projected cluster position
    bool SigOK = true;
    bool HitOK = true;
    // count the number of missed hits on adjacent wires
    short nmissed = 0;
    // count the number of added hits after skipping
    short nHitAfterSkip = 0;
    bool DidaSkip = false;
    bool PostSkip = false;
    std::vector<short>::iterator it = fcl2hits.end() - 1;
    short lasthit = *it;
    short lastwire = plnhits[lasthit]->WireID().Wire;
  if(prt) std::cout<<"cl2FollowUS: last wire "<<lastwire<<" hit "<<lasthit<<std::endl;
    std::vector<float> chifits;
    
    for(short nextwire = lastwire-1; nextwire > fFirstWire-1; --nextwire) {
  if(prt) std::cout<<"cl2FollowUS: next wire "<<nextwire<<std::endl;
      // add hits and check for PH and width consistency
      cl2AddHit(plnhits, nextwire, true, HitOK, SigOK);
  if(prt) std::cout<<"cl2FollowUS: HitOK "<<HitOK<<" SigOK "<<SigOK<<std::endl;
      if(!HitOK) {
        // no hit on this wire. Was there a signal or dead wire?
        if(SigOK) {
          nmissed++;
          if(prt && nmissed > fMaxWirSkip[pass]) std::cout<<"nmissed break"<<std::endl;
          if(nmissed > fMaxWirSkip[pass]) {
            clStopCode = 1;
            break;
          }
          // see if we are in the PostSkip phase and missed more than 1 wire
//  std::cout<<"chk "<<nextwire<<" "<<PostSkip<<" "<<nmissed<<std::endl;
          if(PostSkip && nmissed > 1) {
//  std::cout<<"Missed wire after skip "<<nextwire<<" nmissed "<<nmissed<<std::endl;
            unsigned int newsize = fcl2hits.size() - nHitAfterSkip;
            fcl2hits.resize(newsize);
            std::vector<short>::reverse_iterator ii = fcl2hits.rbegin();
            short jj = *ii;
            wire0 = plnhits[jj]->WireID().Wire;
            cl2Fit(plnhits);
            clStopCode = 2;
            return;
          }
          if(nmissed > 1) {
            DidaSkip = true;
            PostSkip = false;
          }
        } else {
          clStopCode = 0;
          if(prt) std::cout<<"No hit or signal on wire "<<nextwire<<std::endl;
//  std::cout<<"No hit or signal on wire "<<nextwire<<std::endl;
          break;
        }
      } else {
        // update the fit
        // find the origin of the fit
        std::vector<short>::reverse_iterator ii = fcl2hits.rbegin();
        short jj = *ii;
  if(prt) std::cout<<"Add hit W:H "<<nextwire<<":"<<jj<<std::endl;
        wire0 = plnhits[jj]->WireID().Wire;
        cl2Fit(plnhits);
        // update the cluster Begin slope?
        if((short)fcl2hits.size() == fMaxHitsFit[pass] ||
           (short)fcl2hits.size() == fMinHits[pass]) {
          clBeginSlp = (float)clpar[1];
          clBeginSlpErr = (float)clparerr[1];
        }
        // monitor the onset of a kink. Find the average chisq for the fit
        // using the previous 3 - 6 hits. Look for a progressive increase
        // in chisq for the previous 0 - 2 hits.
        chifits.push_back(clChisq);
        if(chifits.size() > 7 && fKinkChiRat[pass] > 0) {
          short chsiz = chifits.size()-1;
  if(prt) {
    std::cout<<"Kink chk "<<chifits[chsiz]<<" "<<chifits[chsiz-1]<<" ";
    std::cout<<chifits[chsiz-2]<<" "<<chifits[chsiz-3]<<" ";
    std::cout<<chifits[chsiz-4]<<" "<<chifits[chsiz-5]<<" ";
    std::cout<<chifits[chsiz-6]<<" "<<chifits[chsiz-7]<<std::endl;
  }
          if( chifits[chsiz-3] > fKinkChiRat[pass] * chifits[chsiz-4] &&
              chifits[chsiz-2] > fKinkChiRat[pass] * chifits[chsiz-3] && 
              chifits[chsiz-1] > fKinkChiRat[pass] * chifits[chsiz-2] &&
              chifits[chsiz]   > fKinkChiRat[pass] * chifits[chsiz-1]) {
            // find the kink angle (crudely) from the 0th and 3rd hit
            std::vector<short>::reverse_iterator i0 = fcl2hits.rbegin();
            short ih0 = *i0;
            std::vector<short>::reverse_iterator i3 = fcl2hits.rbegin() + 3;
            short ih3 = *i3;
            float dt03 = plnhits[ih3]->PeakTime() - plnhits[ih0]->PeakTime();
            float dw03 = plnhits[ih3]->WireID().Wire - plnhits[ih0]->WireID().Wire;
            float th03 = atan( ScaleF * dt03 / dw03);
            // and the 4th and 7th hit
            std::vector<short>::reverse_iterator i4 = fcl2hits.rbegin() + 4;
            short ih4 = *i4;
            std::vector<short>::reverse_iterator i7 = fcl2hits.rbegin() + 7;
            short ih7 = *i7;
            float dt47 = plnhits[ih7]->PeakTime() - plnhits[ih4]->PeakTime();
            float dw47 = plnhits[ih7]->WireID().Wire - plnhits[ih4]->WireID().Wire;
            float th47 = atan(ScaleF * dt47 / dw47);
            float dth = fabs(th03 - th47);
  if(prt) std::cout<<" Kink angle "<<std::setprecision(3)<<dth<<std::endl;
            // cut on the allowed kink angle
            if(dth > fKinkAngCut[pass]) {
  if(prt) std::cout<<"stopped tracking "<<std::endl;
              // kill the last 4 hits, refit and return
              unsigned int newsize = fcl2hits.size() - 4;
              fcl2hits.resize(newsize);
              std::vector<short>::reverse_iterator ii = fcl2hits.rbegin();
              short jj = *ii;
              wire0 = plnhits[jj]->WireID().Wire;
              cl2Fit(plnhits);
              clStopCode = 3;
              break;
            } // kinkang check
          } // chifits check
        } // chifits.size() > 5
        // reset nmissed
        nmissed = 0;
        // start counting hits added after skipping
        if(DidaSkip) {
          // start PostSkip phase with one 
          PostSkip = true;
          DidaSkip = false;
          nHitAfterSkip = 0;
        } // DidaSkip
        // check for PostSkip phase
        if(PostSkip) {
          // end the PostSkip phase if there are enough hits
          nHitAfterSkip++;
          if(nHitAfterSkip == fMinWirAfterSkip[pass]) PostSkip = false;
        } 
        // check for the onset of a kink
        if(clChisq > fChiCut[pass]) {
          // remove the last hit and re-fit
          fcl2hits.erase(fcl2hits.end() - 1);
          // find the origin of the fit
          std::vector<short>::reverse_iterator ii = fcl2hits.rbegin();
          short jj = *ii;
          wire0 = plnhits[jj]->WireID().Wire;
          cl2Fit(plnhits);
          clStopCode = 4;
          break;
        } // clChisq cut
      } // !HitOK check
    } // nextwire
    
    // stopped normally
    if(clStopCode < 0) clStopCode = 0;

    // count the number of consecutive hits on wires at the end
    it = fcl2hits.end() - 1;
    lasthit = *it;
    lastwire = plnhits[lasthit]->WireID().Wire;
//  std::cout<<"Last wire "<<lastwire<<std::endl;
    short nht = 0;
    // number of hits to lop off the end
    short nlop = 0;
    for(std::vector<short>::reverse_iterator ii = fcl2hits.rbegin()+1;
            ii != fcl2hits.rend(); ++ii) {
      short jj = *ii;
      short wir = plnhits[jj]->WireID().Wire;
//  std::cout<<"chk "<<jj<<" wir "<<wir<<" lastwire "<<lastwire<<std::endl;
      if(abs(wir - lastwire) > 1) {
        nlop = nht + 1;
        break;
      }
      nht++;
      if(nht > fMinWirAfterSkip[pass]) {
        break;
      }
      lastwire = wir;
    }
    if(nlop > 0) {
//  std::cout<<"nlop "<<nlop<<std::endl;
      unsigned int newsize = fcl2hits.size() - nlop;
      fcl2hits.resize(newsize);
      clStopCode = 2;
    }
    prt = false;
    return;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2FitMid(art::PtrVector<recob::Hit>& plnhits,
      std::vector<ClusterStore>& tcl, unsigned int it1, short ihtin, short nhit)
  {
    // Fits hits on temp cluster it1 to a line starting at hit ihtin and including
    // nhit hits incrementing towards the hit vector End when nhit > 0 and
    // decrementing towards the hit vector Begin when nhit < 0.
    // The fit params are stashed in the clpar and clparerr arrays. 
    // fAveChg is re-calculated as well.
    
    
    // set chisq bad in case something doesn't work out
    clChisq = 99.;
    
    ClusterStore cls = tcl[it1];
    if(cls.tclhits.size() < 3) return;
    
    Int_t nht = 0;
    if(nhit > 0) {
      nht = nhit;
      // find the first desired hit and move towards the End
      fAveChg = 0.;
      Int_t hitcnt = 0;
      bool UseEm = false;
      for(std::vector<short>::iterator it = cls.tclhits.begin(); 
          it != cls.tclhits.end(); ++it) {
        short ihit = *it;
        // look for the desired first hit. Use this as the origin wire
        if(ihit == ihtin) {
          UseEm = true;
          wire0 = plnhits[ihit]->WireID().Wire;
        }
        // use hits after finding the first desired hit
        if(UseEm) {
          short wire = plnhits[ihit]->WireID().Wire;
          xwir[hitcnt] = (Double_t)(wire - wire0);
          ytim[hitcnt] = plnhits[ihit]->PeakTime();
          ytimerr[hitcnt] = sqrt(fabs(hiterr2[ihit]));
          fAveChg += plnhits[ihit]->Charge();
          if(hitcnt == nht) break;
          hitcnt++;
        }
      }    
      nht = hitcnt;
    } else {
      nht = -nhit;
      // find the first desired hit and move towards the Begin
      fAveChg = 0.;
      Int_t hitcnt = 0;
      bool UseEm = false;
      for(std::vector<short>::reverse_iterator it = cls.tclhits.rbegin(); 
          it != cls.tclhits.rend(); ++it) {
        short ihit = *it;
        // look for the desired first hit. Use this as the origin wire
        if(ihit == ihtin) {
          UseEm = true;
          wire0 = plnhits[ihit]->WireID().Wire;
        }
        // use hits after finding the first desired hit
        if(UseEm) {
          short wire = plnhits[ihit]->WireID().Wire;
          xwir[hitcnt] = (Double_t)(wire - wire0);
          ytim[hitcnt] = plnhits[ihit]->PeakTime();
          ytimerr[hitcnt] = sqrt(fabs(hiterr2[ihit]));
          fAveChg += plnhits[ihit]->Charge();
          if(hitcnt == nht) break;
          hitcnt++;
        }
      }    
      nht = hitcnt;
    }
    
    if(nht < 2) return;

    // use the ROOT linear fitter
    TLinearFitter *lf =  new TLinearFitter(2,"pol1");
    lf->AssignData(nht,1,xwir,ytim,ytimerr);
    // fit the points
    lf->Eval();
    lf->GetParameters(clpar);
    lf->GetErrors(clparerr);
    clChisq = lf->GetChisquare() / (nht - 2);
    fAveChg = fAveChg / float(nht);
    delete lf;
    
    return;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2Fit(art::PtrVector<recob::Hit>& plnhits)
  {
    // Fits the hits on a cluster with origin at wire0 defined by the calling
    // routine. This routine assumes that
    // wires are numbered from lower (upstream) to higher (downstream) and
    // that the hits in the fclhits vector are sorted so that upstream hits
    // are at the end of the vector
    
    Int_t nht;
    // fit all hits or truncate?
    if((short)fcl2hits.size() < fMaxHitsFit[pass]) {
      nht = fcl2hits.size();
    } else {
      nht = fMaxHitsFit[pass];
    }
    
    // apply an angle dependent scale factor. The error should be
    // wire pitch / sqrt(12) for a cluster at 90 degrees. This formula 
    // simply doubles the error, which I think is reasonable for uBooNE and
    // ArgoNeuT
    float angfactor = 1;
    if(fcl2hits.size() > 2) angfactor = 2 - 1/(1 + fabs(clpar[1]));

    // load the hits starting at the back end of the fcl2hits vector.
    // These are the most upstream hits
    short iht = 0;
  if(prt) std::cout<<"cl2Fit wire0 "<<wire0<<" W:H ";
    for(std::vector<short>::reverse_iterator it = fcl2hits.rbegin(); it != fcl2hits.rend(); ++it) {
      short ihit = *it;
      short wire = plnhits[ihit]->WireID().Wire;
  if(prt) std::cout<<wire<<":"<<ihit<<" ";
  if(iht > 20) std::cout<<"OOOPS "<<iht<<std::endl;
        xwir[iht] = (Double_t)(wire - wire0);
        ytim[iht] = plnhits[ihit]->PeakTime();
        ytimerr[iht] = angfactor * sqrt(fabs(hiterr2[ihit]));
        ychg[iht] = plnhits[ihit]->Charge();
        // set charge error to 30%
        ychgerr[iht] = 0.3 * ychg[iht];
        if(iht == nht) break;
        iht++;
    }
    
    nht = iht;
  if(prt) std::cout<<std::endl;
    
    if(nht < 2) {
      clChisq = 999.;
      return;
    }

    // use the ROOT linear fitter
    TLinearFitter *lf =  new TLinearFitter(2,"pol1");
    lf->AssignData(nht,1,xwir,ytim,ytimerr);
    // fit the points
    lf->Eval();
    lf->GetParameters(clpar);
    lf->GetErrors(clparerr);
    float chidof = 0;
    if(nht > 2) chidof = lf->GetChisquare() / (nht - 2);
    clChisq = chidof;

  if(prt) std::cout<<"nht "<<nht<<" fit par "<<(int)clpar[0]<<" "<<clpar[1]<<" chidof "<<chidof;
    
    std::vector<short>::reverse_iterator it0 = fcl2hits.rbegin();
    short ih0 = *it0;
    fAveWid = hitwid[ih0];
    // update the average charge
    if(fNHitsAve[pass] < 2) {
      // simply use the charge and width the last hit
      fAveChg = plnhits[ih0]->Charge();
    } else if(fNHitsAve[pass] == 2) {
      // average the last two points if requested
      std::vector<short>::reverse_iterator it1 = fcl2hits.rbegin()+1;
      short ih1 = *it1;
      fAveChg = (plnhits[ih0]->Charge() + plnhits[ih1]->Charge()) / 2.;
      fAveWid = (hitwid[ih0] + hitwid[ih1]) / 2.;
    } else {
      // fAveChg = -1 for a new cluster. Don't mess with this until
      // there are enough hits to fit. The charge ratio cut in cl2Fit is not
      // in effect until fAveChg > 0.
      if(nht > fNHitsAve[pass]) {
        // fit the charge to a line using the number of hits requested
        // but do not include the Begin hit. It could have a low charge
        // that would skew the fit, if this is a stopping particle
        if((unsigned int)nht == fcl2hits.size()) nht--;
        // don't calculate an average until there are sufficient points
        if(nht > 2) {
          lf->ClearPoints();
          lf->AssignData(nht,1,xwir,ychg,ychgerr);
          lf->Eval();
          lf->GetParameters(clchgpar);
          // clchgpar[0] is the estimated charge at wire0
          fAveChg = clchgpar[0];
  if(prt) std::cout<<" FitChg "<<(int)fAveChg<<" "<<clchgpar[1]<<std::endl;
        }
      }
    }

    delete lf;

    return;
  }


/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2AddHit(art::PtrVector<recob::Hit>& plnhits,
        short kwire, bool hitchk, bool& HitOK, bool& SigOK)
  {
    // Add a hit to the cluster if it meets several criteria:
    // similar pulse height to the cluster (if hitchk true)
    // similar hit width to the cluster (if hitchk true)
    // closest hit to the project cluster position.
    // Return SigOK if there is a nearby hit that was missed due to the cuts
    
    short firsthit = WireHitRange[kwire - fFirstWire].first;
    // skip bad wire, but assume the track was there
    if(firsthit == -1) {
      SigOK = true;
      HitOK = false;
      return;
    }
    // return if no signal and no hit
    if(firsthit == -2) {
      SigOK = false;
      HitOK = false;
      return;
    }
    short lasthit = WireHitRange[kwire - fFirstWire].second;
    // Determine if the last hit added was a large (low) charge hit
    // This will be used to prevent adding large (low) charge hits on two
    // consecutive fits
    std::vector<short>::reverse_iterator it1 = fcl2hits.rbegin();
    short ih1 = *it1;
    float bigchgcut = 2 * fChgCut[pass];
    bool lasthitbig = ( (plnhits[ih1]->Charge() / fAveChg) > bigchgcut);
    float lowchgcut = -1.5 * fChgCut[pass];
    bool lasthitlow = ( (plnhits[ih1]->Charge() / fAveChg) < fChgCut[pass]);
    // find the expected time of a hit on this wire
    float prtime = clpar[0] + (kwire - wire0) * clpar[1];
    // max number of time ticks between projected cluster and hit position
    float best = 50.;
    short imbest = -1;
    SigOK = false;
    for(short khit = firsthit; khit < lasthit; khit++) {
  if(prt) {
    std::cout<<"cl2AddHit chk W:H "<<kwire<<":"<<khit<<" time "<<(int)plnhits[khit]->PeakTime();
    std::cout<<" prtime "<<(short)prtime<<" fAveChg "<<(int)fAveChg;
    std::cout<<" lasthitlow "<<lasthitlow<<" lasthitbig "<<lasthitbig<<std::endl;
  }
      if(hiterr2[khit] < 0) continue;
      // make hit charge and width cuts
      if(hitchk) {
        float chgrat = (plnhits[khit]->Charge() - fAveChg) / fAveChg;
  if(prt) std::cout<<" Chg "<<(short)plnhits[khit]->Charge()<<" chgrat "<<chgrat<<std::endl;
        // Fudge for poorly reconstructed hits with large PH
        if(chgrat > 0.5) {
          if(prtime < plnhits[khit]->PeakTime() + fFudgeBigHits * hitwid[khit] && 
             prtime > plnhits[khit]->PeakTime() - fFudgeBigHits * hitwid[khit]) SigOK = true;
        } else {
          if(prtime < plnhits[khit]->PeakTime() + hitwid[khit] && 
             prtime > plnhits[khit]->PeakTime() - hitwid[khit]) SigOK = true;
        }
        // don't make a charge ratio cut until fAveChg is defined
        if(fAveChg > 0.) {
          if(lasthitlow) {
            // last hit added was low so don't allow another one.
            // We should really be lopping off the previous hit here...
            if(chgrat < -fChgCut[pass]) continue;
          } else {
            // allow a low hit
            if(chgrat < lowchgcut) continue;
          }
          // make the high charge ratio cut. 
          // decide on a high charge cut
          if(lasthitbig) {
            // last hit added was big so don't allow another one.
            // We should really be lopping off the previous hit here...
          if(chgrat > fChgCut[pass]) continue;
          } else {
            // allow a big hit
            if(chgrat > bigchgcut) continue;
          }
        } // fAveChg > 0
/* It isn't clear that this is a worthwhile cut to make. Needs more study
        float widrat = fabs(hitwid[khit] - fAveWid) / fAveWid;
  if(prt) std::cout<<" widrat "<<widrat<<std::endl;
        if(widrat > fWidCut[pass]) continue;
*/
      } // hitchk
      float timediff = (plnhits[khit]->PeakTime() - prtime);
  if(prt) std::cout<<" time diff "<<timediff<<std::endl;
      if(fabs(timediff) < best) {
        best = fabs(timediff);
        imbest = khit;
      } // timediff
    } // khit loop
    if(imbest < 0) {
      HitOK = false;
      return;
    }
    // Found a close hit check the chisq
    float prtimerr2 = clparerr[0]*clparerr[0] + fabs(kwire-wire0)*clparerr[1]*clparerr[1];
  if(prt) {
    std::cout<<"clerr "<<prtimerr2<<" hiterr "<<hiterr2[imbest];
    std::cout<<" best "<<best<<std::endl;
  }
    float err = sqrt(prtimerr2 + fabs(hiterr2[imbest]));
    // increase the error for large angle clusters
    if(abs(clpar[1]) > 30) err *= 3.;
    // (number of sigma)^2 difference
    float numsig2 = best / err;
    if(numsig2 < 10.) {
      fcl2hits.push_back(imbest);
      if(prt) std::cout<<"cl2AddHit ADD W:H "<<kwire<<":"<<imbest<<" best "<<best<<std::endl;
      HitOK = true;
      return;
    } else {
      if(prt) std::cout<<"cl2AddHit bad chisq "<<numsig2<<std::endl;
      HitOK = false;
      return;
    }
  }


} // namespace cluster
