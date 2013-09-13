//////////////////////////////////////////////////////////////////////
///
/// ClusterCrawlerAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithm for crawling along a string of hits to make line clusters
/// Technical note in MicroBooNE docdb #2831
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
#include "cetlib/exception.h"

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

namespace cluster {
//------------------------------------------------------------------------------
  ClusterCrawlerAlg::ClusterCrawlerAlg(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
    cl2Init();
  }

//------------------------------------------------------------------------------
  ClusterCrawlerAlg::~ClusterCrawlerAlg()
  {
  }

  void ClusterCrawlerAlg::reconfigure(fhicl::ParameterSet const& pset)
  { 
    fNumPass            = pset.get<             unsigned short  >("NumPass");
    fMaxHitsFit         = pset.get< std::vector<unsigned short> >("MaxHitsFit");
    fMinHits            = pset.get< std::vector<unsigned short> >("MinHits");
    fNHitsAve           = pset.get< std::vector<unsigned short> >("NHitsAve");
    fChgCut             = pset.get< std::vector<float> >("ChgCut");
    fWidCut             = pset.get< std::vector<float> >("WidCut");
    fChiCut             = pset.get< std::vector<float> >("ChiCut");
    fMaxWirSkip         = pset.get< std::vector<unsigned short> >("MaxWirSkip");
    fMinWirAfterSkip    = pset.get< std::vector<unsigned short> >("MinWirAfterSkip");
    fKinkChiRat         = pset.get< std::vector<float> >("KinkChiRat");
    fKinkAngCut         = pset.get< std::vector<float> >("KinkAngCut");
    fDoMerge            = pset.get< std::vector<bool>  >("DoMerge");
    fTimeDelta          = pset.get< std::vector<float> >("TimeDelta");

    fHitErrFac          = pset.get<             float  >("HitErrFac");
    fHitWidFac          = pset.get<             float  >("HitWidFac");
    fBEChgRat           = pset.get<             float  >("BEChgRat");
    fPairAngCut         = pset.get<             float  >("PairAngCut");
    fCurlyMergeAngCut   = pset.get<             float  >("CurlyMergeAngCut");
    fDoVertex           = pset.get<             bool   >("DoVertex");
    fDebugWire          = pset.get<             short  >("DebugWire");
    fDebugHit           = pset.get<             short  >("DebugHit");

    // some error checking
    bool badinput = false;
    if(fNumPass > fMaxHitsFit.size()) badinput = true;
    if(fNumPass > fMinHits.size()) badinput = true;
    if(fNumPass > fNHitsAve.size()) badinput = true;
    if(fNumPass > fChgCut.size()) badinput = true;
    if(fNumPass > fWidCut.size()) badinput = true;
    if(fNumPass > fChiCut.size()) badinput = true;
    if(fNumPass > fMaxWirSkip.size()) badinput = true;
    if(fNumPass > fMinWirAfterSkip.size()) badinput = true;
    if(fNumPass > fKinkChiRat.size()) badinput = true;
    if(fNumPass > fKinkAngCut.size()) badinput = true;
    if(fNumPass > fDoMerge.size()) badinput = true;
    if(fNumPass > fTimeDelta.size()) badinput = true;

    if(badinput) throw cet::exception("ClusterCrawler")<<"Bad input from fcl file ";

  }

  // used for sorting hits on wires
  bool SortByLowHit(short i, short j) {return ((i > j));}
  // used for sorting clusters by length
//  typedef std::pair<unsigned int, unsigned int> mypair;
//  bool SortByLen(const mypair& L, const mypair& R) {return (L.first > R.first);}

  void ClusterCrawlerAlg::cl2Init() {
    prt = false;
    NClusters = 0;  clBeginSlp = 0; clBeginSlpErr = 0; clBeginTim = 0;
    clBeginWir = 0; clBeginChg = 0; clEndSlp = 0;      clEndSlpErr = 0;
    clEndTim = 0;   clEndWir = 0;   clEndChg = 0;      clChisq = 0;
    clStopCode = 0; clProcCode = 0; clAssn = 0;        fFirstWire = 0;
    fLastWire = 0; fAveWid = 0; fAveChg = 0; pass = 0; ScaleF = 0;
    tcl.clear(); vtx.clear(); WireHitRange.clear(); hiterr2.clear();
    hitwid.clear();
  }


  void ClusterCrawlerAlg::RunCrawler(const art::PtrVector<recob::Hit>& plnhits,int plane)
  {
    // The hit collection is assumed to be sorted in increasing wire order

    cl2Init();

    if(plnhits.size() <= 2) return;

    // get the scale factor to convert dTick/dWire to dX/dU. This is used
    // to make the kink and merging cuts
    uint32_t channel = plnhits[0]->Wire()->RawDigit()->Channel();
    float wirePitch = geom->WirePitch(geom->View(channel));
    float tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
    tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
    ScaleF = tickToDist / wirePitch;

    art::PtrVector<recob::Hit>::const_iterator it = plnhits.begin();
    fFirstWire = (*it)->WireID().Wire;
    it = plnhits.end()-1;
    fLastWire = (*it)->WireID().Wire;

    // return if not enough hits on wires
    if((fLastWire - fFirstWire) < fMinHits[fNumPass - 1]) return;
    
    prt = (fDebugWire > 0);
    
    // try this out
//    WireHitRange.reserve(fLastWire - fFirstWire + 1);
    // initialize WireHitRange to "no hits on wire" condition
    short sflag = -2;
    for(unsigned short wire = fFirstWire; wire <= fLastWire; ++wire) {
      WireHitRange.push_back(std::make_pair(sflag, sflag));
    }

    // find dead wires in this region
    filter::ChannelFilter cf;
    sflag = -1;
    for(unsigned short wire = fFirstWire+1; wire < fLastWire; ++wire) {
      unsigned int pln = plane;
      unsigned int wir = wire;
      uint32_t chan = geom->PlaneWireToChannel(pln,wir);
      // remember to offset references to WireHitRange by the FirstWire
      unsigned short index = wire - fFirstWire;
      if(cf.BadChannel(chan)) WireHitRange[index] = std::make_pair(sflag, sflag);
    }
    
    unsigned short lastwire = fFirstWire;
    unsigned short thishit = 0;
    unsigned short lastfirsthit = 0;
    unsigned short index = 0;
  if(prt) {
    std::cout<<"************* hits in plane "<<plane<<" ******************"<<std::endl;
    std::cout<<" W:H Time Wid Err Chg  Chi Mult"<<std::endl;
  }
    for(auto hitIter = plnhits.begin(); hitIter != plnhits.end(); ++hitIter) {
      unsigned short thiswire = (*hitIter)->WireID().Wire;
      float hitep = (*hitIter)->EndTime()-(*hitIter)->PeakTime();
      // calculate the hit time uncertainty
      float arg = fHitErrFac * hitep;
      hiterr2.push_back(arg * arg);
      // and width for checking that the Signal is OK on a wire
      hitwid.push_back(fHitWidFac * hitep);
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
        index = lastwire - fFirstWire;
        short itmp1 = lastfirsthit;
        short itmp2 = thishit;
        WireHitRange[index] = std::make_pair(itmp1,itmp2);
        lastwire = thiswire;
        lastfirsthit = thishit;
      } else if(thiswire < lastwire) {
        mf::LogError("ClusterCrawler")<<"ERROR: Hits not sorted!!";
        return;
      }
      ++thishit;
    }
    index = fLastWire - fFirstWire;
    short itmp1 = lastfirsthit;
    short itmp2 = thishit;
    WireHitRange[index] = std::make_pair(itmp1,itmp2);
    
    prt = false;
    unsigned short nHitsUsed = 0;
    bool AllDone = false;
    for(unsigned short thispass = 0; thispass < fNumPass; ++thispass) {
      pass = thispass;
  if(fDebugWire != 0) {
    std::cout<<"******** ClusterCrawler plane "<<plane;
    std::cout<<" pass "<<pass<<" ****************"<<std::endl;
  }
      // look for a starting cluster that spans a block of wires
      unsigned short span = 3;
      if(fMinHits[pass] < span) span = fMinHits[pass];
      for(unsigned short iwire = fLastWire; iwire > fFirstWire; iwire--) {
        unsigned short index = iwire - fFirstWire;
        // skip bad wires or no hits on the wire
        if(WireHitRange[index].first < 0) continue;
        unsigned short ifirsthit = WireHitRange[index].first;
        unsigned short ilasthit = WireHitRange[index].second;
        for(unsigned short ihit = ifirsthit; ihit < ilasthit; ++ihit) {
          bool ClusterAdded = false;
          // skip used hits
          if(ihit > plnhits.size()-1) {
            mf::LogError("ClusterCrawler")<<"RunCrawler bad ihit "<<ihit;
            return;
          }
          if(hiterr2[ihit] < 0) continue;
          // skip multiple hits except on the last pass
          if(pass < fNumPass - 1 && plnhits[ihit]->Multiplicity() > 1) continue;
          if((iwire - span + 1) < fFirstWire) continue;
          for(unsigned short jwire = iwire - span + 1; jwire < iwire; ++jwire) {
            unsigned short jindx = jwire - fFirstWire;
            if(WireHitRange[jindx].first < 0) continue;
            unsigned short jfirsthit = WireHitRange[jindx].first;
            unsigned short jlasthit = WireHitRange[jindx].second;
            for(unsigned short jhit = jfirsthit; jhit < jlasthit; ++jhit) {
              if(jhit > plnhits.size()-1) {
                mf::LogError("ClusterCrawler")<<"RunCrawler bad jhit "<<jhit;
                return;
              }
              if(hiterr2[jhit] < 0) continue;
              // skip multiple hits
              if(pass < fNumPass - 1 && plnhits[jhit]->Multiplicity() > 1) continue;
              // start a cluster with these two hits
              fcl2hits.clear();
              fAveWid = -1.;
              fAveChg = -1.;
              clBeginChg = -1.;
              clStopCode = 0;
              clProcCode = pass;
              clAssn = -1; 
              fcl2hits.push_back(ihit);
              fcl2hits.push_back(jhit);
              clpar[0] = plnhits[jhit]->PeakTime();
              clpar[1] = (plnhits[ihit]->PeakTime() - plnhits[jhit]->PeakTime()) / (iwire - jwire);
              clChisq = 0;
              // now look for hits to add on the intervening wires
              bool SigOK = false;
              bool HitOK = false;
              bool clok = true;
              for(unsigned short kwire = jwire+1; kwire < iwire; ++kwire) {
                cl2AddHit(plnhits, kwire, HitOK, SigOK);
                // no hit added and no nearby hit either
                if(!HitOK && !SigOK) {
                  clok = false;
                  break;
                }
              }
              // drop it?
              if(fcl2hits.size() < span || !clok) continue;
              // sort them by decreasing wire number
              // assume that this is the same as sorting by decreasing 
              // hit number. This only needs to be done on the starting cluster
              // hits will be added in the proper order by cl2Follow
              std::sort(fcl2hits.begin(), fcl2hits.end(), SortByLowHit);
              // re-fit
              // define the hit origin
              cl2Fit(plnhits);
              if(clChisq > 10.) continue;
              // check the charge ratio between the DS hit and the next-most
              // DS hit. This ratio is < 2 for a stopping particle. A large
              // ratio indicates that we are starting a cluster just US of a
              // high ionization region
              float chg0 = plnhits[fcl2hits[0]]->Charge();
              float chg1 = plnhits[fcl2hits[1]]->Charge();
              if(chg0 > 2 * chg1) continue;
              // save the cluster begin info
              clBeginWir = iwire;
              clBeginTim = plnhits[ihit]->PeakTime();
              clBeginSlp = clpar[1];
              clBeginSlpErr = clparerr[1];
              // follow a trail of hits upstream
              cl2FollowUS(plnhits);
              if(fcl2hits.size() >= fMinHits[pass]) {
                // it's long enough so save it
                clEndSlp = clpar[1]; // save the slope at the end
                clEndSlpErr = clparerr[1];
                clEndChg = fAveChg;
                cl2TmpStore(plnhits, tcl); // store the cluster
                ClusterAdded = true;
                nHitsUsed += fcl2hits.size();
                AllDone = (nHitsUsed == plnhits.size());
                break;
              }
              if(pass < fNumPass - 1) {
                // Is it long enough for the next pass?
                if(fcl2hits.size() >= fMinHits[pass+1]) {
                  clEndSlp = clpar[1]; // save the slope at the end
                  clEndSlpErr = clparerr[1];
                  clEndChg = fAveChg;
                  // set a special code 
                  clProcCode += 2000;
                  cl2TmpStore(plnhits, tcl);
                  ClusterAdded = true;
                  nHitsUsed += fcl2hits.size();
                  AllDone = (nHitsUsed == plnhits.size());
                  break;
                } // long enough
              } // pass < fNumPass
              // kill it
            } // jhit
            if(ClusterAdded || AllDone) break;
          } // jwire
          if(ClusterAdded || AllDone) break;
        } // ihit
        if(AllDone) break;
      } // iwire
      // try to merge clusters 
      if(fDoMerge[pass]) cl2ChkMerge(plnhits, tcl);
      if(fDoMerge[pass] && pass > 0 && fDebugWire < 0) {
        std::cout<<"After merging in plane = "<<plane<<std::endl;
        cl2Print(plnhits, tcl);
      }
      if(AllDone) break;
    } // pass


  if(fDebugWire != 0 ) {
    std::cout<<"Clustering done in plane = "<<plane<<std::endl;
    cl2Print(plnhits, tcl);
  }


    
    // prepare close pair clusters for 3D matching
//    if(fPairAngCut > 0.) cl2ChkPair(plnhits, tcl);
    
    // try to merge short curly clusters
//    if(fCurlyMergeAngCut > 0.) cl2CurlyMerge(plnhits, tcl);

    short ncl = 0;
    short nht = 0;
    for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
      if(tcl[ii].ID > 0) {
        ++ncl;
        nht += tcl[ii].tclhits.size();
      }
    }
//    std::cout<<"Hits used in plane "<<plane<<" "<<nht<<std::endl;

    if(fDoVertex) {
      // split clusters using vertexing information
//      cl2VtxClusterSplit(plnhits, tcl);
      // find vertices with final set of clusters
      cl2DoVertex(plnhits, tcl, vtx);
    }
    
    // re-define the beginning and end of the cluster using the average charge
    // ratio if it is significant
    if(fBEChgRat > 0.) cl2SetBeginEnd(plnhits, tcl);
    
//    std::cout<<"Clustering done in plane "<<plane<<" nclusters "<<tcl.size()<<std::endl;
    
    WireHitRange.clear();
    hiterr2.clear();
    hitwid.clear();    
//  std::cout<<"Alg out"<<std::endl;
    
  } // RunCrawler


/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2VtxClusterSplit(const art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl)
    {
      // do some preliminary vertexing to look for clusters that span a
      // vertex and should therefore have been split. This routine should
      // only be called after crawling is complete, since it messes with
      // the pass variable
      
      if(tcl.size() < 2) return;
      
      float maxtime = detprop->NumberTimeSamples();
      
      vtx.clear();

  if(fDebugHit < 0) std::cout<<"cl2VtxClusterSplit check "<<std::endl;
      
      for (unsigned short it1 = 0; it1 < tcl.size(); ++it1) {
        // ignore abandoned clusters
        if(tcl[it1].ID < 0) continue;
        float es1 = tcl[it1].EndSlp;
        unsigned short ew1 = tcl[it1].EndWir;
        // don't bother with clusters near the End of the hit collection
        if(ew1 < fFirstWire + 3) continue;
        float et1 = tcl[it1].EndTim;
        // project to the next US wire and look for a nearby hit
        unsigned short uswire = ew1 - 1;
        unsigned short index = uswire - fFirstWire;
        // skip if dead wire or no hits
        if(WireHitRange[index].first < 0) continue;
        float prtime = et1 - es1;
        float bt1 = tcl[it1].BeginTim;
        if(bt1 < 0. || bt1 > maxtime) continue;
        unsigned short firsthit = WireHitRange[index].first;
        unsigned short lasthit = WireHitRange[index].second;
        // look for a hit within 10 ticks
        float best = 10.;
        short imbest = -1;
        for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
          float timediff = (plnhits[khit]->PeakTime() - prtime);
          // find the best time match, ignoring if the hit was already used
          if(fabs(timediff) < best) {
            best = fabs(timediff);
            imbest = khit;
          } // timediff
        } // khit
        if(imbest < 0) continue;
        if(hiterr2[imbest] < 0.) {
  if (fDebugHit < 0)  {
    std::cout<<"Hit US of cluster "<<tcl[it1].ID<<" W:H "<<uswire<<":"<<imbest;
    std::cout<<" is used in another cluster "<<std::endl;
  }
          // The hiterr2 is set < 0 which indicates that this hit is used in
          // a cluster. Look for it
          short inClus = -1;
          unsigned short atPos = 0;
          for(unsigned short it2 = 0; it2 < tcl.size(); ++it2) {
            if(it1 == it2) continue;
            if(tcl[it2].ID < 0) continue;
            for(unsigned short ih2 = 0; ih2 < tcl[it2].tclhits.size(); ++ih2) {
              if(imbest == tcl[it2].tclhits[ih2]) {
                atPos = ih2;
                inClus = it2;
                break;
              } // imbest == tcl...
            } // ih2
            if(inClus > 0) break;
          } // it2
          if(inClus > -1) {
  if(fDebugHit < 0) std::cout<<" in cluster "<<tcl[inClus].ID<<" at position "<<atPos<<std::endl;
            // found the hit on the US wire in a cluster. If atPos = 0,
            // these clusters are end-to-end but weren't merged. Neglect
            // these for now. Only consider the case where the hit is a bit
            // US of the Begin end so we can check the angle difference
            if(atPos > 4) {
              // compare the angle between these two clusters
              // start by doing a local fit of inClus
              unsigned int hit = tcl[inClus].tclhits[atPos];
              cl2FitMid(plnhits, tcl, inClus, hit, -3);
              if(clChisq > 90.) continue;
              float dth = fabs(atan(ScaleF * clpar[1]) - atan(ScaleF * es1));
  if(fDebugHit < 0)  std::cout<<" dth = "<<dth<<std::endl;
              // require a significant angle difference
              if(dth < 0.5) continue;
              // lop the end hits off of cluster inClus
              // first put the cluster hits into fcl2hits
              cl2TmpGet(plnhits, tcl, inClus);
              // The number of hits that we would like to lop off the End
              // of cluster inClus
              unsigned short nlop = tcl[inClus].tclhits.size() - atPos;
              for(unsigned short ii = 0; ii < nlop; ++ii) {
                fcl2hits.pop_back();
                if(fcl2hits.size() < 3) break;
              }
              if(fcl2hits.size() < 3) continue;
              // determine the pass for doing the fit
              pass = clProcCode - 10 * (clProcCode / 10) ;
              if(pass > fNumPass - 1) pass = 0;
              // re-fit the cluster
              cl2Fit(plnhits);
              clEndSlp = clpar[1];
              clEndSlpErr = clparerr[1];
              clProcCode += 200;
              // add the revised cluster
  if(fDebugHit < 0) {
    std::cout<<" Lopped hits off "<<tcl[inClus].ID<<" --> new cluster ";
    std::cout<<tcl[tcl.size()-1].ID<<std::endl;
  }
              cl2TmpStore(plnhits, tcl);
              // declare the old one obsolete
              tcl[inClus].ID = -tcl[inClus].ID;
            } // atPos < tcl[inClus].tclhits.size()
          } // inClus > -1
        } // hiterr2 < 0
      } // it1

//  cl2Print(plnhits, tcl);

    } // cl2VtxClusterSplit


/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2DoVertex(const art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx)
    {
      // try to make 2D vertices
      
      if(tcl.size() < 2) return;
      
      // initialize the begin and end vertex IDs
      for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
        tcl[ii].BeginVtx = -99;
        tcl[ii].EndVtx = -99;
      }
      
      unsigned int plane = plnhits[0]->WireID().Plane;
      float nwires = geom->Nwires(plane);
      float maxtime = detprop->NumberTimeSamples();

  if(fDebugHit < 0) std::cout<<"DoVertex plane "<<plane<<std::endl;

      for (unsigned short it1 = 0; it1 < tcl.size() - 1; ++it1) {
        // ignore abandoned clusters
        if(tcl[it1].ID < 0) continue;
        // ignore already attached clusters
        if(tcl[it1].BeginVtx >= 0 && tcl[it1].EndVtx >= 0) continue;
        float es1 = tcl[it1].EndSlp;
        unsigned short ew1 = tcl[it1].EndWir;
        float et1 = tcl[it1].EndTim;
        float bs1 = tcl[it1].BeginSlp;
        unsigned short bw1 = tcl[it1].BeginWir;
        float bt1 = tcl[it1].BeginTim;
        for (unsigned short it2 = it1 + 1; it2 < tcl.size(); ++it2) {
          // ignore abandoned clusters
          if(tcl[it2].ID < 0) continue;
          // ignore already attached clusters
          if(tcl[it2].BeginVtx >= 0 && tcl[it2].EndVtx >= 0) continue;
          // try to attach cluster to existing vertices at either end
          cl2ClsVertex(plnhits, tcl, vtx, it2);
          // ignore if both clusters are short
          if(tcl[it1].tclhits.size() < 10 &&
             tcl[it2].tclhits.size() < 10) continue;
          float es2 = tcl[it2].EndSlp;
          unsigned short ew2 = tcl[it2].EndWir;
          float et2 = tcl[it2].EndTim;
          float bs2 = tcl[it2].BeginSlp;
          unsigned short bw2 = tcl[it2].BeginWir;
          float bt2 = tcl[it2].BeginTim;
//  if(fDebugHit < 0) std::cout<<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID<<std::endl;
  // topo 1: check for vtx US of the ends of both clusters
          if(tcl[it1].EndVtx < 0 && tcl[it2].EndVtx < 0) {
            float dsl = es2 - es1;
            if(fabs(dsl) < 0.001) dsl = 0.001;
            // find vertex wire and vertex time in float precision (fvw, fvt)
            float fvw = 0.5 + (et1 - ew1 * es1 - et2 + ew2 * es2) / dsl;
            if(fvw > 0. && fvw < nwires) {
              // vertex wire in the detector
              unsigned short vw = fvw;
              // require vtx in the range of wires with hits AND
              // vtx US of both clusters AND
              // vtx not too far US of both clusters
              if(vw >= fFirstWire && 
                 vw <= ew1 && vw <= ew2 &&
                 vw  > ew1 - 10 && vw  > ew2 - 10) {
                float fvt = et1 + (vw - ew1) * es1;
  if(fDebugHit < 0) {
    std::cout<<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID;
    std::cout<<" topo1 vtx wire "<<vw<<" time "<<(int)fvt<<std::endl;
  }
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
              unsigned short vw = fvw;
              if(vw <= ew1 && vw >= bw2) {
                float fvt = et1 + (vw - ew1) * es1;
  if(fDebugHit < 0) {
    std::cout<<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID;
    std::cout<<" topo2 vtx wire "<<vw<<" time "<<(int)fvt<<std::endl;
  }
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
              unsigned short vw = fvw;
              if(vw <= ew2 && vw >= bw1) {
                float fvt = et2 + (vw - ew2) * es2;
  if(fDebugHit < 0) {
    std::cout<<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID;
    std::cout<<" topo3 vtx wire "<<vw<<" time "<<(int)fvt<<std::endl;
  }
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
              unsigned short vw = fvw;
              // require vtx in the range of wires with hits AND
              // vtx DS of both clusters AND
              // vtx not too far DS of both clusters
              if(vw >= bw1 && 
                 vw >= bw2 && vw <= fLastWire &&
                 vw <  bw2 + 10 && vw <  bw1 + 10) {
                float fvt = bt1 + (vw - bw1) * bs1;
  if(fDebugHit < 0) {
    std::cout<<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID;
    std::cout<<" topo4 vtx wire "<<vw<<" time "<<(int)fvt<<std::endl;
  }
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
      for(unsigned short it = 0; it < tcl.size(); ++it) {
        if(tcl[it].ID < 0) continue;
        // cluster weight = number of hits, saturated at 10
        float cw = tcl[it].tclhits.size();
        if(cw > 10) cw = 10;
        if(tcl[it].BeginVtx >=0) vtx[tcl[it].BeginVtx].Wght += cw;
        if(tcl[it].EndVtx >=0) vtx[tcl[it].EndVtx].Wght += cw;
      }
      

  if(fDebugHit < 0) {
    std::cout<<"Vertices "<<vtx.size()<<std::endl;
    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      std::cout<<"vtx "<<iv<<" wire "<<vtx[iv].Wire<<" time "<<(int)vtx[iv].Time<<" wght "<<(int)vtx[iv].Wght;
      std::cout<<" topo "<<vtx[iv].Topo<<std::endl;
    }
    cl2Print(plnhits, tcl);
  }
    }

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ClsVertex(const art::PtrVector<recob::Hit>& plnhits, 
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
        unsigned short it)
    {
      // try to attach cluster it2 to an existing vertex      
      if(vtx.size() == 0) return;
      
      for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
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
    }



/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ChkVertex(const art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
        short vw, float fvt, unsigned short it1, unsigned short it2, short topo)
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
        for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
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
          unsigned short iv = vtx.size() - 1;
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

      }

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ChkSignal(const art::PtrVector<recob::Hit>& plnhits,
      unsigned short wire1, float time1, unsigned short wire2, float time2, bool& SigOK)
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
      if(wiree < fFirstWire || wiree > fLastWire) {
            mf::LogError("ClusterCrawler")<<"ChkSignal bad wiree "<<wiree;
            return;
      }
      if(wireb < fFirstWire || wireb > fLastWire) {
            mf::LogError("ClusterCrawler")<<"ChkSignal bad wireb "<<wireb;
            return;
      }
      short wire0 = wiree;
      // checking a single wire?
      if(wireb == wiree) {
        clpar[1] = 0.;
      } else {
        clpar[1] = (timeb - timee) / (wireb - wiree);
      }
      for(short wire = wiree; wire < wireb + 1; ++wire) {
        // assume there is no signal on this wire
        bool WireSigOK = false;
        float prtime = timee + (wire - wire0) * clpar[1];
        // skip dead wires
        unsigned short index = wire - fFirstWire;
        if(WireHitRange[index].first < 0) continue;
        unsigned short firsthit = WireHitRange[index].first;
        unsigned short lasthit = WireHitRange[index].second;
//  std::cout<<"ChkSignal "<<wiree<<" "<<wire<<" "<<wireb<<" "<<(int)prtime;
//  std::cout<<" first "<<firsthit<<" last "<<lasthit;
        for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
          if(khit < 0 || khit > plnhits.size()-1) {
            mf::LogError("ClusterCrawler")<<"ChkSignal bad khit "<<khit;
            return;
          }
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
    }

/////////////////////////////////////////
/*
    void ClusterCrawlerAlg::cl2CurlyMerge(const art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl)
    {
      // try to merge short curly clusters. The algorithm uses the fact that
      // upstream clusters are added to the tcl struct after downstream 
      // clusters.
      
      if(tcl.size() < 2) return;
      
      unsigned short tclsize = tcl.size();
      for(unsigned short it1 = 0; it1 < tclsize - 1; ++it1) {
        if(tcl[it1].ID < 0) continue;
        // ignore long clusters
        if(tcl[it1].tclhits.size() > 20) continue;
        // ignore close pair clusters
        if(tcl[it1].ProcCode == 1000) continue;
        // ignore longish straight tracks
        unsigned short ew1 = tcl[it1].EndWir;
        float et1 = tcl[it1].EndTim;
        // we will look for a 2nd cluster on the next upstream wire = pw1
        unsigned short pw1 = ew1 - 1;
        // ignore dead wires for now
        unsigned short index = pw1 - fFirstWire;
        if(WireHitRange[index].first == -1) continue;
        // project the time to this wire using the two end hits of cluster 1
        std::vector<unsigned short>::reverse_iterator i4 = tcl[it1].tclhits.rbegin() + 1;
        unsigned short ih4 = *i4;
        float dt34 = plnhits[ih4]->PeakTime() - et1;
        float dw34 = plnhits[ih4]->WireID().Wire - ew1;
        // slope dT/dW
        float sl34 = dt34 / dw34;
        float th34 = atan(ScaleF * sl34);
        // projected time to the next wire (dW = 1)
        float pt1 = et1 - sl34;
        // use the fit array to hold the new(?) cluster
        fcl2hits = tcl[it1].tclhits;
        clBeginSlp = tcl[it1].BeginSlp;
        bool didit = false;
  if(prt) std::cout<<"Curly cl1 "<<tcl[it1].ID<<" pt1 "<<pt1<<std::endl;
        for(unsigned short it2 = it1 + 1; it2 < tclsize; ++it2) {
          if(tcl[it2].ID < 0) continue;
          if(tcl[it2].tclhits.size() > 20) continue;
          // ignore close pair clusters
          if(tcl[it2].ProcCode == 1000) continue;
          unsigned short bw2 = tcl[it2].BeginWir;
          // cluster 2 begins on wire US of the end of cluster 1?
          if(bw2 != pw1) continue;
          float bt2 = tcl[it2].BeginTim;
  if(prt) std::cout<<" >> cl2 "<<tcl[it2].ID<<" dt "<<fabs(bt2 - pt1)<<std::endl;
          // make a rough cut on the time difference
          if(fabs(bt2 - pt1) > 30) continue;
          didit = true;
          // find the begin angle of this cluster
          std::vector<unsigned short>::iterator i2 = tcl[it2].tclhits.begin() + 1;
          unsigned short ih2 = *i2;
          float dt12 = bt2 - plnhits[ih2]->PeakTime();
          short dw12 = bw2 - plnhits[ih2]->WireID().Wire;
          float sl12 = dt12 / dw12;
          float th12 = atan(ScaleF * sl12);
          float dth = fabs(th34 - th12);
  if(prt) std::cout<<" >> dth "<<dth<<" cut "<<fCurlyMergeAngCut<<std::endl;
          if(dth < fCurlyMergeAngCut) {
            // add it to fcl2hits
            for(std::vector<unsigned short>::iterator iht2 = tcl[it2].tclhits.begin();
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
            unsigned short index = pw1 - fFirstWire;
            if(WireHitRange[index].first == - 1) continue;
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
    } // cl2CurlyMerge
*/  
/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2SetBeginEnd(const art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl)
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
          // swap the hit order as well
/*
          std::vector<unsigned short> tvec;
          for(std::vector<short>::reverse_iterator jj = tcl[ii].tclhits.rbegin();
              jj != tcl[ii].tclhits.rend(); ++jj) {
            tvec.push_back(tcl[ii].tclhits[jj]);
          }
          tcl[ii].tclhits = tvec;
*/
        }
      }
    } // cl2SetBeginEnd

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ChkPair(const art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl)
    {
      // break clusters that are a part of a pair of clusters that
      // are less than PairAngCut
      
      if(tcl.size() < 2) return;
      // The size of the ClusterStore vector will increase after merging
      // is done.
      unsigned short tclsize = tcl.size();

      for (unsigned short it1 = 0; it1 < tclsize - 1; ++it1) {
        // ignore abandoned clusters
        if(tcl[it1].ID < 0) continue;
        // ignore short clusters
        if(tcl[it1].tclhits.size() < 10) continue;
        float bs1 = tcl[it1].BeginSlp;
        // convert slope to angle
        float bth1 = atan(ScaleF * bs1);
        unsigned short bw1 = tcl[it1].BeginWir;
        float bt1 = tcl[it1].BeginTim;
        unsigned short ew1 = tcl[it1].EndWir;
        for (unsigned short it2 = it1 + 1; it2 < tclsize; ++it2) {
          // ignore abandoned clusters
          if(tcl[it2].ID < 0) continue;
          if(tcl[it1].ID < 0) continue;
          // ignore short clusters
          if(tcl[it2].tclhits.size() < 10) continue;
          float bs2 = tcl[it2].BeginSlp;
          // convert slope to angle
          float bth2 = atan(ScaleF * bs2);
          unsigned short bw2 = tcl[it2].BeginWir;
          float bt2 = tcl[it2].BeginTim;
          unsigned short ew2 = tcl[it2].EndWir;
          if(fabs(bth1 - bth2) < fPairAngCut) {
            // check for begin angle difference
            if(ew2 > ew1 && ew2 < bw1) {
              // cluster 2 end is in the wire bounds of cluster 1
              // find the vertex position using the begin slope of both
              // clusters. 
              float dsl = bs2 - bs1;
              if(fabs(dsl) < 0.001) continue;
//  std::cout<<"Pair1 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<std::endl;
              unsigned short vw = (int)(0.5 + (bt1 - bw1 * bs1 - bt2 + bw2 * bs2) / dsl);
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
              unsigned short vw = (int)(0.5 + (bt2 - bw2 * bs2 - bt1 + bw1 * bs1) / dsl);
//  std::cout<<"Vtx wire "<<vw<<std::endl;
              // compare the vertex wire with the end of cluster 2
              if(abs(vw - ew2) < 5) cl2DoSplit(plnhits, tcl, it2, it1);
            } // begin/end wire test
          } // PairAngCut test
        } // it2
      } // it1
    }

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2DoSplit(const art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl, unsigned short it1, unsigned short it2)
    {
      // split cluster 1 at the gap in the middle of it. If no gap is found
      // split the cluster at the start of cluster 2

//  std::cout<<"DoSplit cl1 "<<tcl[it1].ID<<" cl2 "<<tcl[it2].ID<<std::endl;

      if(it1 > tcl.size() - 1) return;
      if(it2 > tcl.size() - 1) return;

      short splitwire = -1;
      short lastwire = -1;
      for(unsigned short iht1 = 0; iht1 < tcl[it1].tclhits.size(); ++iht1) {
        unsigned short hit = tcl[it1].tclhits[iht1];
        unsigned short wire = plnhits[hit]->WireID().Wire;
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
      for(unsigned short iht1 = 0; iht1 < tcl[it1].tclhits.size(); ++iht1) {
        unsigned short hit = tcl[it1].tclhits[iht1];
        if(hit > plnhits.size() - 1) {
          std::cout<<"cl2DoSplit bad hit "<<hit<<std::endl;
          return;
        }
        fcl2hits.push_back(hit);
        unsigned short wire = plnhits[hit]->WireID().Wire;
        if(wire == splitwire) {
          clEndTim = plnhits[hit]->PeakTime();
          break;
        }
      }
      // re-fit the end of the cluster
      cl2Fit(plnhits);
      clEndSlp = clpar[1];
      clEndSlpErr = clparerr[1];
      clEndWir = splitwire;
      clEndChg = fAveChg;
      cl2TmpStore(plnhits, tcl);
      // index of the first DownStream cluster
      unsigned short DScl1 = tcl.size() - 1;
      
      // look for an associated cluster with cluster1
      short jj = -1;
      for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
        if(tcl[ii].Assn == (int)it1) {
          jj = ii;
          break;
        }
      }
      unsigned short it3 = 9999;
//  std::cout<<"Assoc cluster? "<<jj<<std::endl;
      
      if(jj >= 0) {
        it3 = jj;
        // cluster 1 has an associated cluster embedded within it.
        // attach the hits on cluster 1 to the embedded cluster
        clStopCode = tcl[it3].StopCode;
        clProcCode = tcl[it3].ProcCode + 1000;
        fcl2hits = tcl[it3].tclhits;
        unsigned short wire = 0;
        for(unsigned short iht1 = 0; iht1 < tcl[it1].tclhits.size(); ++iht1) {
          unsigned short hit = tcl[it1].tclhits[iht1];
          if(hit > plnhits.size() - 1) {
            std::cout<<"cl2DoSplit bad hit "<<hit<<std::endl;
            return;
          }
          wire = plnhits[hit]->WireID().Wire;
          if(wire < splitwire) fcl2hits.push_back(hit);
        }
        cl2Fit(plnhits);
        clEndSlp = clpar[1]; // save the slope at the end
        clEndSlpErr = clparerr[1];
        clEndChg = fAveChg;
        cl2TmpStore(plnhits, tcl);
      } else {
        // create a new cluster using the hits on cluster 1 that are
        // upstream of splitwire. 
        clStopCode = tcl[it1].StopCode;
        clProcCode = tcl[it1].ProcCode + 1000;
        fcl2hits.clear();
        bool didfit = false;
        unsigned short wire = 0;
        for(unsigned short iht1 = 0; iht1 < tcl[it1].tclhits.size(); ++iht1) {
          unsigned short hit = tcl[it1].tclhits[iht1];
          wire = plnhits[hit]->WireID().Wire;
          if(wire < splitwire) {
            fcl2hits.push_back(hit);
            // re-fit the begin end when there are enough hits
            if(fcl2hits.size() == 3) {
//  std::cout<<"refit begin "<<tcl[it1].ID<<std::endl;
              didfit = true;
              cl2Fit(plnhits);
              if(clChisq > 99.) {
  std::cout<<"cl2DoSplit bad fit "<<std::endl;
                return;
              }
              clBeginSlp = clpar[1];
              clBeginSlpErr = clparerr[1];
            } // fcl2hits test
          } // wire < splitwire
        } // iht1 iterator
        if(didfit) {
          // did a fit at the beginning of the cluster
          if(fcl2hits.size() > 4) {
            // long cluster: re-fit the end
            cl2Fit(plnhits);
              if(clChisq > 99.) {
  std::cout<<"cl2DoSplit bad fit "<<std::endl;
                return;
              }
            clEndSlp = clpar[1];
            clEndSlpErr = clparerr[1];
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
      unsigned short UScl = tcl.size() - 1;
      // The first downstream cluster was defined above
      // the second downstream cluster is it2
      unsigned short DScl2 = it2;

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
      for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
        if(tcl[ii].ID > 0) index.push_back(std::make_pair(tcl[ii].tclhits.size(),ii));
      }
      std::sort(index.begin(), index.end(), SortByLen);
      sortindex.clear();
      for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
        sortindex.push_back(index[ii].second);
//    std::cout<<"sortd "<<index[ii].first<<" "<<index[ii].second<<std::endl;
      }
      return; 
    }
*/
/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ChkMerge(const art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl)
    {
      // Try to merge clusters. Clusters that have been subsumed in other
      // clusters, i.e. no longer valid, have ID < 0
      
      if(tcl.size() < 2) return;
      // The size of the ClusterStore vector will increase while merging
      // is done.
      
      unsigned short tclsize = tcl.size();
            
      for(unsigned short it1 = 0; it1 < tclsize - 1; ++it1) {
        // ignore already merged clusters
        if(tcl[it1].ID < 0) continue;
        float bs1 = tcl[it1].BeginSlp;
        // convert slope to angle
        float arg = ScaleF * bs1;
        float bth1 = atan(arg);
        // error on the angle
        float bth1e = ScaleF * tcl[it1].BeginSlpErr / (1 + arg * arg);
        // more compact notation for begin/end, wire/time/chg/slp/theta, 1/2
        unsigned short bw1 = tcl[it1].BeginWir;
        float bt1 = tcl[it1].BeginTim;
        float bc1 = tcl[it1].BeginChg;
        float es1 = tcl[it1].EndSlp;
        // convert slope to angle
        arg = ScaleF * es1;
        float eth1 = atan(arg);
        float eth1e =  ScaleF * tcl[it1].EndSlpErr / (1 + arg * arg);
        unsigned short ew1 = tcl[it1].EndWir;
        float et1 = tcl[it1].EndTim;
        float ec1 = tcl[it1].EndChg;
        unsigned short pass1 = tcl[it1].ProcCode - 10 * (tcl[it1].ProcCode / 10);
        for(unsigned short it2 = it1 + 1; it2 < tclsize; ++it2) {
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
          unsigned short bw2 = tcl[it2].BeginWir;
          float bt2 = tcl[it2].BeginTim;
          float bc2 = tcl[it2].BeginChg;
          float es2 = tcl[it2].EndSlp;
          arg = ScaleF * es2;
          float eth2 = atan(arg);
          float eth2e = ScaleF * tcl[it2].EndSlpErr / (1 + arg * arg);
          unsigned short ew2 = tcl[it2].EndWir;
          float et2 = tcl[it2].EndTim;
          float ec2 = tcl[it2].EndChg;
          unsigned short pass2 = tcl[it2].ProcCode - 10 * (tcl[it2].ProcCode / 10);
          // use the more promiscuous pass for cuts
          float angcut = fKinkAngCut[pass1];
          if(fKinkAngCut[pass2] > angcut) angcut = fKinkAngCut[pass2];
          unsigned short skipcut = fMaxWirSkip[pass1];
          if(fMaxWirSkip[pass2] > skipcut) skipcut = fMaxWirSkip[pass2];
          float chgcut = fChgCut[pass1];
          if(fChgCut[pass2] > chgcut) chgcut = fChgCut[pass2];
          float timecut = fTimeDelta[pass];
          if(fTimeDelta[pass2] > timecut) timecut = fTimeDelta[pass2];
          // increase the time cut for large angle clusters
          timecut *= (2 - 1/(1 + fabs(clpar[1])));
          
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
    }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2ChkMerge12(const art::PtrVector<recob::Hit>& plnhits, 
     std::vector<ClusterStore>& tcl, unsigned short it1, unsigned short it2, bool& didit)
  {
    // Calling routine has done a rough check that cluster it2 is a candidate
    // for merging with cluster it1. The wire range spanned by it2 lies
    // within the wire range of it1 and the clusters are reasonably close
    // together in time.
    
    // assume that no merging was done
    didit = false;
    
  if(fDebugWire < 0) std::cout<<"cl2ChkMerge12 "<<tcl[it1].ID<<" "<<tcl[it2].ID<<std::endl;
    
    ClusterStore& cl1 = tcl[it1];
    // fill a vector spanning the length of cluster 1 and filled with the hit time
    unsigned short ew1 = tcl[it1].EndWir;
    unsigned short bw1 = tcl[it1].BeginWir;
    std::vector<unsigned short> cl1hits;
    // fill the vector with 0s
    for(unsigned short wire = ew1; wire <= bw1; ++wire) {
      cl1hits.push_back(0);
    }
    // put in the hit IDs
    for(unsigned short iht = 0; iht < cl1.tclhits.size(); ++iht) {
      unsigned short hit = cl1.tclhits[iht];
      unsigned short wire = plnhits[hit]->WireID().Wire;
      if(wire - ew1 < 0 || wire - ew1 > (short)cl1hits.size()) {
        mf::LogError("ClusterCrawler")<<"ChkMerge12 bad wire "<<(wire-ew1);
        return;
      }
      cl1hits[wire - ew1] = hit;
    }
    unsigned short ew2 = tcl[it2].EndWir;
    float et2 = tcl[it2].EndTim;
    // look for the closest wire with a hit on cluster 1
    unsigned short wiron1 = 0;
    // count the number of missing hits
    short nmiss = 0;
    for(unsigned short wire = ew2 - 1; wire > ew1; --wire) {
      if(cl1hits[wire - ew1] > 0) {
        wiron1 = wire;
        break;
      }
      ++nmiss;
    } // wire
  if(fDebugWire < 0) std::cout<<"chk next US wire "<<wiron1<<" missed "<<nmiss<<std::endl;
    if(wiron1 == 0) return;
    if(nmiss > fMaxWirSkip[pass]) return;
    
    // compare the wires with hits on cluster 2 with the gap in cluster 1
    // the number of hit wires that fit in the gap
    unsigned short nfit = 0;
    for(unsigned short iht = 0; iht < tcl[it2].tclhits.size(); ++iht) {
      unsigned short hiton2 = tcl[it2].tclhits[iht];
      unsigned short wiron2 = plnhits[hiton2]->WireID().Wire;
      if(wiron2 < ew1 || wiron2 > bw1) return;
      if(cl1hits[wiron2 - ew1] == 0) ++nfit;
    }
    // require complete filling of the gap
    if(nfit < tcl[it2].tclhits.size()) return;
    
    // decode the pass for both clusters and select the matching cuts
    unsigned short pass1 = tcl[it1].ProcCode - 10 * (tcl[it1].ProcCode / 10);
    unsigned short pass2 = tcl[it2].ProcCode - 10 * (tcl[it2].ProcCode / 10);
    unsigned short cpass = pass1;
    // use the tighter cuts
    if(pass2 < pass1) cpass = pass2;
    
    // ***** Check End of Cluster 2 matching with middle of cluster 1
    if(wiron1 - ew1 < 0) return;
    unsigned short hiton1 = cl1hits[wiron1 - ew1];
    if(hiton1 > plnhits.size() - 1) {
      mf::LogError("ClusterCrawler")<<"ChkMerge12 bad hiton1 "<<hiton1;
      return;
    }
    // check the time difference
    float timon1 = plnhits[hiton1]->PeakTime();
    float dtim = fabs(et2 + (wiron1 - ew2) * tcl[it2].EndSlp - timon1);
    if(dtim > fTimeDelta[cpass]) return;
    // check the slope difference. First do a local fit on cluster 1 near
    // the matching point
    cl2FitMid(plnhits, tcl, it1, hiton1, 3);
    if(clChisq > 20.) return;
    // fit parameters are now in clpar. Charge is in fAveChg
    // check for angle consistency
    float dth = fabs(atan(ScaleF * clpar[1]) - atan(ScaleF * tcl[it2].EndSlp));
  if(fDebugWire < 0) std::cout<<"US dtheta "<<dth<<" cut "<<fKinkAngCut[cpass]<<std::endl;
    if(dth > fKinkAngCut[cpass]) return;
    // make a charge ratio cut
    float chgrat = 2 * fabs(fAveChg - tcl[it2].EndChg) / (fAveChg + tcl[it2].EndChg);
  if(fDebugWire < 0)  std::cout<<"US chgrat "<<chgrat<<" cut "<<fChgCut[pass]<<std::endl;
    // ensure that there is a signal on any missing wires at the US end of 1
    bool SigOK = false;
    cl2ChkSignal(plnhits, wiron1, timon1, ew2, et2, SigOK);
  if(fDebugWire < 0)  std::cout<<"US SigOK? "<<SigOK<<std::endl;
    if( !SigOK ) return;


    // ***** Check Begin of Cluster 2 matching with middle of cluster 1
    unsigned short bw2 = tcl[it2].BeginWir;
    float bt2 = tcl[it2].BeginTim;
    nmiss = 0;
    wiron1 = 0;
    for(unsigned short wire = bw2 + 1; wire < bw1; ++wire) {
      if(cl1hits[wire - ew1] > 0) {
        wiron1 = wire;
        break;
      }
      ++nmiss;
    }
    if(wiron1 == 0) return;
    if(nmiss > fMaxWirSkip[pass]) return;
    // fit this section of cluster 1 with 4 hits starting at the hit on the
    // closest wire and moving DS
    hiton1 = cl1hits[wiron1 - ew1];
    if(hiton1 > plnhits.size() - 1) {
      mf::LogError("ClusterCrawler")<<"ChkMerge12 bad hiton1 "<<hiton1;
      return;
    }
    timon1 = plnhits[hiton1]->PeakTime();
    dtim = fabs(bt2 - (wiron1 - bw2) *tcl[it2].BeginSlp - timon1);
    if(dtim > fTimeDelta[cpass]) return;
    cl2FitMid(plnhits, tcl, it1, hiton1, -3);
    if(clChisq > 20.) return;
    // check for angle consistency
    dth = fabs(atan(ScaleF * clpar[1]) - atan(ScaleF * tcl[it2].BeginSlp));
  if(fDebugWire < 0) std::cout<<"DS dtheta "<<dth<<" cut "<<fKinkAngCut[cpass]<<std::endl;
    if(dth > fKinkAngCut[cpass]) return;
    // make a charge ratio cut
    chgrat = 2 * fabs(fAveChg - tcl[it2].BeginChg) / (fAveChg + tcl[it2].BeginChg);
  if(fDebugWire < 0)  std::cout<<"DS chgrat "<<chgrat<<" cut "<<fChgCut[pass]<<std::endl;
    // ensure that there is a signal on any missing wires at the US end of 1
    SigOK = false;
    cl2ChkSignal(plnhits, wiron1, timon1, bw2, bt2, SigOK);
  if(fDebugWire < 0)  std::cout<<"DS SigOK? "<<SigOK<<std::endl;
    if( !SigOK ) return;

  if(fDebugWire < 0)  std::cout<<"Merge em"<<std::endl;
    // success. Merge them
    cl2DoMerge(plnhits, tcl, it1, it2, 100);
    didit = true;
  }


/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2DoMerge(const art::PtrVector<recob::Hit>& plnhits, 
     std::vector<ClusterStore>& tcl, unsigned short it1, unsigned short it2,
     short ProcCode)
  {
    // Merge clusters. Cluster 1 has precedence for assignment of hits
    
    ClusterStore& cl1 = tcl[it1];
    ClusterStore& cl2 = tcl[it2];
    // find the low and high wire for both clusters
    unsigned short hiwire = cl1.BeginWir;
    if(cl2.BeginWir > hiwire) hiwire = cl2.BeginWir;
    if(hiwire > fLastWire) {
      throw cet::exception("ClusterCrawler")<<"DoMerge bad hiwire "<<hiwire;
      return;
    }
    unsigned short lowire = cl1.EndWir;
    if(cl2.EndWir < lowire) lowire = cl2.EndWir;
    if(lowire < fFirstWire) {
      throw cet::exception("ClusterCrawler")<<"DoMerge bad lowire "<<lowire;
      return;
    }
    // make a vector of wire hits
    std::vector<short> wirehit;
    for(unsigned short wire = lowire; wire < hiwire + 2; ++wire) {
      wirehit.push_back(-1);
    }
    unsigned short veclen = hiwire - lowire + 1;
//  std::cout<<"lo/hi "<<lowire<<" "<<hiwire<<" veclen "<<veclen<<std::endl;
//  std::cout<<"tcl size "<<cl1.tclhits.size()<<" "<<cl2.tclhits.size()<<std::endl;
    // put in the hit IDs for cluster 1
    for(unsigned short iht = 0; iht < cl1.tclhits.size(); ++iht) {
      unsigned short hit = cl1.tclhits[iht];
      unsigned short wire = plnhits[hit]->WireID().Wire;
      unsigned short index = wire - lowire;
      if(index > veclen) {
        throw cet::exception("ClusterCrawler")<<"DoMerge bad index "<<index;
        return;
      }
      wirehit[index] = hit;
    }
    // now cluster 2
//  std::cout<<"cl2 "<<std::endl;
    for(unsigned short iht = 0; iht < cl2.tclhits.size(); ++iht) {
      unsigned short hit = cl2.tclhits[iht];
      unsigned short wire = plnhits[hit]->WireID().Wire;
      unsigned short index = wire - lowire;
      if(index > veclen) {
        throw cet::exception("ClusterCrawler")<<"DoMerge bad index "<<index;
        return;
      }
      if(wirehit[index] < 0) {
        wirehit[index] = hit;
      } else {
        // a hit from cluster cl1 is on this wire. Free up cl2 hit for later use
        hiterr2[hit] = fabs(hiterr2[hit]);
      }
    }
    // make the new cluster
    fcl2hits.clear();
    for(std::vector<short>::reverse_iterator ii = wirehit.rbegin();
        ii != wirehit.rend(); ++ii) {
      short hit = *ii;
      if(hit >= 0) {
        unsigned short uhit = hit;
        fcl2hits.push_back(uhit);
        if(fcl2hits.size() == 4) {
          // re-fit the Begin end of the cluster
          cl2Fit(plnhits);
          if(clChisq > 99.) {
            std::cout<<"cl2DoMerge bad Begin fit "<<clChisq<<std::endl;
            return;
          }
          clBeginSlp = clpar[1];
          clBeginSlpErr = clparerr[1];
          clBeginChg = fAveChg;
        }
      }
    }
//  std::cout<<"new size "<<fcl2hits.size()<<std::endl;
    // re-fit the End of the cluster with the current pass params
    cl2Fit(plnhits);
    if(clChisq > 99.) {
      std::cout<<"cl2DoMerge bad End fit "<<clChisq<<std::endl;
      return;
    }
    clEndSlp = clpar[1];
    clEndSlpErr = clparerr[1];
    clEndChg = fAveChg;

    clStopCode = cl1.StopCode;
    clAssn = -1;
    // append it to the tcl vector
    cl2TmpStore(plnhits, tcl);
    unsigned short itnew = tcl.size()-1;
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
    for(unsigned short ii = 0; ii < itnew; ++ii) {
      if(tcl[ii].ID > 0 && tcl[ii].Assn >= 0) {
        if(tcl[ii].Assn == (short)it1 || tcl[ii].Assn == (short)it2) tcl[ii].Assn = itnew;
      }
    }
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2Print(const art::PtrVector<recob::Hit>& plnhits, 
     std::vector<ClusterStore>& tcl)
  {
    // prints clusters to the screen for code development
    std::cout<<"  ID nht Stop  Proc   beg_W:H   begT  bTheta Therr begChg  end_W:H  endT  eTheta Therr endChg";
    std::cout<<" bVx eVx";
    std::cout<<std::endl;
    for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
      std::vector<unsigned short>::const_iterator ihtb = tcl[ii].tclhits.begin();
      unsigned short hitb = *ihtb;
      std::vector<unsigned short>::const_iterator ihte = tcl[ii].tclhits.end()-1;
      unsigned short hite = *ihte;
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
  }

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2TmpGet(const art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl, unsigned short it1)
    {
      // copies temp cluster it1 into the fcl2hits vector, etc. This is 
      // effectively the inverse of cl2TmpStore
      
      if(it1 > tcl.size()) return;


      clBeginSlp = tcl[it1].BeginSlp;
      clBeginSlpErr = tcl[it1].BeginSlpErr;
      clBeginWir = tcl[it1].BeginWir;
      clBeginTim = tcl[it1].BeginTim;
      clBeginChg = tcl[it1].BeginChg;
      clEndSlp = tcl[it1].EndSlp;
      clEndSlpErr = tcl[it1].EndSlpErr;
      clEndWir = tcl[it1].EndWir;
      clEndTim = tcl[it1].EndTim;
      clEndChg = tcl[it1].EndChg;
      clStopCode = tcl[it1].StopCode;
      clProcCode = tcl[it1].ProcCode;
      clAssn = tcl[it1].Assn;
      fcl2hits.clear();
      fcl2hits = tcl[it1].tclhits;
    }


/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2TmpStore(const art::PtrVector<recob::Hit>& plnhits, 
     std::vector<ClusterStore>& tcl)
  {

    if(fcl2hits.size() < 3) {
      mf::LogError("ClusterCrawler")<<"cl2TmpStore trying to store crazy cluster";
      return;
    }
    
    ++NClusters;

    // flag all the hits as used
    for(unsigned short it = 0; it < fcl2hits.size(); ++it) {
      unsigned short hit = fcl2hits[it];
      if(hit > plnhits.size() - 1) {
        std::cout<<"cl2TmpStore bad hit "<<hit<<std::endl;
        return;
      }
      hiterr2[hit] = -fabs(hiterr2[hit]);
    }

    // ensure that the cluster begin/end info is correct

    // define the begin/end charge if it wasnt done already
    if(clEndChg < 0.) {
      // use the next to the last two hits. The End hit may have low charge
      unsigned int ih0 = fcl2hits.size() - 2;
      unsigned int hit = fcl2hits[ih0];
      clEndChg = plnhits[hit]->Charge();
      hit = fcl2hits[ih0 - 1];
      clEndChg += plnhits[hit]->Charge();
      clEndChg = clEndChg / 2.;
    }
    if(clBeginChg < 0.) {
      // use the 2nd and third hit. The Begin hit may have low charge
      unsigned int hit = fcl2hits[1];
      clBeginChg = plnhits[hit]->Charge();
      hit = fcl2hits[2];
      clBeginChg += plnhits[hit]->Charge();
      clBeginChg = clBeginChg / 2.;
    }
    
    std::vector<unsigned short>::const_iterator ibg = fcl2hits.begin();
    unsigned short hitb = *ibg;
    std::vector<unsigned short>::const_iterator iend = fcl2hits.end() - 1;
    unsigned short hite = *iend;

    // store the cluster in the temporary ClusterStore struct
    ClusterStore clstr;
    
    clstr.ID = NClusters;
    clstr.BeginSlp    = clBeginSlp;
    clstr.BeginSlpErr = clBeginSlpErr;
    clstr.BeginWir    = plnhits[hitb]->WireID().Wire;
    clstr.BeginTim    = plnhits[hitb]->PeakTime();
    clstr.BeginChg    = clBeginChg;
    clstr.EndSlp      = clEndSlp;
    clstr.EndSlpErr   = clEndSlpErr;
    clstr.EndWir      = plnhits[hite]->WireID().Wire;
    clstr.EndTim      = plnhits[hite]->PeakTime();
    clstr.EndChg      = clEndChg;
    clstr.StopCode    = clStopCode;
    clstr.ProcCode    = clProcCode;
    clstr.Assn        = clAssn;
    clstr.BeginVtx    = -99;
    clstr.EndVtx      = -99;
    clstr.tclhits     = fcl2hits;
    tcl.push_back(clstr);
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2FollowUS(const art::PtrVector<recob::Hit>& plnhits)
  {
    // follow the cluster upstream

    if(fcl2hits.size() < 2) return;

    std::vector<unsigned short>::const_iterator itt = fcl2hits.begin();
    unsigned short fhit = *itt;
    unsigned short fwir = plnhits[fhit]->WireID().Wire;
  if(fDebugWire > 0 && fDebugHit > 0) {
    prt = ((short)fwir == fDebugWire && (short)fhit == fDebugHit);
  }

  if(prt) {
    std::cout<<"cl2FollowUS PASS "<<pass<<" Hits: ";
    for(std::vector<unsigned short>::reverse_iterator itm = fcl2hits.rbegin(); itm != fcl2hits.rend(); ++itm) {
      std::cout<<*itm<<" ";
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
    unsigned short it = fcl2hits.size() - 1;
    unsigned short lasthit = fcl2hits[it];
    if(lasthit > plnhits.size() - 1) {
      std::cout<<"cl2FollowUS bad lasthit "<<lasthit<<std::endl;
    }
    unsigned short lastwire = plnhits[lasthit]->WireID().Wire;
  if(prt) std::cout<<"cl2FollowUS: last wire "<<lastwire<<" hit "<<lasthit<<std::endl;
    // keep a log of the fit chisq
    std::vector<float> chifits;
    
    for(unsigned short nextwire = lastwire-1; nextwire >= fFirstWire; --nextwire) {
  if(prt) std::cout<<"cl2FollowUS: next wire "<<nextwire<<std::endl;
      // add hits and check for PH and width consistency
      cl2AddHit(plnhits, nextwire, HitOK, SigOK);
  if(prt) std::cout<<"cl2FollowUS: HitOK "<<HitOK<<" SigOK "<<SigOK<<std::endl;
      if(!HitOK) {
        // no hit on this wire. Was there a signal or dead wire?
        if(SigOK) {
          ++nmissed;
          if(prt && nmissed > fMaxWirSkip[pass]) std::cout<<"nmissed break"<<std::endl;
          if(nmissed > fMaxWirSkip[pass]) {
            clStopCode = 1;
            break;
          }
          // see if we are in the PostSkip phase and missed more than 1 wire
//  std::cout<<"chk "<<nextwire<<" "<<PostSkip<<" "<<nmissed<<std::endl;
          if(PostSkip && nmissed > 1) {
//  std::cout<<"Missed wire after skip "<<nextwire<<" nmissed "<<nmissed<<std::endl;
            if((short)(fcl2hits.size() - nHitAfterSkip) < 4) {
              fcl2hits.clear();
              return;
            }
  if(prt) std::cout<<" PostSkip && nmissed = "<<nmissed<<std::endl;
            for(short jj = 0; jj < nHitAfterSkip; ++jj) {
              fcl2hits.pop_back();
            } // pop_back
            cl2Fit(plnhits);
            if(clChisq > 90.) {
              fcl2hits.clear();
              return;
            }
            clStopCode = 2;
            return;
          } // PostSkip && nmissed > 
          if(nmissed > 1) {
            DidaSkip = true;
            PostSkip = false;
          }
        } else {
          // SigOK is false
          clStopCode = 0;
          if(prt) std::cout<<"No hit or signal on wire "<<nextwire<<std::endl;
//  std::cout<<"No hit or signal on wire "<<nextwire<<std::endl;
          break;
        } // else SigOK false
      } else {
        // HitOK is true. Update the fit
        cl2Fit(plnhits);
        if(clChisq > 99.) {
          if(fcl2hits.size() < 3) return;
          // a fit error occurred. Lop off the leading hit and refit
  if(prt) std::cout<<"Fit failed "<<std::endl;
          fcl2hits.pop_back();
          cl2Fit(plnhits);
          if(clChisq > 99.) {
            // something really bad happened. Bail out
            fcl2hits.clear();
            return;
          }
          continue;
        }
        // monitor the onset of a kink. Find the average chisq for the fit
        // using the previous 3 - 6 hits. Look for a progressive increase
        // in chisq for the previous 0 - 2 hits.
        chifits.push_back(clChisq);
        if(chifits.size() > 7 && fKinkChiRat[pass] > 0) {
          unsigned short chsiz = chifits.size()-1;
  if(prt) {
    std::cout<<"Kink chk "<<chifits[chsiz]<<" "<<chifits[chsiz-1]<<" ";
    std::cout<<chifits[chsiz-2]<<" "<<chifits[chsiz-3]<<std::endl;
  }
          if( chifits[chsiz-2] > fKinkChiRat[pass] * chifits[chsiz-3] && 
              chifits[chsiz-1] > fKinkChiRat[pass] * chifits[chsiz-2] &&
              chifits[chsiz]   > fKinkChiRat[pass] * chifits[chsiz-1]) {
            if(fcl2hits.size() < 8) {
    std::cout<<"bad kink check size "<<chifits.size()<<" "<<fcl2hits.size()<<std::endl;
              continue;
            }
            // find the kink angle (crudely) from the 0th and 2nd hit
            unsigned short ih0 = fcl2hits.size() - 1;
            unsigned short hit0 = fcl2hits[ih0];
            if(hit0 > plnhits.size()-1) {
              mf::LogError("ClusterCrawler")<<"FollowUS bad hit0 "<<hit0;
              return;
            }
            unsigned short ih2 = ih0 - 2;
            unsigned short hit2 = fcl2hits[ih2];
            if(hit2 > plnhits.size()-1) {
              mf::LogError("ClusterCrawler")<<"FollowUS bad hit2 "<<hit2;
              return;
            }
            float dt02 = plnhits[hit2]->PeakTime() - plnhits[hit0]->PeakTime();
            float dw02 = plnhits[hit2]->WireID().Wire - plnhits[hit0]->WireID().Wire;
            float th02 = atan( ScaleF * dt02 / dw02);
            // and the 3rd and 5th hit
            unsigned short ih3 = ih2 - 1;
            unsigned short hit3 = fcl2hits[ih3];
            if(hit3 > plnhits.size()-1) {
              mf::LogError("ClusterCrawler")<<"FollowUS bad hit3 "<<hit3;
              return;
            }
            unsigned short ih5 = ih3 - 2;
            unsigned short hit5 = fcl2hits[ih5];
            if(hit5 > plnhits.size()-1) {
              mf::LogError("ClusterCrawler")<<"FollowUS bad hit5 "<<hit5;
              return;
            }
            float dt35 = plnhits[hit5]->PeakTime() - plnhits[hit3]->PeakTime();
            float dw35 = plnhits[hit5]->WireID().Wire - plnhits[hit3]->WireID().Wire;
            float th35 = atan(ScaleF * dt35 / dw35);
            float dth = fabs(th02 - th35);
  if(prt) std::cout<<" Kink angle "<<std::setprecision(3)<<dth<<std::endl;
            // cut on the allowed kink angle
            if(dth > fKinkAngCut[pass]) {
  if(prt) std::cout<<"stopped tracking "<<std::endl;
              // kill the last 3 hits, refit and return
              for(short jj = 0; jj < 3; ++jj) {
                fcl2hits.pop_back();
              }
              cl2Fit(plnhits);
              clStopCode = 3;
              break;
            } // kinkang check
          } // chifits check
        } // chifits.size() > 5
        // done with kink check
        // update the cluster Begin information?
        if(fcl2hits.size() == fMaxHitsFit[pass] ||
           fcl2hits.size() == fMinHits[pass]) {
          clBeginSlp = clpar[1];
          clBeginSlpErr = clparerr[1];
        }
        // set the Begin charge after fAveChg is defined
        if(fAveChg > 0 && clBeginChg < 0) {
          // project the charge to the Begin of the cluster
          clBeginChg = fAveChg + (clBeginWir - nextwire) * fChgSlp;
  if(prt) std::cout<<"Set clBeginChg "<<clBeginChg<<std::endl;
        }
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
          ++nHitAfterSkip;
          if(nHitAfterSkip == fMinWirAfterSkip[pass]) PostSkip = false;
        } 
        // check for bad chisq
        if(clChisq > fChiCut[pass]) {
          // remove the last hit and re-fit
          fcl2hits.pop_back();
          cl2Fit(plnhits);
          if(clChisq > 99.) {
            std::cout<<"cl2FollowUS bad fit after pop_back "<<clChisq<<std::endl;
            return;
          }
          clStopCode = 4;
          break;
        } // clChisq > fChiCut[pass]
      } // !HitOK check
    } // nextwire
    
    // find the US wire
    unsigned short ih0 = fcl2hits.size() - 1;
    unsigned short hit0 = fcl2hits[ih0];
    unsigned short uswir = plnhits[hit0]->WireID().Wire;
    // find the wire fMinWirAfterSkip[pass] indices DS
    unsigned short ihskp = ih0 - fMinWirAfterSkip[pass];
    unsigned short hitskp = fcl2hits[ihskp];
    unsigned short wirskp = plnhits[hitskp]->WireID().Wire;
    unsigned short nlop = wirskp - uswir - fMinWirAfterSkip[pass];

  if(prt) std::cout<<" check nlop "<<nlop<<" wirskp "<<wirskp<<std::endl;
    if(nlop > 0 && nlop < fcl2hits.size() - 4) {
      for(unsigned short ii = 0; ii < nlop; ++ii) {
        fcl2hits.pop_back();
        if(fcl2hits.size() < 3) {
          std::cout<<"cl2FollowUS bad pop_back "<<std::endl;
          return;
        }
      }
      clStopCode = 2;
    }
    prt = false;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2FitMid(const art::PtrVector<recob::Hit>& plnhits,
      std::vector<ClusterStore>& tcl, unsigned short it1, unsigned short ihtin, short nhit)
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

    std::vector<float> xwir;
    std::vector<float> ytim;
    std::vector<float> ytimerr2;
    
    short nht = 0;
    unsigned short wire0 = 0;
    if(nhit > 0) {
      nht = nhit;
      // find the first desired hit and move towards the End
      fAveChg = 0.;
      unsigned short hitcnt = 0;
      bool UseEm = false;
      for(unsigned short it = 0; it < cls.tclhits.size(); ++it) {
        unsigned short ihit = cls.tclhits[it];
        if(ihit > plnhits.size()-1) {
          mf::LogError("ClusterCrawler")<<"FitMid bad ihit "<<ihit;
          return;
        }
        // look for the desired first hit. Use this as the origin wire
        if(ihit == ihtin) {
          UseEm = true;
          wire0 = plnhits[ihit]->WireID().Wire;
        }
        // use hits after finding the first desired hit
        if(UseEm) {
          unsigned short wire = plnhits[ihit]->WireID().Wire;
          xwir.push_back(wire - wire0);
          ytim.push_back(plnhits[ihit]->PeakTime());
          // pass the error^2 to the fitter
	  ytimerr2.push_back(fabs(hiterr2[ihit]));
          fAveChg += plnhits[ihit]->Charge();
          ++hitcnt;
          if(hitcnt == nht) break;
        }
      }    
      nht = hitcnt;
    } else {
      nht = -nhit;
      // find the first desired hit and move towards the Begin
      fAveChg = 0.;
      unsigned short hitcnt = 0;
      bool UseEm = false;
      for(unsigned short it = cls.tclhits.size() - 1; it >= 0; it--) {
        unsigned short ihit = cls.tclhits[it];
        if(ihit > plnhits.size()-1) {
          mf::LogError("ClusterCrawler")<<"FitMid bad ihit "<<ihit;
          return;
        }
        // look for the desired first hit. Use this as the origin wire
        if(ihit == ihtin) {
          UseEm = true;
          wire0 = plnhits[ihit]->WireID().Wire;
        }
        // use hits after finding the first desired hit

        if(UseEm) {
          unsigned short wire = plnhits[ihit]->WireID().Wire;
          xwir.push_back(wire - wire0);
          ytim.push_back(plnhits[ihit]->PeakTime());
	  ytimerr2.push_back(fabs(hiterr2[ihit]));
          fAveChg += plnhits[ihit]->Charge();
          ++hitcnt;
          if(hitcnt == nht) break;
        }
      }    
      nht = hitcnt;
    }
    
    if(nht < 2) return;
    
    float intcpt = 0.;
    float slope = 0.;
    float intcpterr = 0.;
    float slopeerr = 0.;
    float chidof = 0.;
    LinFit(xwir, ytim, ytimerr2, intcpt, slope, intcpterr, slopeerr, chidof);
    clChisq = chidof;
    if(clChisq > 99.) return;
    clpar[0] = intcpt;
    clpar[1] = slope;
    clparerr[0] = intcpterr;
    clparerr[1] = slopeerr;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2Fit(const art::PtrVector<recob::Hit>& plnhits)
  {
    // Fits the hits on the US end of a cluster. This routine assumes that
    // wires are numbered from lower (upstream) to higher (downstream) and
    // that the hits in the fclhits vector are sorted so that upstream hits
    // are at the end of the vector


    clChisq = 999.;
    
    unsigned short nht = 0;
    // fit all hits or truncate?
    if(fcl2hits.size() < fMaxHitsFit[pass]) {
      nht = fcl2hits.size();
    } else {
      nht = fMaxHitsFit[pass];
    }
    if(nht < 2) return;

    std::vector<float> xwir;
    std::vector<float> ytim;
    std::vector<float> ytimerr2;
    // apply an angle dependent scale factor. The error should be
    // wire pitch / sqrt(12) for a cluster at 90 degrees. This formula 
    // simply doubles the error, which I think is reasonable for uBooNE and
    // ArgoNeuT 
    float angfactor = 1;
    if(clpar[1] != 0) angfactor = 2 - 1/(1 + fabs(clpar[1]));

    // load the hits starting at the End of the fcl2hits vector.
    // These are the most upstream hits
    unsigned short iht = 0;

    bool first = true;
    unsigned short wire0 = 0;
    for(std::vector<unsigned short>::reverse_iterator it = fcl2hits.rbegin();
       it != fcl2hits.rend(); ++it) {
      unsigned short ihit = *it;
      unsigned short wire = plnhits[ihit]->WireID().Wire;
      if(first) {
        wire0 = wire;
        first = false;
  if(prt) std::cout<<"cl2Fit W:H ";
      }
  if(prt) std::cout<<wire<<":"<<ihit<<" ";
      xwir.push_back(wire - wire0);
      ytim.push_back(plnhits[ihit]->PeakTime());
      ytimerr2.push_back(angfactor * fabs(hiterr2[ihit]));
      if(iht == nht) break;
      ++iht;
    }
    
    nht = iht;
  if(prt) std::cout<<std::endl;
    
    if(nht < 2) return;

    float intcpt = 0.;
    float slope = 0.;
    float intcpterr = 0.;
    float slopeerr = 0.;
    float chidof = 0.;
    LinFit(xwir, ytim, ytimerr2, intcpt, slope, intcpterr, slopeerr, chidof);
    clChisq = chidof;
    if(chidof > 99.) return;
    clpar[0] = intcpt;
    clpar[1] = slope;
    clparerr[0] = intcpterr;
    clparerr[1] = slopeerr;


  if(prt) {
    std::cout<<"nht "<<nht<<" fit par "<<(int)clpar[0]<<" "<<clpar[1]<<" clChisq "<<clChisq;
    std::cout<<std::endl;
  }
  
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2FitChg(const art::PtrVector<recob::Hit>& plnhits)
  {
    // Fits the charge of hits on the fcl2hits vector to a line, or simply
    // uses the average of 1 or 2 hits as determined by NHitsAve

    unsigned short ih0 = fcl2hits.size() - 1;
    
    // don't find the average charge --> no charge cut is made
    if(fNHitsAve[pass] < 1) return;
    
    if(fNHitsAve[pass] < 2) {
      // simply use the charge and width the last hit
      fAveChg = plnhits[fcl2hits[ih0]]->Charge();
    } else if(fNHitsAve[pass] == 2) {
      // average the last two points if requested
      fAveChg = (plnhits[fcl2hits[ih0]]->Charge() + 
                 plnhits[fcl2hits[ih0 - 1]]->Charge()) / 2.;
    } else if((unsigned short)fcl2hits.size() > fNHitsAve[pass]){
      // do a real fit
      std::vector<float> xwir;
      std::vector<float> ychg;
      std::vector<float> ychgerr2;
      unsigned short wire0 = 0;
      bool first = true;
      unsigned int npt = 0;
      // this loop intentionally ignores the Begin hit
      for(unsigned int ii = fcl2hits.size() - 1; ii > 0; --ii) {
        unsigned short wire = plnhits[fcl2hits[ii]]->WireID().Wire;
        if(first) {
          wire0 = wire;
          first = false;
        }
        xwir.push_back((float)(wire - wire0));
        float chg = plnhits[fcl2hits[ii]]->Charge();
        ychg.push_back(chg);
        ychgerr2.push_back(0.5 * chg);
        if(npt == fNHitsAve[pass]) break;
        ++npt;
      }
      if(ychg.size() < 3) return;
      float intcpt; float slope; float intcpterr;
      float slopeerr; float chidof;
      LinFit(xwir, ychg, ychgerr2, intcpt, slope, intcpterr, slopeerr, chidof);
  if(prt) std::cout<<"cl2FitChg Wire: "<<wire0<<" chidof "<<chidof<<" nht "<<ychg.size();
      if(chidof < 20.) fAveChg = intcpt;
    }
  if(prt) std::cout<<" fAveChg "<<(int)fAveChg<<std::endl;

  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2AddHit(const art::PtrVector<recob::Hit>& plnhits,
        unsigned short kwire, bool& HitOK, bool& SigOK)
  {
    // Add a hit to the cluster if it meets several criteria:
    // similar pulse height to the cluster (if fAveChg is defined)
    // closest hit to the project cluster position.
    // Return SigOK if there is a nearby hit that was missed due to the cuts
    
    if(kwire < fFirstWire || kwire > fLastWire) {
      SigOK = false;
      HitOK = false;
      return;
    }
    unsigned short index = kwire - fFirstWire;
    // skip bad wire, but assume the track was there
    if(WireHitRange[index].first == -1) {
      SigOK = true;
      HitOK = false;
      return;
    }
    // return if no signal and no hit
    if(WireHitRange[index].first == -2) {
      SigOK = false;
      HitOK = false;
      return;
    }
    unsigned short firsthit = WireHitRange[index].first;
    unsigned short lasthit = WireHitRange[index].second;

    // Determine if the last hit added was a large (low) charge hit
    // This will be used to prevent adding large (low) charge hits on two
    // consecutive fits
    std::vector<unsigned short>::reverse_iterator it1 = fcl2hits.rbegin();
    unsigned short ih1 = *it1;
    unsigned short wire0 = plnhits[ih1]->WireID().Wire;
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
    for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
  if(prt) {
    std::cout<<"cl2AddHit: Check W:H "<<kwire<<":"<<khit<<" time "<<(int)plnhits[khit]->PeakTime();
    std::cout<<" prtime "<<(short)prtime<<" fAveChg "<<(int)fAveChg;
    std::cout<<" lasthitlow "<<lasthitlow<<" lasthitbig "<<lasthitbig<<std::endl;
  }
      if(khit > plnhits.size()-1) {
        mf::LogError("ClusterCrawler")<<"AddHit bad khit "<<khit;
        return;
      }
      if(hiterr2[khit] < 0) continue;

      float timediff = (plnhits[khit]->PeakTime() - prtime);

      // make hit charge and width cuts
      // don't make a charge ratio cut until fAveChg is defined
      if(fAveChg > 0.) {
        float chgrat = (plnhits[khit]->Charge() - fAveChg) / fAveChg;
  if(prt) std::cout<<" Chg "<<(int)plnhits[khit]->Charge()<<" chgrat "<<chgrat<<std::endl;
        if(prtime < plnhits[khit]->PeakTime() + hitwid[khit] && 
           prtime > plnhits[khit]->PeakTime() - hitwid[khit]) SigOK = true;
        if(lasthitlow) {
          // last hit added was low and this one is as well. If this is a hit
          // that will likely be sele
          if(chgrat < -fChgCut[pass]) {
            // set SigOK false if this low charge hit will be selected to stop
            // following
            if(fabs(timediff) < 3) {
  if(prt) std::cout<<"cl2addhit: found two low charge hits. Pop_back "<<std::endl;
              SigOK = false;
              HitOK = false;
              // lop off the last low charge hit and re-fit?
              if(fcl2hits.size() > 3) {
                fcl2hits.pop_back();
                cl2Fit(plnhits);
              }
              // stop looking for hits
              return;
            }
          }
        } else {
          // allow a low hit
          if(chgrat < lowchgcut) continue;
        }
        // make the high charge ratio cut. 
        // decide on a high charge cut
        if(lasthitbig) {
          // last hit added was big so don't allow another one.
          if(chgrat > fChgCut[pass]) continue;
        } else {
          // allow a big hit
          if(chgrat > bigchgcut) continue;
        }
      } // fAveChg > 0
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
    std::cout<<" clerr "<<prtimerr2<<" hiterr "<<hiterr2[imbest];
    std::cout<<" best "<<best<<std::endl;
  }
    // apply an angle dependent scale factor to the hit error
    float angfactor = 1;
    if(clpar[1] != 0) angfactor = 2 - 1/(1 + fabs(clpar[1]));
    float err = sqrt(prtimerr2 + angfactor * fabs(hiterr2[imbest]));
    // (number of sigma)^2 difference
    float numsig2 = best / err;
    if(numsig2 < 10.) {
      fcl2hits.push_back(imbest);
      if(prt) {
        std::cout<<" >>ADD W:H "<<kwire<<":"<<imbest<<" best "<<best;
        std::cout<<" numsig2 "<<numsig2<<std::endl;
      }
      HitOK = true;
      // decide whether to define/update fAveChg. Only do this if the last hit
      // charge and the added hit charge is not too high and not too low
      if(fAveChg < 0.) {
        cl2FitChg(plnhits);
      } else if(!lasthitlow && !lasthitbig) {
        // note that we take the absolute value here
        float chgrat = fabs(plnhits[imbest]->Charge() - fAveChg) / fAveChg;
        if(chgrat < fChgCut[pass]) cl2FitChg(plnhits);
      }
    } else {
      if(prt) std::cout<<" >>Bad chisq "<<numsig2<<std::endl;
      HitOK = false;
    }
  }

/////////////////////////////////////////
    void ClusterCrawlerAlg::LinFit(std::vector<float>& x, std::vector<float>& y, 
      std::vector<float>& ey2, float& Intercept, float& Slope, 
      float& InterceptError, float& SlopeError, float& ChiDOF) 
    {
      // fit a line ala Bevington linfit.F. The number of points fit is defined by
      // the size of the y vector. 

      ChiDOF = 999.;

      if(y.size() < 2) return;
      if(x.size() < y.size() || ey2.size() < y.size()) return;
      
      float sum = 0.;
      float sumx = 0.;
      float sumy = 0.;
      float sumxy = 0.;
      float sumx2 = 0.;
      float sumy2 = 0.;

      for(unsigned short ii = 0; ii < y.size(); ++ii) {
        float weight = 1. / ey2[ii];
        sum += weight;
        sumx += weight * x[ii];
        sumy += weight * y[ii];
        sumx2 += weight * x[ii] * x[ii];
        sumxy += weight * x[ii] * y[ii];
        sumy2 += weight * y[ii] * y[ii];
      }
      // calculate coefficients and std dev
      float delta = sum * sumx2 - sumx * sumx;
      if(delta == 0.) return;
      float A = (sumx2 * sumy - sumx * sumxy) / delta;
      float B = (sumxy * sum  - sumx * sumy) / delta;
      Intercept = A;
      Slope = B;
      if(x.size() == 2) {
        ChiDOF = 0.;
        return;
      }
      float ndof = x.size() - 2;
      float varnce = (sumy2 + A*A*sum + B*B*sumx2 - 
                      2 * (A*sumy + B*sumxy -A*B*sumx)) / ndof;
      if(varnce > 0.) {
        InterceptError = sqrt(varnce * sumx2 / delta);
        SlopeError = sqrt(varnce * sum / delta);
      } else {
        InterceptError = 0.;
        SlopeError = 0.;
      }
      sum = 0.;
      // calculate chisq
      for(unsigned short ii = 0; ii < y.size(); ++ii) {
        float arg = y[ii] - A - B * x[ii];
        sum += arg * arg / ey2[ii];
      }
      ChiDOF = sum / ndof;
    }

} // namespace cluster
