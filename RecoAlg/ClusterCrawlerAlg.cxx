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
#include "RecoAlg/ClusterCrawlerAlg.h"
#include "RecoAlg/CCHitFinderAlg.h"

// ROOT Includes 
#include "TGraph.h"
// #include "TMath.h"
#include "TF1.h"

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
    fDoVertex           = pset.get< std::vector<bool>  >("DoVertex");
    fLAClusterFix       = pset.get<             bool   >("LAClusterFix");

    fHitErrFac          = pset.get<             float  >("HitErrFac");
    fHitWidFac          = pset.get<             float  >("HitWidFac");
    fDebugPlane         = pset.get<             short  >("DebugPlane");
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
  typedef std::pair<unsigned int, unsigned int> mypair;
  bool SortByLen(const mypair& L, const mypair& R) {return (L.first > R.first);}

  void ClusterCrawlerAlg::cl2Init() {
    prt = false; vtxprt = false;
    NClusters = 0;  clBeginSlp = 0; clBeginSlpErr = 0; clBeginTim = 0;
    clBeginWir = 0; clBeginChg = 0; clEndSlp = 0;      clEndSlpErr = 0;
    clEndTim = 0;   clEndWir = 0;   clEndChg = 0;      clChisq = 0;
    clStopCode = 0; clProcCode = 0; clAssn = 0;        fFirstWire = 0;
    fLastWire = 0; fAveChg = 0; pass = 0; fScaleF = 0;
    tcl.clear(); vtx.clear(); WireHitRange.clear(); hiterr2.clear();
    hitwid.clear();
  }


  void ClusterCrawlerAlg::RunCrawler(std::vector<CCHitFinderAlg::CCHit>& allhits)
  {
    // Run the ClusterCrawler algorithm - creating seed clusters and following
    // them upstream.

    cl2Init();
    
    if(allhits.size() < 3) return;

    
    for(cstat = 0; cstat < geom->Ncryostats(); ++cstat){
      for(tpc = 0; tpc < geom->Cryostat(cstat).NTPC(); ++tpc){
        for(plane = 0; plane < geom->Cryostat(cstat).TPC(tpc).Nplanes(); ++plane){
          WireHitRange.clear();
          hitwid.clear();
          hiterr2.clear();
          // define a code to ensure clusters are compared within the same plane
          clCTP = 100 * cstat + 10 * tpc + plane;
          // fill the WireHitRange vector with first/last hit on each wire
          // dead wires and wires with no hits are flagged < 0
          GetHitRange(allhits, clCTP, WireHitRange, fFirstWire, fLastWire);
          fFirstHit = WireHitRange[0].first;
          unsigned short lasthit = WireHitRange[fLastWire - fFirstWire].second;
          // fill the hiterr2 and hitwid vectors
          for(unsigned short hit = fFirstHit; hit <= lasthit; ++hit) {
            CCHitFinderAlg::CCHit& theHit = allhits[hit];
            // hit time error
            float arg = fHitErrFac * theHit.RMS;
            hiterr2.push_back(arg * arg);
            // hit width for checking that the Signal is OK on a wire
            hitwid.push_back(fHitWidFac * theHit.RMS);
          }
          // get the scale factor to convert dTick/dWire to dX/dU. This is used
          // to make the kink and merging cuts
          art::Ptr<recob::Wire> theWire = allhits[fFirstHit].Wire;
          uint32_t channel = theWire->RawDigit()->Channel();
          float wirePitch = geom->WirePitch(geom->View(channel));
          float tickToDist = larprop->DriftVelocity(larprop->Efield(),larprop->Temperature());
          tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns
          fScaleF = tickToDist / wirePitch;
      
          fMaxTime = detprop->NumberTimeSamples();
          fNumWires = geom->Nwires(plane);
          
          // look for clusters
          cl2ClusterLoop(allhits);
        } // plane
      } // tpc
    } // cstat
    
    WireHitRange.clear(); 
    hiterr2.clear();
    hitwid.clear();
    
  } // RunCrawler
    
////////////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ClusterLoop(std::vector<CCHitFinderAlg::CCHit>& allhits) {
      // looks for seed clusters in a plane and follows them

      unsigned short nHitsUsed = 0;
      bool AllDone = false;
      for(unsigned short thispass = 0; thispass < fNumPass; ++thispass) {
        pass = thispass;
/*
  if(fDebugPlane == (short)plane) {
    std::cout<<"******** ClusterCrawler plane "<<plane;
    std::cout<<" pass "<<pass<<" ****************"<<std::endl;
  }
*/
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
            if(ihit > allhits.size()-1) {
              mf::LogError("ClusterCrawler")<<"RunCrawler bad ihit "<<ihit;
              return;
            }
            if(hiterr2[ihit - fFirstHit] < 0) continue;
            // skip deleted hits
            if(allhits[ihit].Charge < 0) continue;
            // skip multiple hits except on the last pass
            if(pass < fNumPass - 1 && allhits[ihit].numHits > 1) continue;
            if((iwire - span + 1) < fFirstWire) continue;
            for(unsigned short jwire = iwire - span + 1; jwire < iwire; ++jwire) {
              unsigned short jindx = jwire - fFirstWire;
              if(WireHitRange[jindx].first < 0) continue;
              // Find the hit on wire jwire that best matches a line between
              // a nearby vertex and hit ihit. No constraint if useHit < 0
              unsigned short useHit = 0;
              bool doConstrain = false;
              cl2VtxConstraint(allhits, vtx, iwire, ihit, jwire, useHit, doConstrain);
              unsigned short jfirsthit = WireHitRange[jindx].first;
              unsigned short jlasthit = WireHitRange[jindx].second;
              for(unsigned short jhit = jfirsthit; jhit < jlasthit; ++jhit) {
                if(jhit > allhits.size()-1) {
                  mf::LogError("ClusterCrawler")<<"RunCrawler bad jhit "<<jhit;
                  return;
                }
                if(hiterr2[jhit - fFirstHit] < 0) continue;
                // skip deleted hits
                if(allhits[jhit].Charge < 0) continue;
                // Vertex constraint
                if(doConstrain && jhit != useHit) continue;
                // start a cluster with these two hits
                fcl2hits.clear();
                fAveChg = -1.;
                clBeginChg = -1.;
                clStopCode = 0;
                clProcCode = pass;
                clAssn = -1; 
                fcl2hits.push_back(ihit);
                fcl2hits.push_back(jhit);
                clpar[0] = allhits[jhit].Time;
                clpar[1] = (allhits[ihit].Time - allhits[jhit].Time) / (iwire - jwire);
                clChisq = 0;
                // now look for hits to add on the intervening wires
                bool SigOK = false;
                bool HitOK = false;
                bool clok = true;
                for(unsigned short kwire = jwire+1; kwire < iwire; ++kwire) {
                  cl2AddHit(allhits, kwire, HitOK, SigOK);
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
                // hit number. This only needs to be done on the starting cluster.
                // Hits will be added in the proper order by cl2Follow
                std::sort(fcl2hits.begin(), fcl2hits.end(), SortByLowHit);
                // do a real fit
                cl2Fit(allhits);
                if(clChisq > 10.) continue;
                // check the charge ratio between the DS hit and the next-most
                // DS hit. This ratio is < 2 for a stopping particle. A large
                // ratio indicates that we are starting a cluster just US of a
                // high ionization region
                float chg0 = allhits[fcl2hits[0]].Charge;
                float chg1 = allhits[fcl2hits[1]].Charge;
                if(chg0 > 2 * chg1) continue;
                // save the cluster begin info
                clBeginWir = iwire;
                clBeginTim = allhits[ihit].Time;
                clBeginSlp = clpar[1];
                clBeginSlpErr = clparerr[1];
                // follow a trail of hits upstream
                cl2FollowUS(allhits);
                // do a quality check. fcl2hits size set 0 if bad cluster
                cl2QACheck(allhits);
                // handle large angle cluster hits?
                if(fLAClusterFix) cl2LAClusterFix(allhits, tcl);
                if(fcl2hits.size() >= fMinHits[pass]) {
                  // it's long enough so save it
                  clEndSlp = clpar[1]; // save the slope at the end
                  clEndSlpErr = clparerr[1];
                  clEndChg = fAveChg;
                  cl2TmpStore(allhits, tcl); // store the cluster
                  ClusterAdded = true;
                  nHitsUsed += fcl2hits.size();
                  AllDone = (nHitsUsed == allhits.size());
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
                    cl2TmpStore(allhits, tcl);
                    ClusterAdded = true;
                    nHitsUsed += fcl2hits.size();
                    AllDone = (nHitsUsed == allhits.size());
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
        if(fDoMerge[pass]) cl2ChkMerge(allhits, tcl);
        // form 2D vertices
        if(fDoVertex[pass]) cl2DoVertex(allhits, tcl, vtx);
        if(AllDone) break;
      } // pass

      // split clusters using vertices
      cl2VtxClusterSplit(allhits, tcl, vtx);

  if(fDebugPlane == (short)plane) {
    std::cout<<"Clustering done in plane = "<<plane<<std::endl;
    cl2Print(allhits, tcl);
  }

    
  } // cl2ClusterLoop

//////////////////////////////////////////
    void ClusterCrawlerAlg::cl2QACheck(
      std::vector<CCHitFinderAlg::CCHit>& allhits)
    {
      // ignore short large angle clusters that have multiplicity 1 hits
      if(fcl2hits.size() < 5 && fabs(clpar[1]) > 3.) {
        for(unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
          if(allhits[fcl2hits[ii]].numHits == 1) {
            fcl2hits.clear();
            return;
          }
        } // ii
      } // fcl2hits.size() < 5
    } // cl2QACheck

//////////////////////////////////////////
    void ClusterCrawlerAlg::cl2VtxConstraint(
        std::vector<CCHitFinderAlg::CCHit>& allhits,
        std::vector<VtxStore>& vtx,
        unsigned short iwire, unsigned short ihit, unsigned short jwire,
        unsigned short& useHit, bool& doConstrain)
    {
      // checks hits on wire jwire to see if one is on a line between a US vertex
      // and the hit ihit on wire iwire. If one is found, doConstrain is set true
      // and the hit index is returned
      doConstrain = false;
      if(vtx.size() == 0) return;
      // no vertices made yet on the first pass
      if(pass == 0) return;
      // skip if vertices were not requested to be made on the previous pass
      if( !fDoVertex[pass - 1] ) return;
      
      unsigned short jindx = jwire - fFirstWire;
      unsigned short jfirsthit = WireHitRange[jindx].first;
      unsigned short jlasthit = WireHitRange[jindx].second;
      for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
        if(vtx[iv].CTP != clCTP) continue;
        // vertex must be US of the cluster
        if(vtx[iv].Wire > jwire) continue;
        // but not too far US
        if(vtx[iv].Wire < jwire - 10) continue;
        clpar[0] = allhits[ihit].Time;
        clpar[1] = (vtx[iv].Time - allhits[ihit].Time) / (vtx[iv].Wire - iwire);
        float prtime = clpar[0] + clpar[1] * (jwire - iwire);
        for(unsigned short jhit = jfirsthit; jhit < jlasthit; ++jhit) {
          if(hiterr2[jhit - fFirstHit] < 0) continue;
          if(allhits[jhit].Charge < 0) continue;
          float tdiff = fabs(prtime - allhits[jhit].Time) / allhits[jhit].RMS;
          if(tdiff < 2.5) {
            useHit = jhit;
            doConstrain = true;
            return;
          }
        } // jhit
      } // iv
    } // cl2VtxConstraint


/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2VtxClusterSplit(
        std::vector<CCHitFinderAlg::CCHit>& allhits,
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx)
    {

      // split clusters that cross vertices

      if(vtx.size() == 0) return;
      
      for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        // find max and min times to make rough cuts
        float tmax = tcl[icl].BeginTim;
        if(tcl[icl].EndTim > tmax) tmax = tcl[icl].EndTim;
        float tmin = tcl[icl].BeginTim;
        if(tcl[icl].EndTim < tmin) tmin = tcl[icl].EndTim;
        // pad a bit to allow for cluster wandering
        tmax += 20.;
        tmin -= 20.;
        bool didSplit = false;
        for(unsigned short ivx = 0; ivx < vtx.size(); ++ivx) {
          if(vtx[ivx].CTP != clCTP) continue;
          // ignore vertices near the ends of clusters
          if(vtx[ivx].Wire < tcl[icl].EndWir + 3) continue;
          if(vtx[ivx].Wire > tcl[icl].BeginWir - 3) continue;
          if(vtx[ivx].Time > tmax) continue;
          if(vtx[ivx].Time < tmin) continue;
          for(unsigned short ii = 3; ii < tcl[icl].tclhits.size() - 3; ++ii) {
            unsigned short iht = tcl[icl].tclhits[ii];
            if(vtx[ivx].Wire == allhits[iht].WireNum) {
              if(fabs(vtx[ivx].Time - allhits[iht].Time) < 10.) {
                cl2DoSplit(allhits, tcl, icl, ii, ivx);
                didSplit = true;
                break;
              } //
            } // vtx[iv].Wire == allhits[iht].WireNum
            if(didSplit) break;
          } // iht
          if(didSplit) break;
        } // ivx
      } // icl

      // trim hits from the ends of clusters if they cross a vertex wire
      for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].CTP != clCTP) continue;
        if(tcl[icl].EndVtx >= 0) {
          unsigned short evx = tcl[icl].EndVtx;
          if(tcl[icl].EndWir < vtx[evx].Wire) {
//            std::cout<<"Trim hits? "<<tcl[icl].ID<<" EndWir "<<tcl[icl].EndWir;
//            std::cout<<" Vtx wire "<<vtx[evx].Wire<<std::endl;
          }
        } // tcl[icl].EndVtx >= 0
      } // icl
      

    } // cl2VtxClusterSplit


//////////////////////////////////////////
    void ClusterCrawlerAlg::cl2LAClusterFix(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl)
    {
      // Multiple hits on a single wire associated with a large angle cluster
      // are merged into one hit
      
      // ensure that this is a cluster that will be stored
      bool keeper = false;
      if(fcl2hits.size() >= fMinHits[pass]) keeper = true;
      if(pass < fNumPass - 1) {
        if(fcl2hits.size() >= fMinHits[pass+1]) keeper = true;
      }
      if(!keeper) return;
      
  prt = false;
/*
  for(unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
    unsigned short iht = fcl2hits[iht];
    if(plane == 0 && allhits[iht].WireNum == 1637 && iht == 440) {
      std::cout<<"Found your hit. Mult = "<<allhits[iht].numHits<<std::endl;
      prt = true;
      break;
    }
  }
  if(prt) {
    std::cout<<"LAClusterFix cluster W:H "<<allhits[fcl2hits[0]].WireNum;
    std::cout<<":"<<fcl2hits[0];
  }
*/      
      // require a large angle cluster
      float theta = atan(fScaleF * clBeginSlp);
  if(prt) std::cout<<" Theta "<<theta<<std::endl;
      if(fabs(theta) < 1) return;
      
      // merge the Multiplicity > 1 hits in this cluster into one hit / wire
      for(unsigned short ii = 0; ii < fcl2hits.size(); ++ii) {
        unsigned short iht = fcl2hits[ii];
        if(allhits[iht].numHits > 1) {
          // index of the Hit that we are going to merge into
          unsigned short theHit = 0;
          // need to find the range of time ticks over which to sum the charge
          // and determine the RMS
          short loTime = 9999;
          short hiTime = 0;
          unsigned short nGaus = 0;
          float sumchg = 0.;
          for(unsigned short jj = 0; jj < allhits[iht].numHits; ++jj) {
            unsigned short jht = allhits[iht].LoHitID + jj;
  if(prt) {
    std::cout<<" W:H "<<allhits[jht].WireNum<<":";
    std::cout<<jht<<" time "<<allhits[jht].Time<<" chg "<<(int)allhits[jht].Charge<<std::endl;
  }
            // error checking
            if(allhits[jht].LoHitID != allhits[iht].LoHitID) {
              std::cout<<"Hit multiplet ID error "<<jj<<" "<<allhits[iht].numHits<<std::endl;
              return;
            }
            // hit is not used by another cluster
            if(hiterr2[jht - fFirstHit] < 0.) continue;
            if(theHit == 0) theHit = jht;
            short arg = (short)(allhits[jht].Time - 2.5 * allhits[jht].RMS);
            if(arg < loTime) loTime = arg;
            arg = (short)(allhits[jht].Time + 2.5 * allhits[jht].RMS);
            if(arg > hiTime) hiTime = arg;
            sumchg += allhits[jht].Charge;
            ++nGaus;
          } // jj
          // all hits used?
          if(theHit == 0) return;
          // create a TF1 to define the hit shape
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
          TF1 *Gn = new TF1("gn",eqn.c_str());
          // define the params
          unsigned short iGaus = 0;
          for(unsigned short jj = 0; jj < allhits[iht].numHits; ++jj) {
            unsigned short jht = allhits[iht].LoHitID + jj;
            if(hiterr2[jht - fFirstHit] < 0.) continue;
            unsigned short index = iGaus * 3;
            Gn->SetParameter(index    , allhits[jht].Amplitude);
            Gn->SetParameter(index + 1, allhits[jht].Time);
            Gn->SetParameter(index + 2, allhits[jht].RMS);
            ++iGaus;
          }
          // integrate and get the RMS
          float chgsum = 0.;
          float chgsumt = 0.;
          std::vector<float> signal;
          for(short time = loTime; time <= hiTime; ++time) {
            float arg = Gn->Eval((Double_t)time, 0, 0);
            signal.push_back(arg);
            chgsum   += arg;
            chgsumt  += arg * time;
          }
          delete Gn;
          float aveTime = chgsumt / chgsum;
          // find the RMS
          chgsumt = 0.;
          for(short time = loTime; time <= hiTime; ++time) {
            unsigned short index = time - loTime;
  if(index > signal.size() - 1) {
    std::cout<<"LAClusterFix Bad index "<<std::endl;
    return;
  }
            short dtime = time - (short)aveTime;
            chgsumt += signal[index] * dtime * dtime;
          }
          allhits[theHit].RMS = sqrt(chgsumt / chgsum);
          allhits[theHit].Time = aveTime;
          allhits[theHit].Charge = chgsum;
          // find the amplitude from the integrated charge and the RMS
          allhits[theHit].Amplitude = chgsum / (2.507 * allhits[theHit].RMS);
          // associate theHit with the current cluster
          fcl2hits[ii] = theHit;
  if(prt) {
    std::cout<<"theHit "<<theHit<<" time "<<(int)aveTime<<" RMS "<<allhits[theHit].RMS;
    std::cout<<" chg "<<(int)chgsum<<" Amp "<<allhits[theHit].Amplitude<<std::endl;
  }
          // set the charge of the rest < 0
          for(unsigned short jj = 0; jj < allhits[iht].numHits; ++jj) {
            unsigned short jht = allhits[iht].LoHitID + jj;
            if(hiterr2[jht - fFirstHit] < 0.) continue;
            if(jht == theHit) continue;
            allhits[jht].Charge = -1.;
            hiterr2[jht - fFirstHit] = -1.;
  if(prt) std::cout<<"clobber hit "<<jht<<std::endl;
          } // jj
        } // allhits[iht].numHits > 1
      } // ii

      clProcCode += 200;
      
    } // cl2LAClusterFix

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2DoVertex(std::vector<CCHitFinderAlg::CCHit>& allhits,
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx)
    {
      // try to make 2D vertices
      
      if(tcl.size() < 2) return;

      // form vertices starting with the longest
      std::map<unsigned short, unsigned short> sortindex;
      cl2SortByLength(tcl,sortindex);
      
      float nwires = fNumWires;
      float maxtime = fMaxTime;

  vtxprt = (fDebugPlane == (short)plane && fDebugHit < 0);
  if(vtxprt) {
    std::cout<<"DoVertex plane "<<plane<<" pass "<<pass<<std::endl;
    cl2Print(allhits,tcl);
  }

      for(unsigned short ii1 = 0; ii1 < sortindex.size() - 1; ++ii1) {
        unsigned short it1 = sortindex[ii1];
        // ignore obsolete clusters
        if(tcl[it1].ID < 0) continue;
        // ignore already attached clusters
        if(tcl[it1].BeginVtx >= 0 && tcl[it1].EndVtx >= 0) continue;
        float es1 = tcl[it1].EndSlp;
        float eth1 = atan(fScaleF * tcl[it1].EndSlp);
        unsigned short ew1 = tcl[it1].EndWir;
        float et1 = tcl[it1].EndTim;
        float bs1 = tcl[it1].BeginSlp;
        float bth1 = atan(fScaleF * tcl[it1].BeginSlp);
        unsigned short bw1 = tcl[it1].BeginWir;
        float bt1 = tcl[it1].BeginTim;
        for(unsigned short ii2 = ii1 + 1; ii2 < sortindex.size(); ++ii2) {
          unsigned short it2 = sortindex[ii2];
          if(tcl[it2].ID < 0) continue;
          // ignore already attached clusters
          if(tcl[it2].BeginVtx >= 0 && tcl[it2].EndVtx >= 0) continue;
          // try to attach cluster to existing vertices at either end
          cl2ClsVertex(allhits, tcl, vtx, it2);
          // ignore if both clusters are short
          if(tcl[it1].tclhits.size() < 10 &&
             tcl[it2].tclhits.size() < 10) continue;
          float es2 = tcl[it2].EndSlp;
          float eth2 = atan(fScaleF * tcl[it2].EndSlp);
          unsigned short ew2 = tcl[it2].EndWir;
          float et2 = tcl[it2].EndTim;
          float bs2 = tcl[it2].BeginSlp;
          float bth2 = atan(fScaleF * tcl[it2].BeginSlp);
          unsigned short bw2 = tcl[it2].BeginWir;
          float bt2 = tcl[it2].BeginTim;
//  if(vtxprt) std::cout<<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID<<std::endl;
  // topo 1: check for vtx US of the ends of both clusters
          float dth = fabs(eth1 - eth2);
          if(tcl[it1].EndVtx < 0 && tcl[it2].EndVtx < 0 && dth > 0.3) {
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
  if(vtxprt) {
    std::cout<<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID;
    std::cout<<" topo1 vtx wire "<<vw<<" time "<<(int)fvt<<std::endl;
  }
                if(fvt > 0. && fvt < maxtime) {
                  // vertex wire US of cluster ends and time in the detector
                  // Check this against existing vertices and update
                  cl2ChkVertex(allhits, tcl, vtx, vw, fvt, it1, it2, 1);
                } // fvt in detector
              } // vw topo 1 check
            } // fvw in detector
          } // topo 1
  // topo 2: check for vtx US of it1 and DS of it2
          dth = fabs(eth1 - bth2);
          if(tcl[it1].EndVtx < 0 && tcl[it2].BeginVtx < 0 && dth > 0.3) {
            float dsl = bs2 - es1;
            if(fabs(dsl) < 0.001) dsl = 0.001;
            float fvw = 0.5 + (et1 - ew1 * es1 - bt2 + bw2 * bs2) / dsl;
            if(fvw > 0 && fvw < nwires) {
              // vertex wire in the detector
              unsigned short vw = fvw;
              if(vw <= ew1 && vw >= bw2) {
                float fvt = et1 + (vw - ew1) * es1;
  if(vtxprt) {
    std::cout<<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID;
    std::cout<<" topo2 vtx wire "<<vw<<" time "<<(int)fvt<<std::endl;
  }
                if(fvt > 0. && fvt < maxtime) {
                  cl2ChkVertex(allhits, tcl, vtx, vw, fvt, it1, it2, 2);
                } // fvt in detector
              } // vw topo 2 check
            } // fvw in detector
          } // topo 2
  // topo 3: check for vtx DS of it1 and US of it2
          dth = fabs(bth1 - eth2);
          if(tcl[it1].BeginVtx < 0 && tcl[it2].EndVtx < 0 && dth > 0.3) {
            float dsl = bs1 - es2;
            if(fabs(dsl) < 0.001) dsl = 0.001;
            float fvw = 0.5 + (et2 - ew2 * es2 - bt1 + bw1 * bs1) / dsl;
            if(fvw > 0 && fvw < nwires) {
              unsigned short vw = fvw;
              if(vw <= ew2 && vw >= bw1) {
                float fvt = et2 + (vw - ew2) * es2;
  if(vtxprt) {
    std::cout<<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID;
    std::cout<<" topo3 vtx wire "<<vw<<" time "<<(int)fvt<<std::endl;
  }
                if(fvt > 0. && fvt < maxtime) {
                  cl2ChkVertex(allhits, tcl, vtx, vw, fvt, it1, it2, 3);
                } // fvt in detector
              } // vw topo 3 check
            } // fvw in detector
          } // topo 3
  // topo 4: check for vtx DS of it1 and DS of it2
          dth = fabs(bth1 - bth2);
          if(tcl[it1].BeginVtx < 0 && tcl[it2].BeginVtx < 0 && dth > 0.3) {
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
  if(vtxprt) {
    std::cout<<"Chk clusters "<<tcl[it1].ID<<" "<<tcl[it2].ID;
    std::cout<<" topo4 vtx wire "<<vw<<" time "<<(int)fvt<<std::endl;
  }
                if(fvt > 0. && fvt < maxtime) {
                  // vertex wire US of cluster ends and time in the detector
                  // Check this against existing vertices and update
                  cl2ChkVertex(allhits, tcl, vtx, vw, fvt, it1, it2, 4);
                } // fvt in detector
              } // vw topo 4 check
            } // fvw in detector
          } // it2
        } // topo 4
      } // it1
      

      // fit the vertices
      for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
        if(vtx[iv].CTP != clCTP) continue;
        if(vtx[iv].Wght < 0) continue;
        float ChiDOF = 0.;
        cl2VtxFit(tcl, vtx, iv, ChiDOF);
      }

      // set the vertex weights
      for(unsigned short it = 0; it < tcl.size(); ++it) {
        if(tcl[it].ID < 0) continue;
        if(tcl[it].CTP != clCTP) continue;
        // cluster weight = number of hits, saturated at 10
        float cw = tcl[it].tclhits.size();
        if(cw > 10) cw = 10;
        if(tcl[it].BeginVtx >=0) vtx[tcl[it].BeginVtx].Wght += cw;
        if(tcl[it].EndVtx >=0) vtx[tcl[it].EndVtx].Wght += cw;
      }
      

  if(vtxprt) {
    std::cout<<"Vertices "<<vtx.size()<<std::endl;
    for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
      if(vtx[iv].CTP != clCTP) continue;
      std::cout<<"vtx "<<iv<<" wire "<<vtx[iv].Wire<<" time "<<(int)vtx[iv].Time<<" wght "<<(int)vtx[iv].Wght;
      std::cout<<" topo "<<vtx[iv].Topo<<std::endl;
    }
    cl2Print(allhits, tcl);
  }
    }

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ClsVertex(std::vector<CCHitFinderAlg::CCHit>& allhits, 
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
        unsigned short it)
    {
      // try to attach cluster it to an existing vertex
      
      if(vtx.size() == 0) return;
      
      for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
        // ignore vertices in the wrong cryostat/TPC/Plane
        if(vtx[iv].CTP != clCTP) continue;
        // ignore deleted vertices
        if(vtx[iv].Wght < 0) continue;
        // determine which end to match - begin or end
        if(tcl[it].EndVtx < 0 &&  vtx[iv].Wire <= tcl[it].EndWir + 2) {
          // project cluster to US vertex
          float tdiff = tcl[it].EndTim + (vtx[iv].Wire - tcl[it].EndWir) * tcl[it].EndSlp - vtx[iv].Time;
          if(fabs(tdiff) < 10) {
            bool SigOK = false;
            cl2ChkSignal(allhits, vtx[iv].Wire, vtx[iv].Time, tcl[it].EndWir, tcl[it].EndTim, SigOK);
            if(SigOK) {
              // good match
              tcl[it].EndVtx = iv;
              // re-fit it
              float ChiDOF = 99.;
              cl2VtxFit(tcl, vtx, iv, ChiDOF);
  if(vtxprt) std::cout<<"Add End "<<it<<" to vtx "<<iv<<" chi "<<ChiDOF<<std::endl;
              return;
            } // SigOK
          } // tdiff
        } else if(tcl[it].BeginVtx < 0 &&  vtx[iv].Wire >= tcl[it].BeginWir - 2) {
          // project cluster to DS vertex
          float tdiff = tcl[it].BeginTim + (vtx[iv].Wire - tcl[it].BeginWir) * tcl[it].BeginSlp - vtx[iv].Time;
          if(fabs(tdiff) < 10) {
            bool SigOK = false;
            cl2ChkSignal(allhits, vtx[iv].Wire, vtx[iv].Time, tcl[it].BeginWir, tcl[it].BeginTim, SigOK);
            if(SigOK) {
              // good match
              tcl[it].BeginVtx = iv;
              // re-fit it
              float ChiDOF = 99.;
              cl2VtxFit(tcl, vtx, iv, ChiDOF);
  if(vtxprt) std::cout<<"Add Begin "<<it<<" to vtx "<<iv<<" chi "<<ChiDOF<<std::endl;
              return;
            } // SigOK
          } // tdiff
        } // vtx[iv].Wire
      } // iv
    }



/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ChkVertex(std::vector<CCHitFinderAlg::CCHit>& allhits,
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
        

  if(vtxprt) {
    std::cout<<"ChkVertex "<<tcl[it1].ID<<" "<<tcl[it2].ID<<" topo "<<topo;
    std::cout<<" vw "<<vw<<" vt "<<(int)fvt<<std::endl;
  }
        // check vertex and clusters for proximity to existing vertices
        bool SigOK = false;
        for(unsigned short iv = 0; iv < vtx.size(); ++iv) {
          if(vtx[iv].CTP != clCTP) continue;
          if( abs( vw - vtx[iv].Wire) < 2 &&
             fabs(fvt - vtx[iv].Time) < 10) {
            // got a match. Check the appropriate cluster end and attach
            if( (topo == 1 || topo == 2) && tcl[it1].EndVtx < 0) {
              cl2ChkSignal(allhits, vw, fvt, tcl[it1].EndWir, tcl[it1].EndTim, SigOK);
              if(SigOK) tcl[it1].EndVtx = iv;
  if(vtxprt)  std::cout<<"12 Attach cl "<<tcl[it1].ID<<" to vtx "<<iv<<std::endl;
            } else if( (topo == 3 || topo == 4) && tcl[it1].BeginVtx < 0) {
              cl2ChkSignal(allhits, vw, fvt, tcl[it1].BeginWir, tcl[it1].BeginTim, SigOK);
              if(SigOK) tcl[it1].BeginVtx = iv;
  if(vtxprt)  std::cout<<"34 Attach cl "<<tcl[it1].ID<<" to vtx "<<iv<<std::endl;
            } // cluster it1
            if( (topo == 1 || topo == 3) && tcl[it2].EndVtx < 0) {
              cl2ChkSignal(allhits, vw, fvt, tcl[it2].EndWir, tcl[it2].EndTim, SigOK);
              if(SigOK) tcl[it2].EndVtx = iv;
  if(vtxprt) std::cout<<"13 Attach cl "<<tcl[it2].ID<<" to vtx "<<iv<<std::endl;
            } else if( (topo == 2 || topo == 4) && tcl[it2].BeginVtx < 0) {
              cl2ChkSignal(allhits, vw, fvt, tcl[it2].BeginWir, tcl[it2].BeginTim, SigOK);
              if(SigOK) tcl[it2].BeginVtx = iv;
  if(vtxprt) std::cout<<"24 Attach cl "<<tcl[it2].ID<<" to vtx "<<iv<<std::endl;
            } // cluster it2
            return;
          } // matched vertex
        } // iv

        // no match to existing vertices. Ensure that there is a wire signal between
        // the vertex and the appropriate ends of the clusters
        bool Sig1OK = false;
        bool Sig2OK = false;
        if(topo == 1 || topo == 2) {
          cl2ChkSignal(allhits, vw, fvt, tcl[it1].EndWir, tcl[it1].EndTim, Sig1OK);
        } else {
          cl2ChkSignal(allhits, vw, fvt, tcl[it1].BeginWir, tcl[it1].BeginTim, Sig1OK);
        }
        if(topo == 1 || topo == 3) {
          cl2ChkSignal(allhits, vw, fvt, tcl[it2].EndWir, tcl[it2].EndTim, Sig2OK);
        } else {
          cl2ChkSignal(allhits, vw, fvt, tcl[it2].BeginWir, tcl[it2].BeginTim, Sig2OK);
        }
  if(vtxprt) {
    std::cout<<"Chk new "<<tcl[it1].ID<<" OK "<<Sig1OK<<" "<<tcl[it2].ID<<" OK "<<Sig2OK;
    std::cout<<" Vtx at "<<vw<<" "<<(int)fvt<<std::endl;
  }
        // both clusters must have an OK signal
        if(Sig1OK && Sig2OK) {
          VtxStore newvx;
          newvx.Wire = vw;
          newvx.Time = fvt;
          newvx.Wght = 0;
          newvx.Topo = topo;
          newvx.CTP = clCTP;
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
  if(vtxprt) {
    std::cout<<"New vtx "<<iv<<" in plane "<<plane<<" topo "<<topo<<" cls "<<tcl[it1].ID<<" "<<tcl[it2].ID<<std::endl;
    std::cout<<" time "<<(int)fvt<<" wire "<<vw<<std::endl;
  }
        }

      }

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ChkSignal(std::vector<CCHitFinderAlg::CCHit>& allhits,
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
      if(wiree < fFirstWire || wiree > fLastWire) return;
      if(wireb < fFirstWire || wireb > fLastWire) return;
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
        unsigned short index = wire - fFirstWire;
        // skip dead wires
        if(WireHitRange[index].first == -1) continue;
        // no hits on this wire
        if(WireHitRange[index].first == -2) {
          SigOK = false;
          return;
        }
        unsigned short firsthit = WireHitRange[index].first;
        unsigned short lasthit = WireHitRange[index].second;
//  std::cout<<"ChkSignal "<<wiree<<" "<<wire<<" "<<wireb<<" "<<(int)prtime;
//  std::cout<<" first "<<firsthit<<" last "<<lasthit;
        for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
          // skip obsolete hits
          if(allhits[khit].Charge < 0) continue;
          float tdiff = fabs(prtime - allhits[khit].Time);
          if (tdiff < 2 * hitwid[khit - fFirstHit]) {
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
    void ClusterCrawlerAlg::cl2DoSplit(std::vector<CCHitFinderAlg::CCHit>& allhits,
        std::vector<ClusterStore>& tcl, unsigned short icl, unsigned short pos, 
        unsigned short ivx)
    {
      // split cluster icl into two clusters starting at hit position pos

      // Create the first cluster using the Begin info
      clBeginSlp = tcl[icl].BeginSlp;
      clBeginSlpErr = tcl[icl].BeginSlpErr;
      clBeginWir = tcl[icl].BeginWir;
      clBeginTim = tcl[icl].BeginTim;
      clBeginChg = tcl[icl].BeginChg;
      clStopCode = 5;
      clProcCode = tcl[icl].ProcCode;
      fcl2hits.clear();
      for(unsigned short ii = 0; ii < pos; ++ii) {
        fcl2hits.push_back(tcl[icl].tclhits[ii]);
      }
      // determine the pass in which this cluster was created
      pass = tcl[icl].ProcCode - 10 * (tcl[icl].ProcCode / 10);
      // fit the end hits
      cl2Fit(allhits);
      clEndSlp = clpar[1];
      clEndSlpErr = clparerr[1];
      // find the charge at the end
      cl2FitChg(allhits);
      clEndChg = fAveChg;
      cl2TmpStore(allhits, tcl);
      // associate the End with the supplied vertex
      unsigned short iclnew = tcl.size() - 1;
      tcl[iclnew].EndVtx = ivx;
      tcl[iclnew].BeginVtx = tcl[icl].BeginVtx;

      // now create the second cluster
      clEndSlp = tcl[icl].EndSlp;
      clEndSlpErr = tcl[icl].EndSlpErr;
      clEndWir = tcl[icl].EndWir;
      clEndTim = tcl[icl].EndTim;
      clEndChg = tcl[icl].EndChg;
      clStopCode = 5;
      clProcCode = tcl[icl].ProcCode;
      fcl2hits.clear();
      for(unsigned short ii = pos; ii < tcl[icl].tclhits.size(); ++ii) {
        fcl2hits.push_back(tcl[icl].tclhits[ii]);
        // define the Begin parameters
        if(fcl2hits.size() == fMaxHitsFit[pass] ||
           fcl2hits.size() == fMinHits[pass]) {
          cl2Fit(allhits);
          clBeginSlp = clpar[1];
          clBeginSlpErr = clparerr[1];
        }
        if((unsigned short)fcl2hits.size() == fNHitsAve[pass] + 1) {
          cl2FitChg(allhits);
          clBeginChg = fAveChg;
        }
      } // ii
      cl2TmpStore(allhits, tcl);
      // associate the End with the supplied vertex
      iclnew = tcl.size() - 1;
      tcl[iclnew].BeginVtx = ivx;
      tcl[iclnew].EndVtx = tcl[icl].EndVtx;
      
      // declare icl obsolete
      tcl[icl].ID = -tcl[icl].ID;
      
    } // cl2DoSplit

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ChkMerge(std::vector<CCHitFinderAlg::CCHit>& allhits,
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
        if(tcl[it1].CTP != clCTP) continue;
        float bs1 = tcl[it1].BeginSlp;
        // convert slope to angle
        float arg = fScaleF * bs1;
        float bth1 = atan(arg);
        // error on the angle
        float bth1e = fScaleF * tcl[it1].BeginSlpErr / (1 + arg * arg);
        // more compact notation for begin/end, wire/time/chg/slp/theta, 1/2
        unsigned short bw1 = tcl[it1].BeginWir;
        float bt1 = tcl[it1].BeginTim;
        float bc1 = tcl[it1].BeginChg;
        float es1 = tcl[it1].EndSlp;
        // convert slope to angle
        arg = fScaleF * es1;
        float eth1 = atan(arg);
        float eth1e =  fScaleF * tcl[it1].EndSlpErr / (1 + arg * arg);
        unsigned short ew1 = tcl[it1].EndWir;
        float et1 = tcl[it1].EndTim;
        float ec1 = tcl[it1].EndChg;
        unsigned short pass1 = tcl[it1].ProcCode - 10 * (tcl[it1].ProcCode / 10);
        for(unsigned short it2 = it1 + 1; it2 < tclsize; ++it2) {
          // ignore already merged clusters
          if(tcl[it1].ID < 0) continue;
          if(tcl[it2].ID < 0) continue;
          // only merge if they are in the right cryostat/TPC/plane
          if(tcl[it2].CTP != clCTP) continue;
          float bs2 = tcl[it2].BeginSlp;
          float bs2e = tcl[it2].BeginSlpErr;
          // convert slope to angle
          arg = fScaleF * bs2;
          float bth2 = atan(arg);
          // error on the angle
          float bth2e = fScaleF * tcl[it2].BeginSlpErr / (1 + arg * arg);
          unsigned short bw2 = tcl[it2].BeginWir;
          float bt2 = tcl[it2].BeginTim;
          float bc2 = tcl[it2].BeginChg;
          float es2 = tcl[it2].EndSlp;
          arg = fScaleF * es2;
          float eth2 = atan(arg);
          float eth2e = fScaleF * tcl[it2].EndSlpErr / (1 + arg * arg);
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
              cl2DoMerge(allhits, tcl, vtx, it2, it1, 10);
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
              cl2DoMerge(allhits, tcl, vtx, it1, it2, 10);
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
                cl2ChkMerge12(allhits, tcl, it1, it2, didit);
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
                cl2ChkMerge12(allhits, tcl, it2, it1, didit);
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
  void ClusterCrawlerAlg::cl2ChkMerge12(std::vector<CCHitFinderAlg::CCHit>& allhits, 
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
      unsigned short wire = allhits[hit].WireNum;
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
      unsigned short wiron2 = allhits[hiton2].WireNum;
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
    if(hiton1 > allhits.size() - 1) {
      mf::LogError("ClusterCrawler")<<"ChkMerge12 bad hiton1 "<<hiton1;
      return;
    }
    // check the time difference
    float timon1 = allhits[hiton1].Time;
    float dtim = fabs(et2 + (wiron1 - ew2) * tcl[it2].EndSlp - timon1);
    if(dtim > fTimeDelta[cpass]) return;
    // check the slope difference. First do a local fit on cluster 1 near
    // the matching point
    cl2FitMid(allhits, tcl, it1, hiton1, 3);
    if(clChisq > 20.) return;
    // fit parameters are now in clpar. Charge is in fAveChg
    // check for angle consistency
    float dth = fabs(atan(fScaleF * clpar[1]) - atan(fScaleF * tcl[it2].EndSlp));
  if(fDebugWire < 0) std::cout<<"US dtheta "<<dth<<" cut "<<fKinkAngCut[cpass]<<std::endl;
    if(dth > fKinkAngCut[cpass]) return;
    // make a charge ratio cut
    float chgrat = 2 * fabs(fAveChg - tcl[it2].EndChg) / (fAveChg + tcl[it2].EndChg);
  if(fDebugWire < 0)  std::cout<<"US chgrat "<<chgrat<<" cut "<<fChgCut[pass]<<std::endl;
    // ensure that there is a signal on any missing wires at the US end of 1
    bool SigOK = false;
    cl2ChkSignal(allhits, wiron1, timon1, ew2, et2, SigOK);
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
    if(hiton1 > allhits.size() - 1) {
      mf::LogError("ClusterCrawler")<<"ChkMerge12 bad hiton1 "<<hiton1;
      return;
    }
    timon1 = allhits[hiton1].Time;
    dtim = fabs(bt2 - (wiron1 - bw2) *tcl[it2].BeginSlp - timon1);
    if(dtim > fTimeDelta[cpass]) return;
    cl2FitMid(allhits, tcl, it1, hiton1, -3);
    if(clChisq > 20.) return;
    // check for angle consistency
    dth = fabs(atan(fScaleF * clpar[1]) - atan(fScaleF * tcl[it2].BeginSlp));
  if(fDebugWire < 0) std::cout<<"DS dtheta "<<dth<<" cut "<<fKinkAngCut[cpass]<<std::endl;
    if(dth > fKinkAngCut[cpass]) return;
    // make a charge ratio cut
    chgrat = 2 * fabs(fAveChg - tcl[it2].BeginChg) / (fAveChg + tcl[it2].BeginChg);
  if(fDebugWire < 0)  std::cout<<"DS chgrat "<<chgrat<<" cut "<<fChgCut[pass]<<std::endl;
    // ensure that there is a signal on any missing wires at the US end of 1
    SigOK = false;
    cl2ChkSignal(allhits, wiron1, timon1, bw2, bt2, SigOK);
  if(fDebugWire < 0)  std::cout<<"DS SigOK? "<<SigOK<<std::endl;
    if( !SigOK ) return;

  if(fDebugWire < 0)  std::cout<<"Merge em"<<std::endl;
    // success. Merge them
    cl2DoMerge(allhits, tcl, vtx, it1, it2, 100);
    didit = true;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2DoMerge(
      std::vector<CCHitFinderAlg::CCHit>& allhits, 
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
      unsigned short it1, unsigned short it2, short ProcCode)
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
      unsigned short wire = allhits[hit].WireNum;
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
      unsigned short wire = allhits[hit].WireNum;
      unsigned short index = wire - lowire;
      if(index > veclen) {
        throw cet::exception("ClusterCrawler")<<"DoMerge bad index "<<index;
        return;
      }
      if(wirehit[index] < 0) {
        wirehit[index] = hit;
      } else {
        // a hit from cluster cl1 is on this wire. Free up cl2 hit for later use
        unsigned short index = hit - fFirstHit;
        hiterr2[index] = -fabs(hiterr2[index]);
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
          cl2Fit(allhits);
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
    cl2Fit(allhits);
    if(clChisq > 99.) {
      std::cout<<"cl2DoMerge bad End fit "<<clChisq<<std::endl;
      return;
    }
    clEndSlp = clpar[1];
    clEndSlpErr = clparerr[1];
    clEndChg = fAveChg;

    clStopCode = cl1.StopCode;
    // append it to the tcl vector
    cl2TmpStore(allhits, tcl);
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
    // delete a vertex between these two?
    if(tcl[it1].BeginVtx >= 0) {
      if(tcl[it1].BeginVtx == tcl[it2].EndVtx) {
        vtx[tcl[it2].EndVtx].Wght = -1;
        tcl[it1].BeginVtx = -99;
        tcl[it2].EndVtx = -99;
      }
    } else if(tcl[it1].EndVtx >= 0) {
      if(tcl[it1].EndVtx == tcl[it2].BeginVtx) {
        vtx[tcl[it2].BeginVtx].Wght = -1;
        tcl[it1].EndVtx = -99;
        tcl[it2].BeginVtx = -99;
      }
    }
    // set the vertex assignments negative so that this new cluster
    // can be re-assigned to the appropriate vertex
    tcl[itnew].BeginVtx = -99;
    tcl[itnew].EndVtx = -99;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2Print(std::vector<CCHitFinderAlg::CCHit>& allhits, 
     std::vector<ClusterStore>& tcl)
  {
    // prints clusters to the screen for code development
    std::cout<<"  ID CTP nht Stop  Proc  beg_W:H  begT  bTheta Therr begChg end_W:H  endT  eTheta Therr endChg";
    std::cout<<" bVx eVx";
    std::cout<<std::endl;
    for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
      if(fDebugPlane >= 0 && fDebugPlane != tcl[ii].CTP) continue;
      std::vector<unsigned short>::const_iterator ihtb = tcl[ii].tclhits.begin();
      unsigned short hitb = *ihtb;
      std::vector<unsigned short>::const_iterator ihte = tcl[ii].tclhits.end()-1;
      unsigned short hite = *ihte;
      std::cout<<std::right<<std::setw(4)<<tcl[ii].ID;
      std::cout<<std::right<<std::setw(3)<<tcl[ii].CTP;
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
      float arg = fScaleF * tcl[ii].BeginSlp;
      float theta = atan(arg);
      std::cout<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<theta;
      float thetaerr =  fScaleF * tcl[ii].BeginSlpErr / (1 + arg * arg);
      std::cout<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<thetaerr;
      std::cout<<std::right<<std::setw(5)<<(short)tcl[ii].BeginChg;
      std::cout<<std::right<<std::setw(6)<<tcl[ii].EndWir<<":"<<hite;
      if(hite < 10) {
        std::cout<<"   ";
      } else if(hite < 100) {
        std::cout<<"  ";
      } else if(hite < 1000) std::cout<<" ";
      std::cout<<std::right<<std::setw(6)<<(short)tcl[ii].EndTim;
      arg = fScaleF * tcl[ii].EndSlp;
      theta = atan(arg);
      std::cout<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<theta;
      thetaerr =  fScaleF * tcl[ii].EndSlpErr / (1 + arg * arg);
      std::cout<<std::right<<std::setw(7)<<std::fixed<<std::setprecision(2)<<thetaerr;
      std::cout<<std::right<<std::setw(5)<<(short)tcl[ii].EndChg;
      std::cout<<std::right<<std::setw(5)<<tcl[ii].BeginVtx;
      std::cout<<std::right<<std::setw(5)<<tcl[ii].EndVtx;
      std::cout<<std::endl;
    } // ii
  } // cl2Print

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2TmpGet(std::vector<CCHitFinderAlg::CCHit>& allhits,
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
      clCTP = tcl[it1].CTP;
      fcl2hits = tcl[it1].tclhits;
    }


/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2TmpStore(std::vector<CCHitFinderAlg::CCHit>& allhits, 
     std::vector<ClusterStore>& tcl)
  {

    if(fcl2hits.size() < 3) {
//      mf::LogError("ClusterCrawler")<<"cl2TmpStore trying to store crazy cluster";
      return;
    }
    
    ++NClusters;

    // flag all the hits as used
    for(unsigned short it = 0; it < fcl2hits.size(); ++it) {
      unsigned short hit = fcl2hits[it];
      if(hit > allhits.size() - 1) {
        std::cout<<"cl2TmpStore bad hit "<<hit<<std::endl;
        return;
      }
      if(allhits[hit].Charge < 0) {
        std::cout<<"Trying to use obsolete hit "<<hit;
        std::cout<<" on wire "<<allhits[hit].WireNum;
        std::cout<<" on cluster "<<NClusters;
        std::cout<<" in plane "<<plane<<std::endl;
      }
      unsigned short index = hit - fFirstHit;
      hiterr2[index] = -fabs(hiterr2[index]);
    }

    // ensure that the cluster begin/end info is correct

    // define the begin/end charge if it wasn't done already
    if(clEndChg < 0.) {
      // use the next to the last two hits. The End hit may have low charge
      unsigned int ih0 = fcl2hits.size() - 2;
      unsigned int hit = fcl2hits[ih0];
      clEndChg = allhits[hit].Charge;
      hit = fcl2hits[ih0 - 1];
      clEndChg += allhits[hit].Charge;
      clEndChg = clEndChg / 2.;
    }
    if(clBeginChg < 0.) {
      // use the 2nd and third hit. The Begin hit may have low charge
      unsigned int hit = fcl2hits[1];
      clBeginChg = allhits[hit].Charge;
      hit = fcl2hits[2];
      clBeginChg += allhits[hit].Charge;
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
    clstr.BeginWir    = allhits[hitb].WireNum;
    clstr.BeginTim    = allhits[hitb].Time;
    clstr.BeginChg    = clBeginChg;
    clstr.EndSlp      = clEndSlp;
    clstr.EndSlpErr   = clEndSlpErr;
    clstr.EndWir      = allhits[hite].WireNum;
    clstr.EndTim      = allhits[hite].Time;
    clstr.EndChg      = clEndChg;
    clstr.StopCode    = clStopCode;
    clstr.ProcCode    = clProcCode;
    clstr.Assn        = clAssn;
    clstr.BeginVtx    = -99;
    clstr.EndVtx      = -99;
    clstr.CTP         = clCTP;
    clstr.tclhits     = fcl2hits;
    tcl.push_back(clstr);
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2FollowUS(std::vector<CCHitFinderAlg::CCHit>& allhits)
  {
    // follow the cluster upstream

    if(fcl2hits.size() < 2) return;

    unsigned short dhit = fcl2hits[0];
    unsigned short dwir = allhits[dhit].WireNum;
  if(fDebugWire > 0 && fDebugHit > 0) {
    prt = ((short)dwir == fDebugWire && (short)dhit == fDebugHit);
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
    if(lasthit > allhits.size() - 1) {
      std::cout<<"cl2FollowUS bad lasthit "<<lasthit<<std::endl;
    }
    unsigned short lastwire = allhits[lasthit].WireNum;
  if(prt) std::cout<<"cl2FollowUS: last wire "<<lastwire<<" hit "<<lasthit<<std::endl;
    // keep a log of the fit chisq
    std::vector<float> chifits;
    
    for(unsigned short nextwire = lastwire-1; nextwire >= fFirstWire; --nextwire) {
  if(prt) std::cout<<"cl2FollowUS: next wire "<<nextwire<<std::endl;
      // add hits and check for PH and width consistency
      cl2AddHit(allhits, nextwire, HitOK, SigOK);
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
            cl2Fit(allhits);
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
        cl2Fit(allhits);
        if(clChisq > 99.) {
          if(fcl2hits.size() < 3) return;
          // a fit error occurred. Lop off the leading hit and refit
  if(prt) std::cout<<"Fit failed "<<std::endl;
          fcl2hits.pop_back();
          cl2Fit(allhits);
          if(clChisq > 99.) {
            // something really bad happened. Bail out
            fcl2hits.clear();
            return;
          }
          continue;
        }
        // monitor the onset of a kink. Look for a progressive increase
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
            unsigned short ih2 = ih0 - 2;
            unsigned short hit2 = fcl2hits[ih2];
            float dt02 = allhits[hit2].Time - allhits[hit0].Time;
            float dw02 = allhits[hit2].WireNum - allhits[hit0].WireNum;
            float th02 = atan( fScaleF * dt02 / dw02);
            // and the 3rd and 5th hit
            unsigned short ih3 = ih2 - 1;
            unsigned short hit3 = fcl2hits[ih3];
            unsigned short ih5 = ih3 - 2;
            unsigned short hit5 = fcl2hits[ih5];
            float dt35 = allhits[hit5].Time - allhits[hit3].Time;
            float dw35 = allhits[hit5].WireNum - allhits[hit3].WireNum;
            float th35 = atan(fScaleF * dt35 / dw35);
            float dth = fabs(th02 - th35);
  if(prt) std::cout<<" Kink angle "<<std::setprecision(3)<<dth<<std::endl;
            // cut on the allowed kink angle
            if(dth > fKinkAngCut[pass]) {
  if(prt) std::cout<<"stopped tracking "<<std::endl;
              // kill the last 3 hits, refit and return
              for(short jj = 0; jj < 3; ++jj) {
                fcl2hits.pop_back();
              }
              cl2Fit(allhits);
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
          // remove the last few hits if there is a systematic increase in chisq and re-fit
          // long tracks only
          if(fcl2hits.size() > 10) {
            unsigned short cfsize = chifits.size() - 1;
            for(unsigned short nlop = 0; nlop < 4; ++nlop) {
              float chirat = chifits[cfsize - nlop] / chifits[cfsize - nlop - 1];
              if(chirat > 1.1) fcl2hits.pop_back();
            } // nlop
          } else {
            // just remove the last hit
            fcl2hits.pop_back();
          } // fcl2hits.size
          cl2Fit(allhits);
          if(clChisq > 99.) {
            std::cout<<"cl2FollowUS bad fit after pop_back "<<clChisq<<std::endl;
            return;
          }
  if(prt) std::cout<<"Bad fit chisq - removed hits"<<std::endl;
          clStopCode = 4;
          break;
        } // clChisq > fChiCut[pass]
      } // !HitOK check
    } // nextwire
    
    // find the US wire
    unsigned short ih0 = fcl2hits.size() - 1;
    unsigned short hit0 = fcl2hits[ih0];
    unsigned short uswir = allhits[hit0].WireNum;
    // find the wire fMinWirAfterSkip[pass] indices DS
    unsigned short ihskp = ih0 - fMinWirAfterSkip[pass];
    unsigned short hitskp = fcl2hits[ihskp];
    unsigned short wirskp = allhits[hitskp].WireNum;
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
  void ClusterCrawlerAlg::cl2FitMid(std::vector<CCHitFinderAlg::CCHit>& allhits,
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
        if(ihit > allhits.size()-1) {
          mf::LogError("ClusterCrawler")<<"FitMid bad ihit "<<ihit;
          return;
        }
        // look for the desired first hit. Use this as the origin wire
        if(ihit == ihtin) {
          UseEm = true;
          wire0 = allhits[ihit].WireNum;
        }
        // use hits after finding the first desired hit
        if(UseEm) {
          unsigned short wire = allhits[ihit].WireNum;
          xwir.push_back(wire - wire0);
          ytim.push_back(allhits[ihit].Time);
          // pass the error^2 to the fitter
	  ytimerr2.push_back(fabs(hiterr2[ihit - fFirstHit]));
          fAveChg += allhits[ihit].Charge;
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
        if(ihit > allhits.size()-1) {
          mf::LogError("ClusterCrawler")<<"FitMid bad ihit "<<ihit;
          return;
        }
        // look for the desired first hit. Use this as the origin wire
        if(ihit == ihtin) {
          UseEm = true;
          wire0 = allhits[ihit].WireNum;
        }
        // use hits after finding the first desired hit

        if(UseEm) {
          unsigned short wire = allhits[ihit].WireNum;
          xwir.push_back(wire - wire0);
          ytim.push_back(allhits[ihit].Time);
	  ytimerr2.push_back(fabs(hiterr2[ihit - fFirstHit]));
          fAveChg += allhits[ihit].Charge;
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
  void ClusterCrawlerAlg::cl2Fit(std::vector<CCHitFinderAlg::CCHit>& allhits)
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
      unsigned short wire = allhits[ihit].WireNum;
      if(first) {
        wire0 = wire;
        first = false;
  if(prt) std::cout<<"cl2Fit W:H ";
      }
  if(prt) std::cout<<wire<<":"<<ihit<<" ";
      xwir.push_back(wire - wire0);
      ytim.push_back(allhits[ihit].Time);
      ytimerr2.push_back(angfactor * fabs(hiterr2[ihit - fFirstHit]));
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
  void ClusterCrawlerAlg::cl2FitChg(std::vector<CCHitFinderAlg::CCHit>& allhits)
  {
    // Fits the charge of hits on the fcl2hits vector to a line, or simply
    // uses the average of 1 or 2 hits as determined by NHitsAve

    unsigned short ih0 = fcl2hits.size() - 1;
    
    // don't find the average charge --> no charge cut is made
    if(fNHitsAve[pass] < 1) return;
    
    if(fNHitsAve[pass] < 2) {
      // simply use the charge and width the last hit
      fAveChg = allhits[fcl2hits[ih0]].Charge;
    } else if(fNHitsAve[pass] == 2) {
      // average the last two points if requested
      fAveChg = (allhits[fcl2hits[ih0]].Charge + 
                 allhits[fcl2hits[ih0 - 1]].Charge) / 2.;
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
        unsigned short wire = allhits[fcl2hits[ii]].WireNum;
        if(first) {
          wire0 = wire;
          first = false;
        }
        xwir.push_back((float)(wire - wire0));
        float chg = allhits[fcl2hits[ii]].Charge;
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
  void ClusterCrawlerAlg::cl2AddHit(std::vector<CCHitFinderAlg::CCHit>& allhits,
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
    unsigned short wire0 = allhits[ih1].WireNum;
    float bigchgcut = 2 * fChgCut[pass];
    bool lasthitbig = ( (allhits[ih1].Charge / fAveChg) > bigchgcut);
    float lowchgcut = -1.5 * fChgCut[pass];
    bool lasthitlow = ( (allhits[ih1].Charge / fAveChg) < fChgCut[pass]);

    // find the expected time of a hit on this wire
    float prtime = clpar[0] + (kwire - wire0) * clpar[1];
    // max number of time ticks between projected cluster and hit position
    float best = 50.;
    short imbest = -1;
    SigOK = false;
    for(unsigned short khit = firsthit; khit < lasthit; ++khit) {
  if(prt) {
    std::cout<<"cl2AddHit: Check W:H "<<kwire<<":"<<khit<<" time "<<(int)allhits[khit].Time;
    std::cout<<" prtime "<<(short)prtime<<" fAveChg "<<(int)fAveChg;
    std::cout<<" lasthitlow "<<lasthitlow<<" lasthitbig "<<lasthitbig;
    std::cout<<" hitwid "<<hitwid[khit - fFirstHit]<<std::endl;
  }

      if(prtime < allhits[khit].Time + hitwid[khit - fFirstHit] && 
         prtime > allhits[khit].Time - hitwid[khit - fFirstHit]) SigOK = true;
         
      // ignore obsolete hits
      if(allhits[khit].Charge < 0) continue;
      // ignore used hits
      if(hiterr2[khit - fFirstHit] < 0) continue;

      float timediff = (allhits[khit].Time - prtime);

      // make hit charge and width cuts
      // don't make a charge ratio cut until fAveChg is defined
      if(fAveChg > 0.) {
        float chgrat = (allhits[khit].Charge - fAveChg) / fAveChg;
  if(prt) std::cout<<" Chg "<<(int)allhits[khit].Charge<<" chgrat "<<chgrat<<" SigOK "<<SigOK<<std::endl;
        if(lasthitlow) {
          // last hit added was low and this one is as well. If this is a hit
          // that will likely be selected
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
                cl2Fit(allhits);
              }
              // stop looking for hits
              return;
            } // fabs(timediff) < 3
          } // chgrat < -fChgCut[pass]s
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
    // Found a close hit. Check the chisq
    float prtimerr2 = clparerr[0]*clparerr[0] + fabs(kwire-wire0)*clparerr[1]*clparerr[1];
  if(prt) {
    std::cout<<" clerr "<<prtimerr2<<" hiterr "<<hiterr2[imbest - fFirstHit];
    std::cout<<" best "<<best<<std::endl;
  }
    // apply an angle dependent scale factor to the hit error
    float angfactor = 1;
    if(clpar[1] != 0) angfactor = 2 - 1/(1 + fabs(clpar[1]));
    float err = sqrt(prtimerr2 + angfactor * fabs(hiterr2[imbest - fFirstHit]));
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
        cl2FitChg(allhits);
      } else if(!lasthitlow && !lasthitbig) {
        // note that we take the absolute value here
        float chgrat = fabs(allhits[imbest].Charge - fAveChg) / fAveChg;
        if(chgrat < fChgCut[pass]) cl2FitChg(allhits);
      }
    } else {
      if(prt) std::cout<<" >>Bad chisq "<<numsig2<<std::endl;
      HitOK = false;
    }
  }

//////////////////////////////////////
    void ClusterCrawlerAlg::cl2VtxFit(std::vector<ClusterStore>& tcl,
        std::vector<VtxStore>& vtx, unsigned short iv, float& ChiDOF)
    {
      
      std::vector<float> x;
      std::vector<float> y;
      std::vector<float> ey2;
      
      for(unsigned short icl = 0; icl < tcl.size(); ++icl) {
        if(tcl[icl].ID < 0) continue;
        if(tcl[icl].EndVtx == iv) {
          x.push_back(tcl[icl].EndSlp);
          float arg = tcl[icl].EndSlp * tcl[icl].EndWir - tcl[icl].EndTim;
          y.push_back(arg);
          if(tcl[icl].EndSlpErr > 0.) {
            arg = tcl[icl].EndSlpErr * tcl[icl].EndWir;
          } else {
            arg = .01 * tcl[icl].EndWir;
          }
          ey2.push_back(arg * arg);
        } else if(tcl[icl].BeginVtx == iv) {
          x.push_back(tcl[icl].BeginSlp);
          float arg = tcl[icl].BeginSlp * tcl[icl].BeginWir - tcl[icl].BeginTim;
          y.push_back(arg);
          if(tcl[icl].BeginSlpErr > 0.) {
            arg = tcl[icl].BeginSlpErr * tcl[icl].BeginWir;
          } else {
            arg = .01 * tcl[icl].BeginWir;
          }
          ey2.push_back(arg * arg);
        }
      } // ii
      if(x.size() < 2) {
//        std::cout<<"cl2VtxFit: Cluster-Vertex assigment error "<<plane;
//        std::cout<<" iv "<<iv<<std::endl;
        vtx[iv].Wght = -1;
        return;
      }
      
      float tv = 0.;
      float tverr = 0.;
      float wv = 0.;
      float wverr = 0.;
      LinFit(x, y, ey2, tv, wv, tverr, wverr, ChiDOF);
      if(ChiDOF < 5) {
        vtx[iv].Wire = (int)(wv + 0.5);
        vtx[iv].Time = -tv;
      }
//  std::cout<<"VtxFit wv "<<wv<<" tv "<<-tv<<" chi "<<ChiDOF<<" size "<<x.size();
//  std::cout<<" Fit err "<<wverr<<" "<<tverr<<std::endl;
    } // cl2VtxFit

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
                      2 * (A*sumy + B*sumxy - A*B*sumx)) / ndof;
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

//////////////////////////////////
    void ClusterCrawlerAlg::GetHitRange(std::vector<CCHitFinderAlg::CCHit>& allhits,
      unsigned short CTP, 
      std::vector< std::pair<short, short> >& WireHitRange,
      unsigned short& firstwire, unsigned short& lastwire)
    {
      art::ServiceHandle<geo::Geometry> geom;
      // fills the WireHitRange vector for the supplied Cryostat/TPC/Plane code
      bool first = true;
      lastwire = 0;
      unsigned short firsthit = 0;
      unsigned int cst = CTP / 100;
      unsigned int tpc = (CTP - 100 * cst) / 10;
      unsigned int pla = CTP - 100 * cst - 10 * tpc;
      unsigned short lasthit = 0;
      // find the first and last wire with a hit
      for(unsigned short hit = 0; hit < allhits.size(); ++hit) {
        art::Ptr<recob::Wire> theWire = allhits[hit].Wire;
        uint32_t channel = theWire->Channel();
        std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
        if(wids[0].Plane != pla) continue;
        if(wids[0].TPC != tpc) continue;
        if(wids[0].Cryostat != cst) continue;
        unsigned short theWireNum = allhits[hit].WireNum;
        if(theWireNum != allhits[hit].WireNum) {
          std::cout<<"GetHitRange WireNum screwup "<<std::endl;
        }
        if(first) {
          firsthit = hit;
          firstwire = theWireNum;
          first = false;
        }
        lastwire = theWireNum;
        lasthit = hit;
      } //hit

      // now we can define the WireHitRange vector.
      // start by defining the "no hits on wire" condition
      short sflag = -2;
      for(unsigned short wire = firstwire; wire <= lastwire; ++wire) {
        WireHitRange.push_back(std::make_pair(sflag, sflag));
      }
      // overwrite with the "dead wires" condition
      filter::ChannelFilter cf;
      sflag = -1;
      for(unsigned short wire = firstwire+1; wire < lastwire; ++wire) {
        uint32_t chan = geom->PlaneWireToChannel((int)pla,(int)wire,(int)tpc,(int)cst);
        // remember to offset references to WireHitRange by the FirstWire
        unsigned short index = wire - firstwire;
        if(cf.BadChannel(chan)) WireHitRange[index] = std::make_pair(sflag, sflag);
      }
          
      lastwire = firstwire;
      unsigned short thishit = firsthit;
      unsigned short lastfirsthit = firsthit;
      // next overwrite with the index of the first/last hit on each wire
      for(unsigned short hit = firsthit; hit <= lasthit; ++hit) {
        CCHitFinderAlg::CCHit& theHit = allhits[hit];
        unsigned short thiswire = theHit.WireNum;
        if(thiswire > lastwire) {
          unsigned short index = lastwire - firstwire;
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
      } //hit
      // define for the last wire
      unsigned short index = lastwire - firstwire;
      short itmp1 = lastfirsthit;
      short itmp2 = thishit;
      WireHitRange[index] = std::make_pair(itmp1,itmp2);
    } // GetHitRange

/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2SortByLength(std::vector<ClusterStore>& tcl,
        std::map<unsigned short, unsigned short>& sortindex)
    {
      // sorts the temporary cluster vector by decreasing number of hits,
      // while ignoring abandoned clusters. Returns index map with the
      // sort order
      
      // form a vector of pairs of the number of hits and the index
      std::vector< std::pair<unsigned short, unsigned short> > index;
      for(unsigned short ii = 0; ii < tcl.size(); ++ii) {
        if(tcl[ii].ID > 0 && tcl[ii].CTP == clCTP) 
          index.push_back(std::make_pair(tcl[ii].tclhits.size(),ii));
      }
      std::sort(index.begin(), index.end(), SortByLen);
      sortindex.clear();
      for(unsigned short ii = 0; ii < index.size(); ++ii) {
       sortindex[ii]=index[ii].second;
      }
      return; 
    }


} // namespace cluster
