/////////////////////////////////////////////////////////////////////
///
/// ClusterCrawlerAlg class
///
/// Bruce Baller, baller@fnal.gov
///
/// Algorithms for crawling along a string of hits to make a cluster
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
//#include "art/Framework/Services/Optional/TFileService.h" 
//#include "art/Framework/Services/Optional/TFileDirectory.h" 
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
  }

  void ClusterCrawlerAlg::reconfigure(fhicl::ParameterSet const& pset)
  { 
    fNumPass = 2;
    if(fNumPass > 5) {
      mf::LogError("ClusterCrawler") << "Too many passes specified";
      return;
    }
    // first pass settings
    fMaxHitsFit[0] = 20; // First pass - lots of fitted hits --> high momentum tracks
    fMinHits[0] = 10; // min number of hits on the cluster
    // second pass settings
    fMaxHitsFit[1] = 3; // Second pass - fewer fitted hits --> low momentum tracks
    fMinHits[1] = 3; // min number of hits on the cluster
    // settings used for all passes
    fHitErrFac = 0.3; // hit time error = fHitErrFac * (EndTime - PeakTime)
    fChiCut  = 2.;   // stop adding hits to clusters if fChiCut is reached
    fSigCut = 0.5;  // max fractional hit width difference for adding cluster hits
    fChgCut = 0.7;  // max fractional hit charge difference for adding cluster hits
    fMaxWirSkip = 10;  // max number of dead/occupied wires that can be skipped
                      // while following the cluster
    // merging cuts
    fDoMerge = false; // run the merging code?
    fWirDelta = 3;  // max wire difference between start/end for merging
    fSlpDelta = 10; // max slope difference for merging
    fChgDelta = 1.0;  // max charge ratio difference for merging
    fTimDelta = 4;  // max time difference for matching
  }

  void ClusterCrawlerAlg::beginJob()
  {
  }

  bool SortByLowHit(int i, int j) {return ((i > j));}

  void ClusterCrawlerAlg::RunCrawler(art::PtrVector<recob::Hit>& plnhits,int plane,
                std::vector<ClusterStore>& tcl)
  {
    // The hit collection is assumed to be sorted in increasing wire order

    if(plnhits.size() <= 2) return;
    
    art::PtrVector<recob::Hit>::const_iterator it = plnhits.begin();
    fFirstWire = (*it)->WireID().Wire;
    std::cout<<"FirstWire "<<fFirstWire;
    it = plnhits.end()-1;
    int LastWire = (*it)->WireID().Wire;
    std::cout<<" LastWire "<<LastWire<<std::endl;
    
    prt = false;

    // find dead wires in this region
    filter::ChannelFilter cf;
    art::ServiceHandle<geo::Geometry> geom;
    for(int wire = fFirstWire+1; wire < LastWire; wire++) {
      unsigned int pln = plane;
      unsigned int wir = wire;
      uint32_t chan = geom->PlaneWireToChannel(pln,wir);
      // set the FirstWirHit negative
      if(cf.BadChannel(chan)) FirstWirHit[wire] = -1;
    }
    
    int wire = -1;
    int iht = 0;
    for(art::PtrVector<recob::Hit>::const_iterator hitIter = plnhits.begin();
          hitIter < plnhits.end(); ++hitIter) {
      int thiswire = (*hitIter)->WireID().Wire;
      // calculate the hit position uncertainty
      double arg = fHitErrFac *((*hitIter)->EndTime() - (*hitIter)->PeakTime());
      hiterr2[iht] = arg * arg;
  if(prt) {
  std::cout<<"Wire "<<thiswire<<" hit "<<iht<<" Time "<<(int)(*hitIter)->PeakTime();
  std::cout<<" Chg "<<(int)(*hitIter)->Charge();
  std::cout<<" Wid "<<(*hitIter)->EndTime()-(*hitIter)->PeakTime();
  std::cout<<" Mult "<<(*hitIter)->Multiplicity()<<std::endl;
  }
      if(thiswire > wire) {
        FirstWirHit[thiswire] = iht;
        wire = thiswire;
      } else if(thiswire < wire) {
        mf::LogError("ClusterCrawler")<<"ERROR: Hits not sorted!!";
        return;
      }
      iht++;
    }
    FirstWirHit[LastWire+1] = plnhits.size();
    
    prt = false;
    
    // use the ROOT linear fitter. Define 1st order polynomial w 2 params
    TLinearFitter *lf =  new TLinearFitter(2);
    lf->SetFormula("pol1");
    // determine the length of the TLinearFitter arrays
    int maxhits = 0;
    for(int ii = 0; ii < fNumPass; ii++) {
      if(fMaxHitsFit[ii] > maxhits) maxhits = fMaxHitsFit[ii];
    }
    xwir = new Double_t[maxhits];
    ytim = new Double_t[maxhits];
    ytimerr = new Double_t[maxhits];

    unsigned int nHitsUsed = 0;
    bool AllDone = false;
    for(int thispass = 0; thispass < fNumPass; thispass++) {
      pass = thispass;
      int nClusters = 0;
      // look for a starting cluster that spans a block of wires
      int span = 3;
      if(fMinHits[pass] < span) span = fMinHits[pass];
      for(int iwire = LastWire; iwire > fFirstWire; iwire--) {
        int ifirsthit = FirstWirHit[iwire];
        // skip bad wires
        if(ifirsthit <= 0) continue;
        int ilasthit = FirstWirHit[iwire+1];
        // find the index of the last hit on this wire if it is adjacent
        // to a dead wire
        if(ilasthit <= 0) {
          int kwire = iwire + 2;
          while(ilasthit <= 0) {
            ilasthit = FirstWirHit[kwire];
            kwire++;
          }
        }
        for(int ihit = ifirsthit; ihit < ilasthit; ihit++) {
          bool ClusterAdded = false;
          if(hiterr2[ihit] <= 0) continue;
          // Start a cluster if it spans span wires
          for(int jwire = iwire - span + 1; jwire < iwire; jwire++) {
            int jfirsthit = FirstWirHit[jwire];
            if(jfirsthit <= 0) continue;
            int jlasthit = FirstWirHit[jwire+1];
            if(jlasthit <= 0) {
              int kwire = jwire + 2;
              while(jlasthit <= 0) {
                jlasthit = FirstWirHit[kwire];
                kwire++;
              }
            }
  if(prt) {
    std::cout<<"Last wire "<<iwire<<" Hits "<<ifirsthit;
    std::cout<<" to "<<ilasthit-1;
    std::cout<<" >> First wire "<<jwire<<" Hits "<<jfirsthit;
    std::cout<<" to "<<jlasthit-1<<std::endl;
  }
           for(int jhit = jfirsthit; jhit < jlasthit; jhit++) {
              if(hiterr2[jhit] < 0) continue;
  if(prt) std::cout<<"jhit "<<jhit<<" "<<plnhits[jhit]->PeakTime()<<std::endl;
              // start a cluster with these two hits
              fcl2hits.clear();
              // set to nonsensical values
              fAveWid = -1.;
              fAveChg = -1.;
//  prt = (ihit == 218);
              fcl2hits.push_back(ihit);
              fcl2hits.push_back(jhit);
              // set hiterr2 negative to indicate the hit is used
              hiterr2[ihit] = -hiterr2[ihit];
              hiterr2[jhit] = -hiterr2[jhit];
              // define the fit origin. Use the upstream hit
              wire0 = jwire;
              cl2Fit(plnhits);
              // kill it if something bad happened in the fitter
              if(clchisq > 5) {
                cl2Kill(plnhits);
                continue;
              }
              // now look for hits to add on the intervening wires
              bool SigOK = true;
              bool clok = true;
              for(int kwire = jwire+1; kwire < iwire; kwire++) {
                int nadd = cl2AddHit(plnhits, kwire, false, SigOK);
                // no hit added and no nearby hit either
                if(nadd == 0 && !SigOK) {
                  clok = false;
                  break;
                }
              }
  if(prt) {
    std::cout<<">>>>> Starting cluster hits ";
    for(std::vector<int>::iterator it = fcl2hits.begin(); it != fcl2hits.end(); ++it) {
      std::cout<<*it<<" ";
    }
    std::cout<<" nhits "<<fcl2hits.size()<<" span "<<span;
    std::cout<<" clok "<<clok<<std::endl;
  }
              // kill it?
              if((int)fcl2hits.size() < span || !clok) {
                cl2Kill(plnhits);
                continue;
              }
              // sort them by decreasing wire number
              // assume that this is the same as sorting by decreasing 
              // hit number. This only needs to be done on the starting cluster
              // hits will be added in the proper order by cl2Follow
              std::sort(fcl2hits.begin(), fcl2hits.end(), SortByLowHit);
  if(prt) {
    std::cout<<">>>>> Sorted ";
    for(std::vector<int>::iterator it = fcl2hits.begin(); it != fcl2hits.end(); ++it) {
      std::cout<<*it<<" ";
    }
    std::cout<<std::endl;
  }
              // re-fit
              // define the hit origin
              std::vector<int>::reverse_iterator ii = fcl2hits.rbegin();
              int jj = *ii;
              wire0 = plnhits[jj]->WireID().Wire;
              cl2Fit(plnhits);
              if(clchisq > 5) {
                cl2Kill(plnhits);
                continue;
              }
              // save the slope at the start
              clslpstart = (double)clpar[1];
              // follow a trail of hits upstream
              cl2FollowUS(plnhits);
              if((int)fcl2hits.size() < fMinHits[pass]) {
                // kill the cluster if it is too short
                if(prt) std::cout<<"Cluster too short - kill"<<std::endl;
                cl2Kill(plnhits);
              } else {
                // Re-follow the cluster using parameters from the next pass
                // This picks up hits on lower momentum tracks that were missed
                // by the tighter cuts on the current pass
                if(pass < fNumPass) cl2ReFollow(plnhits);
                // store the cluster
                clslpend = (double)clpar[1]; // save the slope at the end
                cl2TmpStore(plnhits, tcl, pass, nClusters); // store the cluster
                nClusters++;
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
    } // pass
    
//    std::cout<<"Before Merging"<<std::endl;
//    cl2Print(plnhits, tcl);
    // try to merge clusters
    if(fDoMerge) cl2ChkMerge(plnhits, tcl);
//    std::cout<<"After Merging"<<std::endl;
//    cl2Print(plnhits, tcl);
    
    return;
  } // RunCrawler
  
/////////////////////////////////////////
    void ClusterCrawlerAlg::cl2ChkMerge(art::PtrVector<recob::Hit>& plnhits,
        std::vector<ClusterStore>& tcl)
    {
      // Try to merge clusters. Clusters that have been subsumed in other
      // clusters, i.e. no longer valid, have ID < 0
      
      if(tcl.size() < 2) return;
      // The size of the ClusterStore vector will increase after merging
      // is done.
      unsigned int tclsize = tcl.size();
      
      for (unsigned int it1 = 0; it1 < tclsize - 1; it1++) {
        bool didmerge = false;
        ClusterStore& cl1 = tcl[it1];
        // ignore already merged clusters
        if(cl1.ID < 0) continue;
        std::vector<int>::const_iterator ihts1 = cl1.tclhits.begin();
        int hits1 = *ihts1;
        int sw1 = plnhits[hits1]->WireID().Wire;
        int st1 = plnhits[hits1]->PeakTime();
        int sc1 = plnhits[hits1]->Charge();
        int ss1 = cl1.slpstart;
        std::vector<int>::const_iterator ihte1 = cl1.tclhits.end()-1;
        int hite1 = *ihte1;
        int ew1 = plnhits[hite1]->WireID().Wire;
        int et1 = plnhits[hite1]->PeakTime();
        int ec1 = plnhits[hite1]->Charge();
        int es1 = cl1.slpend;
        for (unsigned int it2 = it1 + 1; it2 < tclsize; it2++) {
          ClusterStore& cl2 = tcl[it2];
          // ignore already merged clusters
          if(cl2.ID < 0) continue;
          std::vector<int>::const_iterator ihts2 = cl2.tclhits.begin();
          int hits2 = *ihts2;
          int sw2 = plnhits[hits2]->WireID().Wire;
          int st2 = plnhits[hits2]->PeakTime();
          int sc2 = plnhits[hits2]->Charge();
          int ss2 = cl2.slpstart;
          std::vector<int>::const_iterator ihte2 = cl2.tclhits.end()-1;
          int hite2 = *ihte2;
          int ew2 = plnhits[hite2]->WireID().Wire;
          int et2 = plnhits[hite2]->PeakTime();
          int ec2 = plnhits[hite2]->Charge();
          int es2 = cl2.slpend;
          // look for US and DS broken clusters
          // US cluster 1 merge with DS cluster 2?
          if(abs(sw1 - ew2) < fWirDelta) {
            double dchg = abs((double)(sc1 - ec2) / (double)(sc1 + ec2));
            int dslp = abs(ss1 - es2);
            // project sw1,st1 to ew2
            int dtim = abs(st1 + (ew2-sw1)*ss1 - et2);
            if(dchg < fChgDelta && dslp < fSlpDelta && dtim < fTimDelta) {
              cl2DoMerge(plnhits, tcl, it1, it2);
              tclsize = tcl.size();
              didmerge = true;
              continue;
            }
          // look for US and DS broken clusters
          // US cluster 2 merge with DS cluster 1?
          } else if(abs(sw2 - ew1) < fWirDelta) {
            double dchg = abs((double)(sc2 - ec1) / (double)(sc2 + ec1));
            int dslp = abs(ss2 - es1);
            // project sw2,st2 to ew1
            int dtim = abs(st2 + (ew1-sw2)*ss2 - et1);
            if(dchg < fChgDelta && dslp < fSlpDelta && dtim < fTimDelta) {
              cl2DoMerge(plnhits, tcl, it2, it1);
              tclsize = tcl.size();
              didmerge = true;
              continue;
            }
          // look for small cl1 within the wire boundary of cl2
          } else if(sw2 < sw1 && ew2 > ew1) {
            int dslp = abs(ss2 - ss1);
            std::cout<<"cl2: "<<tcl[it2].ID<<" within cl1 "<<tcl[it1].ID;
            std::cout<<" ? dslp "<<dslp<<std::endl;
          // look for small cl2 within the wire boundary of cl1
          } else if(sw1 < sw2 && ew1 > ew2) {
            int dslp = abs(ss2 - ss1);
            std::cout<<"cl1: "<<tcl[it1].ID<<" within cl2 "<<tcl[it2].ID;
            std::cout<<" ? dslp "<<dslp<<std::endl;
          } // start/end wire check
        } // cluster 2
        if(didmerge) continue;
      } // cluster 1
      
      return;
    }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2DoMerge(art::PtrVector<recob::Hit>& plnhits, 
     std::vector<ClusterStore>& tcl, unsigned int it1, unsigned int it2)
  {
    // Merge clusters. Cluster 1 has precedence for assignment of hits
    ClusterStore& cl1 = tcl[it1];
    ClusterStore& cl2 = tcl[it2];
    std::cout<<"DoMerge "<<cl1.ID<<" "<<cl2.ID<<std::endl;
    // ensure that there is only one hit/wire on both clusters
    std::map<int, int> wirehit;
    int hiwire = -99999;
    int lowire = 99999;
    for(std::vector<int>::const_iterator iht = cl1.tclhits.begin();
          iht != cl1.tclhits.end(); ++iht) {
      int hit = *iht;
      int wire = plnhits[hit]->WireID().Wire;
      if(wire < lowire) lowire = wire;
      if(wire > hiwire) hiwire = wire;
      wirehit[wire] = hit;
    }
    for(std::vector<int>::const_iterator iht = cl2.tclhits.begin();
          iht != cl2.tclhits.end(); ++iht) {
      int hit = *iht;
      int wire = plnhits[hit]->WireID().Wire;
      if(wire < lowire) lowire = wire;
      if(wire > hiwire) hiwire = wire;
      if(wirehit[wire] == 0) {
        wirehit[wire] = hit;
      } else {
        // a hit from cluster cl1 is on this wire. Free it up for later use
        int freehit = wirehit[wire];
        hiterr2[freehit] = abs(hiterr2[freehit]);
        wirehit[wire] = hit;
      }
    }
    // make a new cluster
    ClusterStore clnew;
    clnew.ID = cl1.ID + 10000;
    // Transfer into fit vector
    fcl2hits.clear();
    for(int wire = hiwire; wire >= lowire; wire--) {
      if(wirehit[wire] > 0) fcl2hits.push_back(wirehit[wire]);
    }
    // re-fit the cluster using the parameters for the higher pass
    int pass1 = cl1.ID / 1000;
    int pass2 = cl2.ID / 1000;
    if(pass1 > pass2) {
      pass = pass1;
    } else {
      pass = pass2;
    }
    // clean up the cluster
//    cl2Clean(plnhits);
    // re-fit the end of the cluster
    std::vector<int>::reverse_iterator iend = fcl2hits.rbegin();
    int jend = *iend;
    wire0 = plnhits[jend]->WireID().Wire;
    cl2Fit(plnhits);
    clnew.slpend = clslpend;
    // re-fit the start of the cluster
    std::vector<int>::iterator istrt = fcl2hits.begin();
    int jstrt = *istrt;
    wire0 = plnhits[jstrt]->WireID().Wire;
    cl2Fit(plnhits);
    clnew.slpstart = clslpstart;
    clnew.tclhits = fcl2hits;
    tcl.push_back(clnew);
    // mark cl1 and cl2 obsolete
    for(unsigned int ii = 0; ii < tcl.size(); ii++) {
      ClusterStore& clstr = tcl[ii];
      if(clstr.ID == cl1.ID) clstr.ID = -clstr.ID;
      if(clstr.ID == cl2.ID) clstr.ID = -clstr.ID;
    }
    return;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2Clean(art::PtrVector<recob::Hit>& plnhits)
  {
    // clean up the cluster in the fcl2hit array
    // form a map of the hits used on each wire
    std::map<int,int> wirehit;
    int lowire = 99999;
    int hiwire = -99999;
    for(std::vector<int>::iterator it = fcl2hits.begin(); it != fcl2hits.end(); ++it) {
      int hit = *it;
      int wire = plnhits[hit]->WireID().Wire;
      wirehit[wire] = hit;
      if(wire < lowire) lowire = wire;
      if(wire > hiwire) hiwire = wire;
    }
    // form a map of the projected cluster position on each wire
  }


/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2Print(art::PtrVector<recob::Hit>& plnhits, 
     std::vector<ClusterStore>& tcl)
  {
    // prints clusters to the screen for code development
    std::cout<<"Cluster  ID  startW  startH startT startSlp endW endH endT EndSlp" <<std::endl;
    for(unsigned int ii = 0; ii < tcl.size(); ii++) {
      std::vector<int>::const_iterator ihts = tcl[ii].tclhits.begin();
      int hits = *ihts;
      int wirs = plnhits[hits]->WireID().Wire;
      std::vector<int>::const_iterator ihte = tcl[ii].tclhits.end()-1;
      int hite = *ihte;
      int wire = plnhits[hite]->WireID().Wire;
      std::cout<<std::setw(4)<<ii<<std::setw(6)<<tcl[ii].ID;
      std::cout<<std::setw(5)<<wirs<<std::setw(5)<<hits;
      std::cout<<std::setw(6)<<(int)plnhits[hits]->PeakTime();
      std::cout<<std::setw(7)<<std::setprecision(3)<<tcl[ii].slpstart;
      std::cout<<std::setw(4)<<wire<<std::setw(6)<<hite;
      std::cout<<std::setw(6)<<(int)plnhits[hite]->PeakTime();
      std::cout<<std::setw(7)<<tcl[ii].slpend<<std::endl;
    }
  }


/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2TmpStore(art::PtrVector<recob::Hit>& plnhits, 
     std::vector<ClusterStore>& tcl, int pass, int nClusters)
  {
    // store the cluster in the temporary ClusterStore struct
    ClusterStore clstr;
    
    clstr.ID = 1000 * pass + nClusters;
    clstr.slpstart = clslpstart;
    clstr.slpend   = clslpend;
    clstr.tclhits = fcl2hits;
    tcl.push_back(clstr);
    return;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2FollowUS(art::PtrVector<recob::Hit>& plnhits)
  {
    // follow the cluster upstream

  if(prt) {
    std::cout<<"cl2FollowUS hits: ";
    for(std::vector<int>::const_iterator it = fcl2hits.begin(); it != fcl2hits.end(); ++it) {
      std::cout<<*it<<" ";
    }
    std::cout<<std::endl;
  }
    
    // SigOK = true if there is a ADC signal near the projected cluster position
    bool SigOK = true;
    // count the number of missed hits on adjacent wires
    int nmissed = 0;
    std::vector<int>::iterator it = fcl2hits.end() - 1;
    int lasthit = *it;
    int lastwire = plnhits[lasthit]->WireID().Wire;
  if(prt) std::cout<<"cl2FollowUS: last wire "<<lastwire<<" hit "<<lasthit<<std::endl;
    std::vector<double> chifits;
    
    for(int nextwire = lastwire-1; nextwire > fFirstWire-1; --nextwire) {
  if(prt) std::cout<<"cl2FollowUS: next wire "<<nextwire<<std::endl;
      // add hits and check for PH and width consistency
      int nadd = cl2AddHit(plnhits, nextwire, true, SigOK);
  if(prt) std::cout<<"cl2FollowUS: nadd "<<nadd<<" SigOK "<<SigOK<<std::endl;
      // no hit on this wire. Was there a signal or dead wire?
      if(nadd == 0) {
        if(SigOK) {
          nmissed++;
          if(prt && nmissed > fMaxWirSkip) std::cout<<"nmissed break"<<std::endl;
          if(nmissed > fMaxWirSkip) return;
        } else {
          if(prt) std::cout<<"No hit or signal on wire "<<nextwire<<std::endl;
          return;
        }
      } else {
        // update the fit
        // find the origin of the fit
        std::vector<int>::reverse_iterator ii = fcl2hits.rbegin();
        int jj = *ii;
        wire0 = plnhits[jj]->WireID().Wire;
        cl2Fit(plnhits);
        chifits.push_back(clchisq);
        // reset nmissed
        nmissed = 0;
        // check for the onset of a kink
        if(clchisq > fChiCut) {
/*
          int nchk = 0;
          std::cout<<"Badchi ";
          for(std::vector<double>::reverse_iterator ifit = chifits.rbegin();
              ifit != chifits.rend(); ++ifit) {
            if(nchk < 5) std::cout<<(*ifit)<<" ";
          }
          std::cout<<std::endl;
*/
          // remove the last hit and re-fit
          std::vector<double>::iterator imlast = chifits.end() - 1;
          int iht = *imlast;
          hiterr2[iht] = abs(hiterr2[iht]);
          fcl2hits.erase(fcl2hits.end() - 1);
          // find the origin of the fit
          std::vector<int>::reverse_iterator ii = fcl2hits.rbegin();
          int jj = *ii;
          wire0 = plnhits[jj]->WireID().Wire;
          cl2Fit(plnhits);
          return;
        }
      }
    }
    return;
  }

/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2ReFollow(art::PtrVector<recob::Hit>& plnhits)
  {
    // Try to extend a cluster using parameters from the next pass

    if(pass == fNumPass) return;
    
    // store the fcl2hits vector
    std::vector<int> savhits;
    savhits = fcl2hits;
    double savstart = clslpstart;
    double savend = clslpend;
    fcl2hits.clear();
    // look for a missing hit
    std::vector<int>::iterator it = savhits.begin();
    int lasthit = *it;
    int lastwire = plnhits[lasthit]->WireID().Wire;
//    std::cout<<"ReFollow start size= "<<savhits.size();
//    std::cout<<" lastwire "<<lastwire<<std::endl;
//    std::cout<<"lastwire "<<lastwire<<std::endl;
    for(std::vector<int>::iterator it = savhits.begin(); it != savhits.end()-1; ++it) {
      int hit = *it;
      fcl2hits.push_back(hit);
      int wire = plnhits[hit]->WireID().Wire;
//  std::cout<<"wire "<<wire<<" "<<hit<<std::endl;
      if(wire < lastwire - 1) {
//  std::cout<<"missing hit on wire "<<wire<<std::endl;
        // release the hits to the end of the vector
        for(std::vector<int>::iterator ii = it+1; ii != savhits.end(); ++ii) {
          int iht = *ii;
          hiterr2[iht] = abs(hiterr2[iht]);
        }
        break;
      }
      lastwire = wire;
      lasthit = hit;
    }
//    std::cout<<"Start on wire "<<lastwire<<" "<<lasthit<<std::endl;
    // define the params for the next pass
    pass++;
    // re-fit and follow
    std::vector<int>::reverse_iterator ii = fcl2hits.rbegin();
    int jj = *ii;
    wire0 = plnhits[jj]->WireID().Wire;
    cl2Fit(plnhits);
    cl2FollowUS(plnhits);
    pass--;
    std::vector<int>::reverse_iterator kk = savhits.rbegin();
    int iit = *kk;
    int wir = plnhits[iit]->WireID().Wire;
//    std::cout<<"ReFollow done. size= "<<fcl2hits.size();
//    std::cout<<" last wire "<<wir<<std::endl;
    if(fcl2hits.size() < savhits.size()) {
      // failure
      clslpstart = savstart;
      clslpend = savend;
      fcl2hits = savhits;
    }
  }


/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2Fit(art::PtrVector<recob::Hit>& plnhits)
  {
    // Fits the hits on a cluster with origin at wire0 defined by the calling
    // routine. This routine assumes that
    // wires are numbered from lower (upstream) to higher (downstream) and
    // that the hits in the fclhits vector are sorted so that upstream hits
    // are at the end of the vector
    
    // use the ROOT linear fitter
    TLinearFitter *lf =  new TLinearFitter(2,"pol1");
    Int_t nht;
    // fit all hits or truncate?
    if((int)fcl2hits.size() < fMaxHitsFit[pass]) {
      nht = fcl2hits.size();
    } else {
      nht = fMaxHitsFit[pass];
    }

    // load the hits starting at the back end of the fcl2hits vector.
    // These are the most upstream hits
    int iht = 0;
  if(prt) std::cout<<"cl2Fit nhit "<<fcl2hits.size()<<" fit nht "<<nht<<" hits ";
    for(std::vector<int>::reverse_iterator it = fcl2hits.rbegin(); it != fcl2hits.rend(); ++it) {
      int ihit = *it;
  if(prt) std::cout<<ihit<<" ";
      int wire = plnhits[ihit]->WireID().Wire;
      if(wire >= wire0) {
        xwir[iht] = (Double_t)(wire - wire0);
        ytim[iht] = plnhits[ihit]->PeakTime();
        ytimerr[iht] = sqrt(abs(hiterr2[ihit]));
        if(iht == nht) break;
        iht++;
      }
    }
  if(prt) std::cout<<std::endl;

    if(nht < 2) {
      clchisq = 999.;
      return;
    }

    lf->AssignData(nht,1,xwir,ytim,ytimerr);
    // fit the points
    lf->Eval();
    lf->GetParameters(clpar);
    lf->GetErrors(clparerr);
    double chidof = 0;
    if(nht > 2) chidof = lf->GetChisquare() / (nht - 2);

    if(prt) std::cout<<"Fit par "<<clpar[0]<<" "<<clpar[1]<<" chidof "<<chidof;

    // find average PH and width if there are enough hits
    if(fcl2hits.size() > 1) {
      double chgsum = 0.;
      double widsum = 0.;
      int num = 0;
      for(std::vector<int>::reverse_iterator it = fcl2hits.rbegin(); it != fcl2hits.rend(); ++it) {
        int ihit = *it;
        chgsum += plnhits[ihit]->Charge(); // PH using hit area
        widsum += plnhits[ihit]->EndTime() - plnhits[ihit]->PeakTime();
        num++;
        if(num == 4) break;
      }
      fAveWid = widsum / (double)num;
      fAveChg = chgsum / (double)num;
      if(prt) std::cout<<" Ave Wid/Chg "<<fAveWid<<" "<<fAveChg;
    }
    if(prt) std::cout<<std::endl;
    clchisq = chidof;
    return;
  }


/////////////////////////////////////////
  int ClusterCrawlerAlg::cl2AddHit(art::PtrVector<recob::Hit>& plnhits,
        int kwire, bool hitchk, bool SigOK)
  {
    // Add a hit to the cluster if it meets several criteria:
    // similar pulse height to the cluster (if hitchk true)
    // similar hit width to the cluster (if hitchk true)
    // closest hit to the project cluster position.
    // Return SigOK if there is a nearby hit that was missed due to the cuts
    
    int firsthit = FirstWirHit[kwire];
    // skip bad wire, but assume the track was there
    if(firsthit < 0) {
      SigOK = true;
      return 0;
    }
    int lasthit = FirstWirHit[kwire+1];
    // check for dead wire
    if(lasthit <= 0) {
      int lwire = kwire + 2;
      while(lasthit <= 0) {
        lasthit = FirstWirHit[lwire];
        lwire++;
      }
    }
    // find the expected time of a hit on this wire
    double prtime = clpar[0] + (kwire - wire0) * clpar[1];
    // max number of time ticks between projected cluster and hit position
    double best = 10.;
    int imbest = -1;
    SigOK = false;
    for(int khit = firsthit; khit < lasthit; khit++) {
      if(prtime < plnhits[khit]->EndTime() && 
         prtime > plnhits[khit]->StartTime()) SigOK = true;
      if(hiterr2[khit] < 0) continue;
      // make hit charge and width cuts
  if(prt) std::cout<<"cl2AddHit chk wire "<<kwire<<" hit "<<khit;
      if(hitchk) {
        double chgrat = abs(plnhits[khit]->Charge() - fAveChg) / fAveChg;
  if(prt) {
    std::cout<<" Chg "<<(int)plnhits[khit]->Charge()<<" chgrat "<<chgrat;
    if(chgrat > fChgCut) std::cout<<std::endl;
  }
        if(chgrat > fChgCut) continue;
        double sigrat = abs(plnhits[khit]->EndTime() - plnhits[khit]->PeakTime()
             - fAveWid) / fAveWid;
  if(prt) {
    std::cout<<" sigrat "<<sigrat;
    if(sigrat > fSigCut) std::cout<<std::endl;
  }
        if(sigrat > fSigCut) continue;
      }
      double timediff = abs(plnhits[khit]->PeakTime() - prtime);
  if(prt) std::cout<<" time diff "<<timediff;
      if(timediff < best) {
        best = timediff;
        imbest = khit;
      }
  if(prt) std::cout<<" imbest "<<imbest<<std::endl;
    }
    if(imbest < 0) return 0;
    // Found a close hit check the chisq
    double prtimerr2 = clparerr[0]*clparerr[0] + fabs(kwire-wire0)*clparerr[1]*clparerr[1];
  if(prt) {
    std::cout<<"clerr "<<prtimerr2<<" hiterr "<<hiterr2[imbest];
    std::cout<<" best "<<best<<std::endl;
  }
    double err2 = prtimerr2 + hiterr2[imbest];
    // (number of sigma)^2 difference
    double numsig2 = best * best / err2;
    // equivalent to a 3 sigma cut
    if(numsig2 < 9.) {
      fcl2hits.push_back(imbest);
      hiterr2[imbest] = -hiterr2[imbest];
      if(prt) std::cout<<"cl2AddHit ADD wire "<<kwire<<" hit "<<imbest<<" best "<<best<<std::endl;
      return 1;
    } else {
      if(prt) std::cout<<"cl2AddHit bad chisq "<<numsig2<<std::endl;
      return 0;
    }
  }


/////////////////////////////////////////
  void ClusterCrawlerAlg::cl2Kill(art::PtrVector<recob::Hit>& plnhits)
  {
    // kill the current cluster
    for(std::vector<int>::iterator it = fcl2hits.begin(); it != fcl2hits.end(); ++it) {
      int ihit = *it;
      hiterr2[ihit] = abs(hiterr2[ihit]);
    }
    fcl2hits.clear();
    return;
  }

} // namespace cluster
