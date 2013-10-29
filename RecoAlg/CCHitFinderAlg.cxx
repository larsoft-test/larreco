//////////////////////////////////////////////////////////////////////
///
/// CCHitFinder class
///
/// Bruce Baller, baller@fnal.gov
///
/// Find hits for ClusterCrawler and put them in a temporary struct.
/// These hits may be modified by ClusterCrawler before saving them
/// in the event
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

#include "RecoAlg/CCHitFinderAlg.h"

namespace cluster{

//------------------------------------------------------------------------------
  CCHitFinderAlg::CCHitFinderAlg(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  void CCHitFinderAlg::reconfigure(fhicl::ParameterSet const& pset)
  {
    fCalDataModuleLabel = pset.get< std::string  >("CalDataModuleLabel");
    fMinSigInd          = pset.get< float       >("MinSigInd");
    fMinSigCol          = pset.get< float       >("MinSigCol");
    fMinRMSInd          = pset.get< float       >("MinRMSInd");
    fMinRMSCol          = pset.get< float       >("MinRMSCol");
    fMaxBumps           = pset.get< unsigned short >("MaxBumps");
    fMaxXtraHits        = pset.get< unsigned short >("MaxXtraHits");
    fChiSplit           = pset.get< float       >("ChiSplit");
    fChiNorms           = pset.get< std::vector< float > >("ChiNorms");
    fTimeOffsets        = pset.get< std::vector< float > >("TimeOffsets");
    fChgNorms           = pset.get< std::vector< float > >("ChgNorms");

    // stuff these parameters into the hitcut struct so they can be accessed
    // by other CC algs
    hitcuts.MinSigInd = fMinSigInd;
    hitcuts.MinSigCol = fMinSigCol;
    hitcuts.MinRMSInd = fMinRMSInd;
    hitcuts.MinRMSCol = fMinRMSCol;
    hitcuts.ChiSplit  = fChiSplit;
    hitcuts.ChiNorms  = fChiNorms;
    hitcuts.TimeOffsets = fTimeOffsets;
    hitcuts.ChgNorms  = fChgNorms;
  }

//------------------------------------------------------------------------------
  CCHitFinderAlg::~CCHitFinderAlg()
  {
  }
  
  void CCHitFinderAlg::RunCCHitFinder(art::Event & evt) {
  

    allhits.clear();

    // make this accessible to ClusterCrawler_module
    art::Handle< std::vector<recob::Wire> > wireVecHandle;
    evt.getByLabel(fCalDataModuleLabel,wireVecHandle);
    
  std::cout<<"CCHitFinder "<<std::endl;

    // don't expect more than 50% of the max time to have a signal
    unsigned short maxticks = detprop->NumberTimeSamples() / 2 - 4;
    float *ticks = new float[maxticks];
    // define the ticks array used for fitting 
    for(unsigned short ii = 0; ii < maxticks; ++ii) {
      ticks[ii] = ii;
    }
    float *signl = new float[maxticks];

    prt = false;

    for(size_t wireIter = 0; wireIter < wireVecHandle->size(); wireIter++){

      art::Ptr<recob::Wire> theWire(wireVecHandle, wireIter);
      theChannel = theWire->Channel();
      geo::SigType_t SigType = geom->SignalType(theChannel);
      minSig = 0.;
      minRMS = 0.;
      if(SigType == geo::kInduction){
        minSig = fMinSigInd;
        minRMS = fMinRMSInd;
      }//<-- End if Induction Plane
      else if(SigType == geo::kCollection){
        minSig = fMinSigCol;
        minRMS  = fMinRMSCol;
      }//<-- End if Collection Plane


      // minimum number of time samples
      unsigned short minSamples = 2 * minRMS;

      std::vector<geo::WireID> wids = geom->ChannelToWire(theChannel);
      thePlane = wids[0].Plane;
      theWireNum = wids[0].Wire;

      // factor used to normalize the chi/dof fits for each plane
      chinorm = fChiNorms[thePlane];
      timeoff = fTimeOffsets[thePlane];
      ChgNorm = fChgNorms[thePlane];

      // debugging
//  prt = (thePlane == 2 && theWireNum == 1184);

      std::vector<float> signal(theWire->Signal());

      unsigned short nabove = 0;
      unsigned short tstart = 0;
      unsigned short maxtime = signal.size() - 2;
      float maxSig = 0.;
      unsigned short maxSigT = 0;
      for(unsigned short time = 3; time < maxtime; ++time) {
        if(signal[time] > minSig) {
          if(nabove == 0) tstart = time;
          // monitor the max signal and it's time
          if(signal[time] > maxSig) {
            maxSig = signal[time];
            maxSigT = time;
          }
          ++nabove;
        } else {
          // check for a wide enough signal above threshold
          if(nabove > minSamples) {
            // skip this wire if the RAT is too long
            if(nabove > maxticks) break;
            // skip this wire if the max signal is at the beginning
            // --> deconvolution failure
            if(maxSigT < tstart + 2) break;
            unsigned short npt = 0;
            // look for bumps to inform the fit
            bumps.clear();
            for(unsigned short ii = tstart; ii < time; ++ii) {
              signl[npt] = signal[ii];
              if(signal[ii    ] > signal[ii - 1] &&
                 signal[ii - 1] > signal[ii - 2] &&
                 signal[ii    ] > signal[ii + 1] &&
                 signal[ii + 1] > signal[ii + 2]) bumps.push_back(npt);
  if(prt) std::cout<<"signl "<<ii<<" "<<signl[npt]<<std::endl;
              ++npt;
            }
            // just make a crude hit if too many bumps
            if(bumps.size() > fMaxBumps) {
              MakeCrudeHit(npt, ticks, signl);
              StoreHits(tstart, theWire);
              nabove = 0;
              maxSig = 0.;
              continue;
            }
            // start looking for hits with the found bumps
            unsigned short nHitsFit = bumps.size();
            unsigned short nfit = 0;
            chidof = 0.;
            bool HitStored = false;
            unsigned short nMaxFit = bumps.size() + fMaxXtraHits;
            while(nHitsFit <= nMaxFit) {
              FitNG(nHitsFit, npt, ticks, signl);
              // good chisq so store it
              if(chidof < fChiSplit) {
                StoreHits(tstart, theWire);
                HitStored = true;
                break;
              }
              // the previous fit was better, so revert to it and
              // store it
              ++nHitsFit;
              ++nfit;
            } // nHitsFit < fMaxXtraHits
            if( !HitStored && npt < maxticks) {
              // failed all fitting. Make a crude hit
              MakeCrudeHit(npt, ticks, signl);
              StoreHits(tstart, theWire);
            }
          } // nabove > minSamples
          nabove = 0;
          maxSig = 0.;
        } // signal < minSig
      } // time
    } // wireIter

    delete ticks;
    delete signl;

  } //RunCCHitFinder


/////////////////////////////////////////
  void CCHitFinderAlg::FitNG(unsigned short nGaus, unsigned short npt, 
    float *ticks, float *signl)
  {
    // Fit the signal to n Gaussians

    short ndof = npt - 3 * nGaus;
    
    chidof = 9999.;

    if(ndof < 3) return;
    if(bumps.size() == 0) return;

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
    
    TGraph *fitn = new TGraph(npt, ticks, signl);
    TF1 *Gn = new TF1("gn",eqn.c_str());

  if(prt) std::cout<<"FitNG nGaus "<<nGaus<<" nBumps "<<bumps.size()<<std::endl;

    // put in the bump parameters. Assume that nGaus >= bumps.size()
    for(unsigned short ii = 0; ii < bumps.size(); ++ii) {
      unsigned short index = ii * 3;
      unsigned short bumptime = bumps[ii];
      double amp = signl[bumptime];
      Gn->SetParameter(index    , amp);
      Gn->SetParLimits(index, 0., 9999.);
      Gn->SetParameter(index + 1, (double)bumptime);
      Gn->SetParLimits(index + 1, 0, (double)npt);
      Gn->SetParameter(index + 2, (double)minRMS);
      Gn->SetParLimits(index + 2, 1., 3*(double)minRMS);
  if(prt) std::cout<<"Bump params "<<ii<<" "<<(short)amp<<" "<<(int)bumptime<<" "<<(int)minRMS<<std::endl;
    } // ii bumps

    // search for other bumps that may be hidden by the already found ones
    for(unsigned short ii = bumps.size(); ii < nGaus; ++ii) {
      // bump height must exceed minSig
      float big = minSig;
      unsigned short imbig = 0;
      for(unsigned short jj = 0; jj < npt; ++jj) {
        float diff = signl[jj] - Gn->Eval((Double_t)jj, 0, 0, 0);
        if(diff > big) {
          big = diff;
          imbig = jj;
        }
      } // jj
      if(imbig > 0) {
  if(prt) std::cout<<"Found bump "<<ii<<" "<<(short)big<<" "<<imbig<<std::endl;
        // set the parameters for the bump
        unsigned short index = ii * 3;
        Gn->SetParameter(index    , (double)big);
        Gn->SetParLimits(index, 0., 9999.);
        Gn->SetParameter(index + 1, (double)imbig);
        Gn->SetParLimits(index + 1, 0, (double)npt);
        Gn->SetParameter(index + 2, (double)minRMS);
        Gn->SetParLimits(index + 2, 1., 3*(double)minRMS);
      } // imbig > 0
    } // ii 
    
    // W = set weights to 1, N = no drawing or storing, Q = quiet
    // B = bounded parameters
    fitn->Fit(Gn,"WNQB");
    
    // load the fit into a temp vector
    std::vector<double> partmp;
    std::vector<double> partmperr;

    for(unsigned short ipar = 0; ipar < 3 * nGaus; ++ipar) {
      partmp.push_back(Gn->GetParameter(ipar));
      partmperr.push_back(Gn->GetParError(ipar));
    }
    chidof = Gn->GetChisquare() / ( ndof * chinorm);
    

    if(prt) {
      std::cout<<"Fit "<<nGaus<<" chi "<<chidof<<" npars "<<partmp.size()<<std::endl;
      std::cout<<"pars    errs "<<std::endl;
      for(unsigned short ii = 0; ii < partmp.size(); ++ii) {
        std::cout<<ii<<" "<<partmp[ii]<<" "<<partmperr[ii]<<std::endl;
      }
    }

    // ensure that the fit is reasonable
    bool fitok = true;
    for(unsigned short ii = 0; ii < nGaus; ++ii) {
      unsigned short index = ii * 3;
      // ensure that the fitted time is within the signal bounds
      short fittime = partmp[index + 1];
      if(fittime < 0 || fittime > npt - 1) {
        fitok = false;
        break;
      }
      // ensure that the signal peak is large enough
      if(partmp[index] < minSig) {
        fitok = false;
        break;
      }
      // ensure that the RMS is large enough but not too large
      float rms = partmp[index + 2];
      if(rms < minRMS || rms > 4 * minRMS) {
        fitok = false;
        break;
      }
      // ensure that the hits are not too similar in time (< 2 ticks)
      for(unsigned short jj = 0; jj < nGaus; ++jj) {
        if(jj == ii) continue;
        unsigned short jndex = jj * 3;
        float timediff = fabs(partmp[jndex + 1] - partmp[index + 1]);
        if(timediff < 2.) {
          fitok = false;
          break;
        }
      }
      if(!fitok) break;
    }

    if(fitok) {
      par = partmp;
      parerr = partmperr;
    } else {
      chidof = 9999.;
      if(prt) std::cout<<"Bad fit parameters"<<std::endl;
    }
    
    delete fitn;
    delete Gn;
    
    return;
  }

/////////////////////////////////////////
  void CCHitFinderAlg::MakeCrudeHit(unsigned short npt, 
    float *ticks, float *signl)
  {
    // make a single crude hit if fitting failed
    float sumS = 0.;
    float sumST = 0.;
    for(unsigned short ii = 0; ii < npt; ++ii) {
      sumS  += signl[ii];
      sumST += signl[ii] * ticks[ii];
    }
    float mean = sumST / sumS;
    float rms = 0.;
    for(unsigned short ii = 0; ii < npt; ++ii) {
      float arg = ticks[ii] - mean;
      rms += signl[ii] * arg * arg;
    }
    rms = sqrt(rms / sumS);
    float amp = sumS / (Sqrt2Pi * rms);
    par.clear();
  if(prt) std::cout<<"Crude hit Amp "<<(int)amp<<" mean "<<(int)mean<<" rms "<<rms<<std::endl;
    par.push_back(amp);
    par.push_back(mean);
    par.push_back(rms);
    // need to do the errors better
    parerr.clear();
    float amperr = npt;
    float meanerr = sqrt(1/sumS);
    float rmserr = 0.2 * rms;
    parerr.push_back(amperr);
    parerr.push_back(meanerr);
    parerr.push_back(rmserr);
  if(prt) std::cout<<" errors Amp "<<amperr<<" mean "<<meanerr<<" rms "<<rmserr<<std::endl;
    chidof = 9999.;
  }


/////////////////////////////////////////
  void CCHitFinderAlg::StoreHits(unsigned short TStart,
    art::Ptr<recob::Wire>& theWire)
  {
    // store the hits in the struct
    unsigned short nhits = par.size() / 3;
    
    if(nhits == 0) return;
    
    CCHit onehit;
    // lohitid is the index of the first hit that will be added. Hits with
    // Multiplicity > 1 will reside in a block from
    // lohitid to lohitid + numHits - 1
    unsigned short lohitid = allhits.size();
    for(unsigned short hit = 0; hit < nhits; ++hit) {
      unsigned short index = 3 * hit;
      
      onehit.Charge = Sqrt2Pi * par[index] * par[index + 2] / ChgNorm;
      onehit.ChargeErr = SqrtPi * (parerr[index] * par[index + 2] +
                                  par[index] * parerr[index + 2]);
      onehit.Amplitude = par[index];
      onehit.AmplitudeErr = parerr[index];
      onehit.Time = par[index + 1] + TStart + timeoff;
      onehit.TimeErr = parerr[index + 1];
      onehit.RMS = par[index + 2];
      onehit.RMSErr = parerr[index + 2];
      onehit.ChiDOF = chidof;
      onehit.Wire = theWire;
      onehit.WireNum = theWireNum;
      onehit.numHits = nhits;
      onehit.LoHitID = lohitid;

  if(prt) {
    std::cout<<"W:H "<<theWireNum;
    std::cout<<":"<<allhits.size()<<" Chg "<<(short)onehit.Charge;
    std::cout<<" Time "<<(short)onehit.Time<<" RMS "<<onehit.RMS;
    std::cout<<" lo/hi ID "<<onehit.LoHitID;
    std::cout<<" chidof "<<chidof;
    std::cout<<std::endl;
  }

      allhits.push_back(onehit);
    } // hit
  } // StoreHits


} // namespace cluster

