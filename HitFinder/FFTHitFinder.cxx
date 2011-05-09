////////////////////////////////////////////////////////////////////////
//
// FFTHitFinder class
//
// pagebri3@msu.edu
//
//  This algorithm is designed to find hits on wires after deconvolution
//  with an average shape used as the input response.
////////////////////////////////////////////////////////////////////////

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

// Framework includes
#include "art/Framework/Core/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Core/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft Includes
#include "HitFinder/FFTHitFinder.h"
#include "Geometry/geo.h"
#include "RecoBase/recobase.h"

// ROOT 
#include "TH1D.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TF1.h"

namespace hit{

  //-------------------------------------------------
  FFTHitFinder::FFTHitFinder(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
    produces< std::vector<recob::Hit> >();
  }


  //-------------------------------------------------
  FFTHitFinder::~FFTHitFinder()
  {
  }
  
  //-------------------------------------------------
  void FFTHitFinder::reconfigure(fhicl::ParameterSet p)
  {
    fCalDataModuleLabel = p.get< std::string  >("CalDataModuleLabel");
    fMinSigInd          = p.get< double       >("MinSigInd");
    fMinSigCol          = p.get< double       >("MinSigCol"); 
    fIndWidth           = p.get< double       >("IndWidth");  
    fColWidth           = p.get< double       >("ColWidth");  	  	  
    fMaxMultiHit        = p.get< int          >("MaxMultiHit");
  }  
  //-------------------------------------------------
  void FFTHitFinder::beginJob()
  {
  }
  void FFTHitFinder::endJob()
  {
  }

  //  This algorithm uses the fact that deconvoluted signals are very smooth 
  //  and looks for hits 
  //  as areas between local minima that have signal above threshold.
  //-------------------------------------------------
  void FFTHitFinder::produce(art::Event& evt)
  { 
    
    std::auto_ptr<std::vector<recob::Hit> > hcol(new std::vector<recob::Hit>);
    // Read in the wire List object(s).
    art::Handle< std::vector<recob::Wire> > wireVecHandle;
    evt.getByLabel(fCalDataModuleLabel,wireVecHandle);
    
    art::ServiceHandle<geo::Geometry> geom;
    
    std::vector<double> signal;            
    std::vector<int> startTimes;             // stores time of 1st local minimum
    std::vector<int> maxTimes;    	     // stores time of local maximum    
    std::vector<int> endTimes;    	     // stores time of 2nd local minimum
    std::vector<double>::iterator timeIter;  // iterator for time bins
    int minTimeHolder    = 0;                // hold position of minTime for hit region
    int time             = 0;                // current time bin
    unsigned int channel = 0;                // channel number
    unsigned int wire    = 0;                // wire number
    unsigned int plane   = 0;                // plane number
    bool maxFound        = false;            // Flag for whether a value > threshold 
                                             // has been found
    double threshold     = 0.;               // minimum signal size for id'ing a hit
    std::string eqn      = "gaus(0)";        // string for equation to fit
    std::stringstream numConv;
    geo::SigType_t sigType = geo::kInduction;// type of plane we are looking at

    //loop over wires

    for(unsigned int wireIter=0; wireIter<wireVecHandle->size(); wireIter++) {
      art::Ptr<recob::Wire> wireVec(wireVecHandle, wireIter);
      startTimes.clear();
      maxTimes.clear();
      endTimes.clear();
      signal        = wireVec->fSignal;
      time          = 0;
      minTimeHolder = 0;
      maxFound      = false;
      channel       = wireVec->RawDigit()->Channel();
      geom->ChannelToWire(channel,plane,wire);
      sigType       = geom->Plane(plane).SignalType();
      threshold     = fMinSigInd;
      if(sigType == geo::kCollection) threshold = fMinSigCol;

      // loop over signal
      for(timeIter = signal.begin();timeIter+2!=signal.end();timeIter++) {    
	//test if timeIter2 is a local minimum
	if(*timeIter>*(timeIter+1) && *(timeIter+1)<*(timeIter+2)) {
	  //only add points if we've already found a local max above threshold.
	  if(maxFound) {
	    endTimes.push_back(time+1);
	    maxFound = false;
	    //keep these in case new hit starts right away
	    minTimeHolder = time+2; 
	  }
	  else minTimeHolder = time+1; 
	}
	//if not a minimum test if we are at a local maximum 
	//if so and the max value is above threshold add it and proceed.
	else if(*timeIter<*(timeIter+1) && 
		*(timeIter+1)>*(timeIter+2) && 
		*(timeIter+1) > threshold) { 
	  maxFound = true;
	  maxTimes.push_back(time+1);
	  startTimes.push_back(minTimeHolder);          
	}
	time++;
      }//end loop over signal vec

      //if no inflection found before end, add end point
      if(maxFound&&maxTimes.size()>endTimes.size()) 
	endTimes.push_back(signal.size()-1); 	  
    
      //All code below does the fitting, adding of hits
      //to the hit vector and when all wires are complete 
      //saving them 
      mf::LogDebug("FFTHitFinder") <<channel<<" "<< startTimes.size() <<" "<<maxTimes.size() 
      << " "<< endTimes.size();

      double totSig(0); //stoes the total hit signal
      double startT(0); //stores the start time
      double endT(0);  //stores the end time
      int numHits(0);  //number of consecutive hits being fitted
      int size(0);     //size of data vector for fit
      int hitIndex(0);  //index of current hit in sequence
      double amplitude(0), position(0), width(0);  //fit parameters
      double amplitudeErr(0), positionErr(0), widthErr(0);  //fit errors
      double goodnessOfFit(0), chargeErr(0);  //Chi2/NDF and error on charge
     
      //stores gaussian paramters first index is the hit number
      //the second refers to height, position, and width respectively
      std::vector<double>  hitSig;
      //add found hits to hit vector
      while(hitIndex<(signed)startTimes.size()) {
	eqn="gaus(0)";
	if(sigType==geo::kInduction) width=fIndWidth;
	else width= fColWidth;
	startT=endT=0;
	numHits=1;
	while(numHits < fMaxMultiHit 
	      && numHits+hitIndex < (signed)endTimes.size() 
	      && signal[endTimes[hitIndex+numHits-1]] >threshold/2.0 
	      && startTimes[hitIndex+numHits] - endTimes[hitIndex+numHits-1] < 2) 
	  numHits++;

	//finds the first point >0
	startT=startTimes[hitIndex];
	while(signal[(int)startT] < 0) startT++;
	//finds the first point from the end >0
	endT=endTimes[hitIndex+numHits-1];
	while(signal[(int)endT] <0) endT--;
	size = (int)(endT-startT);
	TH1D hitSignal("hitSignal","",size,startT,endT);

	for(int i = (int)startT; i < (int)endT; i++)
	  hitSignal.Fill(i,signal[i]);
	

	for(int i = 3; i < numHits*3; i+=3) {
	  eqn.append("+gaus(");
	  numConv.str("");
	  numConv << i;
	  eqn.append(numConv.str());
	  eqn.append(")");
	}
	TF1 gSum("gSum",eqn.c_str(),0,size);
	if(numHits > 1) {
	  TArrayD data(numHits*numHits);
	  TVectorD amps(numHits); 
	  for(int i = 0; i < numHits; i++) {
	    amps[i]=signal[maxTimes[hitIndex+i]];
	    mf::LogDebug("FFTHitFinder") <<" ai: " <<amps[i] ;
	    for(int j = 0; j < numHits;j++) 
	      data[i+numHits*j] = TMath::Gaus(maxTimes[hitIndex+j],
					      maxTimes[hitIndex+i],
					      width);
	  }//end loop over hits
      //This section uses a linear approximation in order to get an
      //initial value of the individual hit amplitudes 
      try
      {
	  TMatrixD h(numHits,numHits);
	  h.Use(numHits,numHits,data.GetArray());
	  TDecompSVD a(h);
	  a.Solve(amps);
      }
      catch(...){mf::LogInfo("FFTHitFinder")<<"TDcompSVD failed";hitIndex+=numHits;continue;}
      
	  for(int i = 0;i < numHits; i++) {
	    gSum.SetParameter(3*i, amps[i]);
	    gSum.SetParameter(1+3*i, maxTimes[hitIndex+i]);
	    gSum.SetParameter(2+3*i, width);
	    gSum.SetParLimits(3*i, 0.0, 10.0*amps[i]);
	    gSum.SetParLimits(1+3*i, startT , endT);
	    gSum.SetParLimits(2+3*i, 0.0, 10.0*width);
	  }//end loop over hits
	}//end if numHits > 1
	else {
	  gSum.SetParameters(signal[maxTimes[hitIndex]],maxTimes[hitIndex],width);
	  gSum.SetParLimits(0,0.0,1.5*signal[maxTimes[hitIndex]]);
	  gSum.SetParLimits(1, startT , endT);
	  gSum.SetParLimits(2,0.0,10.0*width);
	}

	/// \todo - just get the integral from the fit for totSig
	hitSignal.Fit(&gSum,"WQNR","", startT, endT);
	for(int hitNumber = 0; hitNumber < numHits; hitNumber++) {
	  if(gSum.GetParameter(3*hitNumber) > threshold/2.0) { 
	    amplitude = gSum.GetParameter(3*hitNumber);
	    position = gSum.GetParameter(3*hitNumber+1);
	    width = gSum.GetParameter(3*hitNumber+2);
            amplitudeErr = gSum.GetParError(3*hitNumber);
	    positionErr = gSum.GetParError(3*hitNumber+1);
	    widthErr = gSum.GetParError(3*hitNumber+2);
            goodnessOfFit= gSum.GetChisquare()/(double)gSum.GetNDF();
            chargeErr=TMath::Sqrt(TMath::Pi())*(amplitudeErr*width+widthErr*amplitude);   //estimate from area of Gaussian
	    hitSig.resize(size);
	    for(int sigPos = 0; sigPos<size; sigPos++){
	      hitSig[sigPos] = amplitude*TMath::Gaus(sigPos+startT,position, width);
	      totSig+=hitSig[(int)sigPos];
              
	    }              	    

	    // make the hit
	    recob::Hit hit(wireVec, 
			   position-width, 
                           widthErr,
			   position+width, 
                           widthErr,
			   position,
                           positionErr,
			   totSig,         
                           chargeErr,                
			   amplitude,
                           amplitudeErr,
			   1,                  /// \todo - mulitplicity has to be determined
			   goodnessOfFit);               	    
	    hcol->push_back(hit);
	    
	  }//end if over threshold
	}//end loop over hits
	hitIndex+=numHits;

      } // end while on hitIndex<(signed)startTimes.size()
    } // while on Wires

    evt.put(hcol);

  } // End of produce()
  
} // end of hit namespace
