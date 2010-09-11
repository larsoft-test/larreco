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
#include <sstream>
#include <math.h>
#include <algorithm>

// Framework includes
//////////////////////////////////////////////////////
// Of course these are all different
//////////////////////////////////////////////////////
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/ParameterSet/interface/ParameterSet.h" 
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h" 
#include "DataFormats/Common/interface/Handle.h" 
#include "DataFormats/Common/interface/Ptr.h" 
#include "DataFormats/Common/interface/PtrVector.h" 
#include "FWCore/Framework/interface/MakerMacros.h" 
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "FWCore/Services/interface/TFileService.h" 
#include "FWCore/Framework/interface/TFileDirectory.h" 
#include "FWCore/MessageLogger/interface/MessageLogger.h" 
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h" 

#include "HitFinder/inc/FFTHitFinder.h"
#include "Geometry/inc/geo.h"
#include "Geometry//inc/WireGeo.h"
#include "RecoBase/inc/Hit.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TDecompSVD.h"

namespace hit{

  DEFINE_FWK_MODULE(FFTHitFinder);


//-------------------------------------------------


  FFTHitFinder::FFTHitFinder(edm::ParameterSet const& pset) :
    fSpacer (false),
    fCalDataModuleLabel(pset.getParameter< std::string  >("CalDataModuleLabel")),
    fMinSigInd(pset.getParameter< double >("MinSigInd")),
    fMinSigCol(pset.getParameter< double >("MinSigCol")),
    fIndWidth(pset.getParameter< double >("IndWidth")),
    fColWidth(pset.getParameter< double >("ColWidth")),
    fDrift(pset.getParameter< double >("Drift")),
    fPOffset(pset.getParameter< double >("POffset")),
    fOOffset(pset.getParameter< double >("OOffset")),
    fMaxMultiHit(pset.getParameter< int >("MaxMultiHit"))
  {
    produces< std::vector<recob::Hit> >();
  }


//-------------------------------------------------
  FFTHitFinder::~FFTHitFinder()
{
}

//-------------------------------------------------
void FFTHitFinder::beginJob(const edm::EventSetup&)
{
}
void FFTHitFinder::endJob()
{
}
//  This algorithm uses the fact that deconvoluted signals are very smooth 
//  and looks for hits 
//  as areas between local minima that have signal above threshold.
//-------------------------------------------------
  void FFTHitFinder::produce(edm::Event& evt, edm::EventSetup const&)
{ 

  //////////////////////////////////////////////////////
  // Make a std::auto_ptr<> for the thing you want to put into the event
  // because that handles the memory management for you
  //////////////////////////////////////////////////////
  std::auto_ptr<std::vector<recob::Hit> > hcol(new std::vector<recob::Hit>);
 // Read in the wire List object(s).
  edm::Handle< std::vector<recob::Wire> > wireVecHandle;
  evt.getByLabel(fCalDataModuleLabel,wireVecHandle);

  edm::Service<geo::Geometry> geom;
  edm::Service<edm::TFileService> tfs;
	
  std::vector<double> signal;            
  std::vector<int> startTimes;  //stores time of 1st local minimum
  std::vector<int> maxTimes;    //stores time of local maximum
  std::vector<int> endTimes;    //stores time of 2nd local minimum
  //  std::vector<recob::Hit *> hitVector;  //holds hits to be saved
  std::vector<double>::iterator timeIter;  //iterator for time bins
  int minTimeHolder;  //hold position of minTime for hit region
  int time;  //current time bin
  unsigned int channel;  //channel number
  unsigned int wire;              //wire number
  unsigned int plane;             //plane number
  bool maxFound; //Flag for whether a value>threshold has been found
  double threshold;
  std::string eqn;
  std::stringstream numConv;
  TF1 * gSum;
  TH1D * hitSignal;
  geo::SigType_t sigType;
  geo::View_t view;     
  //loop over wires

  for(unsigned int wireIter=0; wireIter<wireVecHandle->size(); wireIter++) {
    edm::Ptr<recob::Wire> wireVec(wireVecHandle, wireIter);
    recob::Wire* wireVec2 = new recob::Wire();
    *wireVec2=*wireVec;
    startTimes.clear();
    maxTimes.clear();
    endTimes.clear();
    signal = wireVec->fSignal;
    time=0;
    minTimeHolder = 0;
    maxFound = false;
    channel=wireVec->RawDigit()->Channel();
    //std::cout << "Channel: " << channel << std::endl;
    geom->ChannelToWire(channel,plane,wire);
    sigType = geom->Plane(plane).SignalType();
    view = geom->Plane(plane).View();
    if(sigType == geo::kInduction) threshold=fMinSigInd;
    else threshold= fMinSigCol;
    //loop over signal
    for(timeIter = signal.begin();timeIter+2!=signal.end();timeIter++) 
    {    
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
    }
    //if no inflection found before end, add end point
    if(maxFound&&maxTimes.size()>endTimes.size()) endTimes.push_back(signal.size()-1); 	  
    //All code below does the fitting, adding of hits
    //to the hit vector and when all wires are complete 
    //saving them 
    //std::cout << startTimes.size() <<" "<<maxTimes.size() << " "<< endTimes.size()<<std::endl;
    //for(int i = 0; i < maxTimes.size(); i++) std::cout << "times: " <<startTimes[i] << " " << maxTimes[i] << " " << endTimes[i] << std::endl;
    double totSig(0); //stoes the total hit signal
    double startT(0); //stores the start time
    double endT(0);  //stores the end time
    int numHits;
    int size(0); 
    int hitIndex(0);
    double amplitude, position,width;
    double XYZ[3]={0,0,0};
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
      while(numHits < fMaxMultiHit && numHits+hitIndex < (signed)endTimes.size() && signal[endTimes[hitIndex+numHits-1]] >threshold/2.0 && startTimes[hitIndex+numHits] - endTimes[hitIndex+numHits-1] < 2) numHits++;
      /*std::cout <<"Hit times: "<< startTimes[hitIndex] << " " 
		<< maxTimes[hitIndex] << " "
		<< endTimes[hitIndex] << " "
		<< startTimes[hitIndex+1] << " "
		<< maxTimes[hitIndex+1] << " "
		<< endTimes[hitIndex+1] << " # hits: "<< numHits << std::endl;
      */	  
      //finds the first point >0
      startT=startTimes[hitIndex];
      while(signal[(int)startT] < 0) startT++;
      //finds the first point from the end >0
      endT=endTimes[hitIndex+numHits-1];
      while(signal[(int)endT] <0) endT--;
      size = (int)(endT-startT);
      hitSignal = tfs->make<TH1D>("hitSignal","",size,startT,endT);
      for(int i = (int)startT; i < (int)endT; i++) hitSignal->Fill(i,signal[i]);
      //std::cout << wireVec->RawDigit()->Channel() << " size " << size << std::endl;
      for(int i = 3; i < numHits*3; i+=3) {
        eqn.append("+gaus(");
        numConv.str("");
        numConv << i;
        eqn.append(numConv.str());
        eqn.append(")");
      }
      gSum = tfs->make<TF1>("gSum",eqn.c_str(),0,size);
      if(numHits > 1) {
	TArrayD data(numHits*numHits);
	TVectorD amps(numHits); 
	for(int i = 0; i < numHits; i++) {
	  amps[i]=signal[maxTimes[hitIndex+i]];
	  //std::cout <<" ai: " <<amps[i] ;
	  for(int j = 0; j < numHits;j++) 
	    data[i+numHits*j]=TMath::Gaus(maxTimes[hitIndex+j],maxTimes[hitIndex+i],width);
	}
        
	//std::cout <<"channel"<<channel<< std::endl;
	  TMatrixD h(numHits,numHits);
	  h.Use(numHits,numHits,data.GetArray());
	  TDecompSVD a(h);
	  a.Solve(amps);
         
	for(int i = 0;i < numHits; i++) {
	  //std::cout <<" af: " << amps[i] <<" "<< maxTimes[hitIndex+i]<<" "<<width<<std::endl;
	  gSum->SetParameter(3*i, amps[i]);
	  gSum->SetParameter(1+3*i, maxTimes[hitIndex+i]);
	  gSum->SetParameter(2+3*i, width);
          gSum->SetParLimits(3*i, 0.0, 10.0*amps[i]);
          gSum->SetParLimits(1+3*i, startT , endT);
          gSum->SetParLimits(2+3*i, 0.0, 10.0*width);
	}
      }
      else {
        gSum->SetParameters(signal[maxTimes[hitIndex]],maxTimes[hitIndex],width);
        gSum->SetParLimits(0,0.0,1.5*signal[maxTimes[hitIndex]]);
        gSum->SetParLimits(1, startT , endT);
        gSum->SetParLimits(2,0.0,10.0*width);
        //std::cout <<" af: " << signal[maxTimes[hitIndex]] <<" "<< maxTimes[hitIndex]<<" "<<width<<std::endl;
      }
      hitSignal->Fit(gSum,"QNR","", startT, endT);
      for(int hitNumber = 0; hitNumber < numHits; hitNumber++) {
	if(gSum->GetParameter(3*hitNumber) > threshold/2.0) { 
	  // Make this an object, not ptr in FMWK->ART port. EC, 10-Sep-2010.
	  recob::Hit hit(wireVec);// = new recob::Hit(wireVec); 
	  //recob::Hit hit(wireVec2);
          amplitude = gSum->GetParameter(3*hitNumber);
          position = gSum->GetParameter(3*hitNumber+1);
          width = gSum->GetParameter(3*hitNumber+2);
          hitSig.clear();
          hitSig.resize(size);
	  for(double sigPos = 0; sigPos<size; sigPos++){
	    hitSig[(int)sigPos] = amplitude*TMath::Gaus(sigPos,position, width);
	    totSig+=hitSig[(int)sigPos]; 
	  }              	    
	  hit.fHitSignal=hitSig;            
	  hit.SetCrossingTime(position);             
	  hit.SetUpTime(2*width);
	  hit.SetStartTime(position-width);
	  hit.SetEndTime(position+width);
	  hit.SetUpADC(totSig);
	  hit.SetMIPs(amplitude);
	  XYZ[0]=XYZ[1]=XYZ[2]=0;
	  if(sigType == geo::kInduction) XYZ[0] = fPOffset*fDrift;
	  XYZ[0]+=fDrift*(position) - fOOffset;
	  hit.SetXYZ(XYZ);
	  hit.SetView(view);

	  hcol->push_back(hit);

	}
      }
      delete hitSignal;
      delete gSum;
      hitIndex+=numHits;


    }  // end while on hitIndex<(signed)startTimes.size()
  }    // while on Wires

      evt.put(hcol);

}      // End of produce()
  



} // end of hit namespace
