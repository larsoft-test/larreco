////////////////////////////////////////////////////////////////////////
// Class:       HitCheater
// Module Type: producer
// File:        HitCheater.cxx
//
// Generated at Tue Nov  8 09:41:20 2011 by Brian Rebel using artmod
// from art v1_00_02.
////////////////////////////////////////////////////////////////////////

#include "HitFinder/HitCheater.h"
#include "Geometry/geo.h"
#include "Simulation/sim.h"
#include "RecoBase/recobase.h"
#include "Utilities/DetectorProperties.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Event.h"

#include "TH1.h"
#include "TString.h"

//-------------------------------------------------------------------
hit::HitCheater::HitCheater(fhicl::ParameterSet const & p)
{
  this->reconfigure(p);

  produces< std::vector<recob::Hit> >();
}

//-------------------------------------------------------------------
hit::HitCheater::~HitCheater()
{
}

//-------------------------------------------------------------------
void hit::HitCheater::produce(art::Event & e)
{
  // make the auto_ptr for the hits
  std::auto_ptr< std::vector<recob::Hit> > hits(new std::vector<recob::Hit>);

  // Read in the wire List object(s).
  art::Handle< std::vector<recob::Wire> > wHandle;
  e.getByLabel(fWireModuleLabel,wHandle);
  
  // make a map of wires to channel numbers
  std::map<unsigned int, art::Ptr<recob::Wire> > wireMap;
  
  for(size_t wc = 0; wc < wHandle->size(); ++wc){
    art::Ptr<recob::Wire> wire(wHandle, wc);
    wireMap[wire->RawDigit()->Channel()] = wire;
  }

  // get the sim::SimChannels out of the event
  std::vector<const sim::SimChannel*> sccol;
  e.getView(fG4ModuleLabel, sccol);

  for(size_t sc = 0; sc < sccol.size(); ++sc){

    // loop over all the ides for this channel
    const std::map<unsigned short, std::vector<sim::IDE> >& idemap = sccol[sc]->TDCIDEMap();
    
    FindHitsOnChannel(idemap, *hits, wireMap.find(sccol[sc]->Channel())->second);

  }// end loop over SimChannels

  // put the cheated hits into the event
  LOG_DEBUG("HitCheater") << "putting " << hits->size() << " hits into the event";
  e.put(hits);

  return;
}

//-------------------------------------------------------------------
void hit::HitCheater::FindHitsOnChannel(std::map<unsigned short, std::vector<sim::IDE> > const& idemap,
					std::vector<recob::Hit>& hits,
					art::Ptr<recob::Wire>&   wire)
{
  // make a map of the total signal at each tdc value
  // if there is a gap of >= 1 tdc count between signals
  // make separate hits.  If there is a dip between 
  // two peaks of order 600 electrons, make multiple hits
  std::map<unsigned short, double> tdcVsE;
  std::vector<unsigned short> tdcEnds;

  std::map<unsigned short, std::vector<sim::IDE> >::const_iterator mapitr = idemap.begin();
  unsigned short prev = mapitr->first;
  for(mapitr = idemap.begin(); mapitr != idemap.end(); mapitr++){
    unsigned short tdc = mapitr->first;

    // more than a one tdc gap between times with 
    // signal, start a new hit
    if(tdc - prev > 1) tdcEnds.push_back(prev);

    double total = 0.;
    // loop over the ides
    std::vector<sim::IDE>::const_iterator ideitr = mapitr->second.begin();
    for(ideitr = mapitr->second.begin(); ideitr != mapitr->second.end(); ideitr++){
      total += ideitr->numElectrons;
    }
    
    tdcVsE[tdc] = total;

    // reset the starts variable
    prev = tdc;
  }// end loop over the ide map
  tdcEnds.push_back(prev);

  // now loop over the tdcVsE map and make some hits
  std::map<unsigned short, double>::iterator tdcVsEitr;  
  std::vector<unsigned short>::iterator endItr = tdcEnds.begin();
  for(tdcVsEitr = tdcVsE.begin(); tdcVsEitr != tdcVsE.end() && endItr != tdcEnds.end(); tdcVsEitr++){
    double startTime =  1.*tdcVsEitr->first;
    double endTime   =  0.;
    double peakTime  =  0.;
    double maxCharge = -1.;
    double totCharge =  0.;
    int multiplicity =  1 ;
    while(tdcVsEitr->first <  *endItr){
      double adc = fElectronsToADC*tdcVsEitr->first;
      totCharge += adc;
      if(adc > maxCharge){
	maxCharge = adc;
	peakTime = 1.*tdcVsEitr->first;
      }
      endTime = 1.*tdcVsEitr->first;
      tdcVsEitr++;
      if(tdcVsEitr == tdcVsE.end()) break;
    } // end loop for current hit

    // make the new hit object
    // use 1. as the uncertainty on the times, as we know them to within 1 tdc 
    // tick when cheating, make the uncertainty on the charges just be the sqrt 
    // of their values.
    // don't make an object if the total charge is less than 5 adc count
    if(totCharge > 5.){
      hits.push_back(recob::Hit(wire, 
				startTime, 1.,
				endTime,   1.,
				peakTime,  1.,
				totCharge, sqrt(totCharge),
				maxCharge, sqrt(maxCharge),
				multiplicity,
				1.)
		     );
      
      LOG_DEBUG("HitCheater") << "new hit is " << hits.back();
    }

    endItr++;
  }// end loop over tdc values

  return;
}

//-------------------------------------------------------------------
void hit::HitCheater::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  // get the geometry to determine how many channels
  // are in the induction and collection planes...
  // only worry about 1 induction plane for now
  // and the first tpc
  art::ServiceHandle<geo::Geometry> geo;

  fChannelEnergyDepsInd.clear();
  fChannelEnergyDepsCol.clear();

  // plane 0 is always induction, the last plane is always collection,
  // just worry about tpc 0.
  fChannelEnergyDepsInd.resize(geo->TPC(0).Plane(0).Nwires());
  for(unsigned int w = 0; w < geo->TPC(0).Plane(0).Nwires(); ++w){
    TString name = "InductionChannel_";
    name += w;
    fChannelEnergyDepsInd[w] = tfs->make<TH1D>(name, ";TDC;Electrons;", 4096, 0., 4095);
  }

  fChannelEnergyDepsCol.resize(geo->TPC(0).Plane(geo->TPC(0).Nplanes()-1).Nwires());
  for(unsigned int w = 0; w < geo->TPC(0).Plane(geo->TPC(0).Nplanes()-1).Nwires(); ++w){
    TString name = "CollectionChannel_";
    name += w;
    fChannelEnergyDepsCol[w] = tfs->make<TH1D>(name, ";TDC;Electrons;", 4096, 0., 4095);
  }
  
  return;
}

//-------------------------------------------------------------------
void hit::HitCheater::reconfigure(fhicl::ParameterSet const & p)
{
  fG4ModuleLabel   = p.get< std::string >("G4ModuleLabel",   "largeant");
  fWireModuleLabel = p.get< std::string >("WireModuleLabel", "caldata" );

  art::ServiceHandle<util::DetectorProperties> detprop;
  fElectronsToADC = detprop->ElectronsToADC();

  return;
}

