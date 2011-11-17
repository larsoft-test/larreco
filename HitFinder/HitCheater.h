#ifndef HitCheater_h
#define HitCheater_h
////////////////////////////////////////////////////////////////////////
// Class:       HitCheater
// Module Type: producer
// File:        HitCheater.h
//
// Generated at Tue Nov  8 09:41:20 2011 by Brian Rebel using artmod
// from art v1_00_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"

class TH1D;

namespace hit {
  class HitCheater;
}

namespace recob {
  class Wire;
  class Hit;
}

namespace sim {
  class IDE;
}

class hit::HitCheater : public art::EDProducer {
public:
  explicit HitCheater(fhicl::ParameterSet const & p);
  virtual ~HitCheater();

  virtual void produce(art::Event & e);

  virtual void beginJob();
  virtual void reconfigure(fhicl::ParameterSet const & p);

private:

  void FindHitsOnChannel(std::map<unsigned short, std::vector<sim::IDE> > const& idemap,
			 std::vector<recob::Hit>& hits,
			 art::Ptr<recob::Wire>&   wire);

  std::string         fG4ModuleLabel;        ///< label name for module making sim::SimChannels
  std::string         fWireModuleLabel;      ///< label name for module making recob::Wires
  std::vector<TH1D *> fChannelEnergyDepsInd; ///< Energy depositions vs time for each induction channel
  std::vector<TH1D *> fChannelEnergyDepsCol; ///< Energy depositions vs time for each collection channel

  double              fElectronsToADC;       ///< Conversion factor of electrons to ADC counts

};
#endif /* HitCheater_h */
