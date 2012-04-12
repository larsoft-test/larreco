#ifndef PROPTEST2_H
#define PROPTEST2_H

//
// Name:  PropTest.h
//
// Purpose: PropTest module.  Test propagators.
//
// Created:  12-Apr-2012  H. Greenlee

#include "art/Framework/Core/EDAnalyzer.h"

namespace trkf
{
  class PropTest : public art::EDAnalyzer
  {
  public:

    // Constructor, destructor.

    explicit PropTest(fhicl::ParameterSet const& pset);
    ~PropTest();

    // Overrides.

    void beginJob();
    void analyze(const art::Event& evt);

  };
}

#endif
