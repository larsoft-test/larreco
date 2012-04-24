#ifndef KALMANFILTERTEST2_H
#define KALMANFILTERTEST2_H

//
// Name:  KalmanFilterTest.h
//
// Purpose: KalmanFilterTest module.
//
// Created:  12-Apr-2012  H. Greenlee

#include "art/Framework/Core/EDAnalyzer.h"
#include "TrackFinder/KalmanFilterAlg.h"

namespace trkf
{
  class KalmanFilterTest : public art::EDAnalyzer
  {
  public:

    // Constructor, destructor.

    explicit KalmanFilterTest(fhicl::ParameterSet const& pset);
    ~KalmanFilterTest();

    // Overrides.

    void beginJob();
    void analyze(const art::Event& evt);

  private:

    // Attributes.

    trkf::KalmanFilterAlg fKFAlg;

  };
}

#endif
