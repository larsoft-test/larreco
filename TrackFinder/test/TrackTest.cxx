//
// File: TrackTest.cxx
//
// Purpose: Single source executable with tests for track classes 
//

#include <iostream>
#include <cassert>
#include "TrackFinder/KTrack.h"
#include "TrackFinder/KETrack.h"
#include "TrackFinder/KFitTrack.h"

int main()
{
  // Make sure assert is enabled.

  bool assert_flag = false;
  assert((assert_flag = true, assert_flag));
  if ( ! assert_flag ) {
    std::cerr << "Assert is disabled" << std::endl;
    return 1;
  }

  // Make some tracks.

  trkf::KTrack trk;
  trkf::KETrack tre;
  trkf::KFitTrack trf;

  // Some simple tests.

  assert(!trk.isValid());
  assert(trf.stat() == trkf::KFitTrack::INVALID);

  // Done (success).

  std::cout << "TrackTest: All tests passed." << std::endl;

  return 0;
}
