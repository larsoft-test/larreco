///////////////////////////////////////////////////////////////////////
///
/// \file   KETrack.cxx
///
/// \brief  Basic Kalman filter track class, with error.
///
/// \author H. Greenlee
///
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "TrackFinder/KETrack.h"

namespace trkf {

  /// Default constructor.
  KETrack::KETrack()
  {}

  /// Constructor - specify surface only.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  ///
  KETrack::KETrack(const boost::shared_ptr<const Surface>& psurf) :
    KTrack(psurf)
  {}

  /// Constructor - surface + track parameters + error matrix.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  /// vec   - Track state vector.
  /// err   - Track error matrix.
  /// dir   - Track direction.
  /// pdg   - Pdg code.
  ///
  KETrack::KETrack(const boost::shared_ptr<const Surface>& psurf,
		   const TrackVector& vec,
		   const TrackError& err, 
		   Surface::TrackDirection dir,
		   int pdg) :
    KTrack(psurf, vec, dir, pdg),
    fErr(err)
  {}

  /// Constructor - KTrack + error matrix.
  ///
  /// Arguments:
  ///
  /// trk - KTrack.
  /// err - Track error matrix.
  ///
  KETrack::KETrack(const KTrack& trk, const TrackError& err) :
    KTrack(trk),
    fErr(err)
  {}

  /// Destructor.
  KETrack::~KETrack()
  {}

  /// Printout
  std::ostream& KETrack::Print(std::ostream& out, bool doTitle) const
  {
    if(doTitle)
      out << "KETrack:\n";

    // Print base class.

    KTrack::Print(out, false);

    // Print diagonal errors.

    out << "  Diagonal errors:\n"
	<< "  [";
    for(unsigned int i = 0; i < fErr.size1(); ++i) {
      if(i != 0)
	out << ", ";
      double err = fErr(i,i);
      err = (err >= 0. ? std::sqrt(err) : -std::sqrt(-err));
      out << err;
    }
    out << "]\n";

    // Print correlations.

    out << "  Correlation matrix:";
    for(unsigned int i = 0; i < fErr.size1(); ++i) {
      if(i == 0)
	out << "\n  [";
      else
	out << "\n   ";
      for(unsigned int j = 0; j <= i; ++j) {
	if(j != 0)
	  out << ", ";
	if(i == j)
	  out << 1.;
	else {
	  double eiijj = fErr(i,i) * fErr(j,j);
	  double eij = fErr(i,j);
	  if(eiijj != 0.)
	    eij /= std::sqrt(std::abs(eiijj));
	  else
	    eij = 0.;
	  out << eij;
	}
      }
    }
    out << "]\n";
    return out;
  }

} // end namespace trkf
