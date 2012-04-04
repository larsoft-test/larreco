////////////////////////////////////////////////////////////////////////
///
/// \file   KHit.h
///
/// \brief  Kalman filter measurement class template.
///
/// \author H. Greenlee 
///
/// Class KHit represents a general measurement on a surface.  It is
/// specialized compared to base class KHitBase by specifying the
/// dimension of the measurement vector by an integer template
/// parameter N.
///
/// KHit class inherits the following attribute from KHitBase.
///
/// 1.  Measurement surface.
///
/// KHit adds the following attributes.
///
/// 2.  Measurement vector.
/// 3.  Measurement error matrix.
/// 4.  Prediction vector.
/// 5.  Prediction error matrix.
/// 6.  Residual vector.
/// 7.  Residual error matrix.
/// 8.  Inverse of residual error matrix.
/// 9.  Kalman H-matrix.
/// 10.  Incremental chisquare.
///
/// The first two attributes (measurement vector + error matrix) are
/// filled during construction, and the remaining attributes are left
/// in a default state.
///
/// The remaining attributes are filled by the prediction method.  The
/// attributes that are filled in by the prediction method are
/// mutable, and prediction method is const.  The actual calculation
/// of the prediction vector, prediction error matrix, and H-matrix
/// are left to be implemented by a derived class, which must override
/// the pure virtual method subpredict.  The remaining unfilled
/// attributes are calculated locally in this class in the prediction
/// method inherited from KHitBase.
///
/// The intended use case is as follows.
///
/// 1.  Track (KETrack) is propagated to measurement surface.
///     Methods predict and update will throw an exception of track
///     surface does not match measurement surface.
/// 2.  Prediction is updated by calling method KHit::predict.
/// 3.  At this point the calling program can make a cut on the
///     incremental chisquare, returned by method KHit::chisq.
/// 4.  If the chisquare cut passes, update track by calling method
///     KHit::update.
///
////////////////////////////////////////////////////////////////////////

#ifndef KHIT_H
#define KHIT_H

#include "TrackFinder/KHitBase.h"
#include "cetlib/exception.h"

namespace trkf {

  template <int N>
  class KHit : public KHitBase
  {
  public:

    /// Default constructor.
    KHit();

    /// Initializing Constructor -- surface only.
    KHit(const boost::shared_ptr<const Surface>& psurf);

    /// Fully Initializing Constructor.
    KHit(const boost::shared_ptr<const Surface>& psurf,
	 const typename KVector<N>::type& mvec,
	 const typename KSymMatrix<N>::type& merr);

    /// Destructor.
    virtual ~KHit();

    // Modifiers.

    /// Set measurement vector.
    void setMeasVector(const typename KVector<N>::type& mvec) {fMvec = mvec;}

    /// Set measurement error.
    void setMeasError(const typename KSymMatrix<N>::type& merr) {fMerr = merr;}

    // Accessors.

    /// Measurement vector.
    const typename KVector<N>::type& getMeasVector() const {return fMvec;}

    /// Measurement error matrix.
    const typename KSymMatrix<N>::type& getMeasError() const {return fMerr;}

    /// Prediction vector.
    const typename KVector<N>::type& getPredVector() const {return fPvec;}

    /// Prediction matrix.
    const typename KSymMatrix<N>::type& getPredError() const {return fPerr;}

    /// Residual vector.
    const typename KVector<N>::type& getResVector() const {return fRvec;}

    /// Residual error matrix.
    const typename KSymMatrix<N>::type& getResError() const {return fRerr;}

    /// Residual inv. error matrix.
    const typename KSymMatrix<N>::type& getResInvError() const {return fRinv;}

    /// Kalman H-matrix.
    const typename KHMatrix<N>::type& getH() const {return fH;}

    /// Incremental chisquare.
    double getChisq() const {return fChisq;}

    // Overrides.
    // Implementation of overrides is found at the bottom of this header.

    /// Prediction method (return false if fail).
    bool predict(const KETrack& tre) const;

    /// Update track method.
    void update(KETrack& tre) const;

    // Pure virtual methods.

    /// Calculate prediction function (return via arguments).
    virtual bool subpredict(const KETrack& tre,
			    typename KVector<N>::type& pvec,
			    typename KSymMatrix<N>::type& perr,
			    typename KHMatrix<N>::type& hmatrix) const = 0;

  private:

    // Attributes.

    typename KVector<N>::type fMvec;             ///< Measurement vector.
    typename KSymMatrix<N>::type fMerr;          ///< Measurement error matrix.
    mutable typename KVector<N>::type fPvec;     ///< Prediction vector.
    mutable typename KSymMatrix<N>::type fPerr;  ///< Prediction  error matrix.
    mutable typename KVector<N>::type fRvec;     ///< Residual vector.
    mutable typename KSymMatrix<N>::type fRerr;  ///< Residual error matrix.
    mutable typename KSymMatrix<N>::type fRinv;  ///< Residual inverse error matrix.
    mutable typename KHMatrix<N>::type fH;       ///< Kalman H-matrix.
    mutable double fChisq;                       ///< Incremental chisquare.
  };

  // Method implementations.

  /// Default constructor.
  template <int N>
  KHit<N>::KHit() :
    fChisq(0.)
  {}

  /// Initializing Constructor -- surface only.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  ///
  template<int N>
  KHit<N>::KHit(const boost::shared_ptr<const Surface>& psurf) :
    KHitBase(psurf),
    fChisq(0.)
  {}

  /// Fully Initializing Constructor.
  ///
  /// Arguments:
  ///
  /// psurf - Surface pointer.
  /// mvec  - Measurement vector.
  /// merr  - Measurement error.
  ///
  template<int N>
  KHit<N>::KHit(const boost::shared_ptr<const Surface>& psurf,
		const typename KVector<N>::type& mvec,
		const typename KSymMatrix<N>::type& merr) :
    KHitBase(psurf),
    fMvec(mvec),
    fMerr(merr),
    fChisq(0.)
  {}

  /// Destructor.
  template <int N>
  KHit<N>::~KHit()
  {}

  /// Prediction method.
  ///
  /// Arguments;
  ///
  /// tre - Track prediction.
  ///
  template <int N>
  bool KHit<N>::predict(const KETrack& tre) const
  {
    // Update prediction vector, error matrox, and H-matrix.
    // Method subpredict should throw an exception if track
    // surface does not match measurement surface.

    bool ok = subpredict(tre, fPvec, fPerr, fH);
    if(!ok)
      return ok;

    // Update residual

    fRvec = fMvec - fPvec;
    fRerr = fMerr + fPerr;
    fRinv = fRerr;
    ok = syminvert(fRinv);
    if(!ok)
      return ok;

    // Calculate incremental chisquare.

    typename KVector<N>::type rtemp = prod(fRinv, fRvec);
    fChisq = inner_prod(fRvec, rtemp);

    return true;
  }

  /// Update track method.
  ///
  /// Arguments:
  ///
  /// tre - Track to be updated.
  ///
  template <int N>
  void KHit<N>::update(KETrack& tre) const
  {
    // Make sure that the track surface and the measurement surface are the same.
    // Throw an exception if they are not.

    if(!getSurface()->isEqual(*tre.getSurface()))
      throw cet::exception("KHit") << "Track surface not the same as measurement surface.\n";

    const TrackVector& tvec = tre.getVector();
    const TrackError& terr = tre.getError();
    TrackVector::size_type size = tvec.size();

    // Calculate gain matrix.

    typename KGMatrix<N>::type temp(size, N);
    typename KGMatrix<N>::type gain(size, N);
    temp = prod(trans(fH), fRinv);
    gain = prod(terr, temp);

    // Calculate updated track state.

    TrackVector newvec = tre.getVector() + prod(gain, fRvec);

    // Calculate updated error matrix.

    TrackMatrix fact = ublas::identity_matrix<TrackVector::value_type>(size);
    fact -= prod(gain, fH);
    TrackError newerr = prod(fact, terr);

    // Update track.

    tre.setVector(newvec);
    tre.setError(newerr);
  }
}

#endif
