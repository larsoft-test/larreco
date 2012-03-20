////////////////////////////////////////////////////////////////////////
///
/// \file   KalmanLinearAlgebra.h
///
/// \brief  Kalman filter linear algebra classes.
///
/// \author H. Greenlee 
///
/// There are six Kalman filter linear algebra classes defined
/// in this header:
///
/// 1. TrackVector - Track state vector, dimension 5
/// 2. TrackError - Track error matrix, dimension 5x5.
/// 3. MeasVector - Measurement vector, dimension N.
/// 4. MeasError - Measurement error matrix, dimension NxN.
/// 5. HMatrix - Kalman H matrix, dimension Nx5.
/// 6. KMatrix - Kalman gain matrix, dimension 5XN.
///
/// The above classes derive from some ublas (boost/numeric/ublas)
/// linear algebra class, nonvirtually inheriting all methods, 
/// except for locally defined default and initial-value constructors,
/// plus compiler-generated methods.
///
/// We do this because of the unwieldy syntax of the underlying ublas 
/// linear algebra package.  It would have been reasonable to use
/// typedefs instead of derived classes, except that the last four
/// classes are templates, and c++ does not allow templated typedefs.
///
/// All linear algebra objects use the following storage model.
///
/// 1.  Matrices are stored in row major order (normal c/c++ convention).
/// 2.  Symmetric matrices are stored in lower triangular format.
/// 3.  Size of objects is specified at compilation time and
///     memory is preallocated on the stack (using ublass bounded_array
///     container).
/// 4.  Virtual functions and allocators are avoided entirely.
///
////////////////////////////////////////////////////////////////////////

#ifndef KALMANLINEARALGEBRA_H
#define KALMANLINEARALGEBRA_H

#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/symmetric.hpp"

namespace trkf {

  /// Define a shortened alias for ublas namespace.
  namespace ublas = boost::numeric::ublas;

  /// Track state vector, dimension 5.
  class TrackVector : public ublas::vector<double, ublas::bounded_array<double, 5> >
  {
  public:

    /// Default constructor.
    TrackVector() : ublas::vector<double, ublas::bounded_array<double, 5> >(5) {}

    /// Initial value constructor.
    TrackVector(double val) : ublas::vector<double, ublas::bounded_array<double, 5> >(5, val) {}
  };

  /// Track error matrix, dimension 5x5.
  class TrackError : public ublas::symmetric_matrix<double, ublas::lower, ublas::row_major, ublas::bounded_array<double, 25> >
  {
  public:

    /// Default constructor.
    TrackError() : ublas::symmetric_matrix<double, ublas::lower, ublas::row_major, ublas::bounded_array<double, 25> >(5) {}
  };

  /// Measurement vector, dimension N.
  template<int N>
  class MeasVector : public ublas::vector<double, ublas::bounded_array<double, N> >
  {
  public:

    /// Default constructor.
    MeasVector() : ublas::vector<double, ublas::bounded_array<double, N> >(N) {}

    /// Initial value constructor.
    MeasVector(double val) : ublas::vector<double, ublas::bounded_array<double, N> >(N, val) {}
  };

  /// Measurement error matrix, dimension NxN.
  template<int N>
  class MeasError : public ublas::symmetric_matrix<double, ublas::lower, ublas::row_major, ublas::bounded_array<double, N*N> >
  {
  public:

    /// Default constructor.
    MeasError() : ublas::symmetric_matrix<double, ublas::lower, ublas::row_major, ublas::bounded_array<double, N*N> >(N) {}
  };

  /// Kalman H-matrix matrix, dimension Nx5.
  template<int N>
  class HMatrix : public ublas::matrix<double, ublas::row_major, ublas::bounded_array<double, 5*N> >
  {
  public:

    /// Default constructor.
    HMatrix() : ublas::matrix<double, ublas::row_major, ublas::bounded_array<double, 5*N> >(N, 5) {}

    /// Initial value constructor.
    HMatrix(double val) : ublas::matrix<double, ublas::row_major, ublas::bounded_array<double, 5*N> >(N, 5, val) {}
  };

  /// Kalman gain matrix matrix, dimension 5xN.
  template<int N>
  class KMatrix : public ublas::matrix<double, ublas::row_major, ublas::bounded_array<double, 5*N> >
  {
  public:

    /// Default constructor.
    KMatrix() : ublas::matrix<double, ublas::row_major, ublas::bounded_array<double, 5*N> >(5, N) {}

    /// Initial value constructor.
    KMatrix(double val) : ublas::matrix<double, ublas::row_major, ublas::bounded_array<double, 5*N> >(5, N, val) {}
  };
}

#endif
