//
// File: LATest.cxx
//
// Purpose: Test linear algebra.
//

#include <iostream>
#include <cassert>
#include <cmath>
#include "TrackFinder/KalmanLinearAlgebra.h"
#include "boost/numeric/ublas/io.hpp"

int main()
{
  // Make sure assert is enabled.

  bool assert_flag = false;
  assert((assert_flag = true, assert_flag));
  if ( ! assert_flag ) {
    std::cerr << "Assert is disabled" << std::endl;
    return 1;
  }

  // Test matrix inversion.

  // 1x1 KSymMatrix.

  trkf::KSymMatrix<1>::type m1;
  for(unsigned int i = 0; i < m1.size1(); ++i) {
    for(unsigned int j = 0; j <= i; ++j)
      m1(i,j) = i+j+1.;
  }
  trkf::KSymMatrix<1>::type minv1(m1);
  bool ok = trkf::syminvert(minv1);
  assert(ok);
  trkf::ublas::matrix<double> unit1 = prod(m1, minv1);
  std::cout << m1 << std::endl;
  std::cout << minv1 << std::endl;
  std::cout << unit1 << std::endl;
  for(unsigned int i = 0; i < m1.size1(); ++i) {
    for(unsigned int j = 0; j <= i; ++j) {
      double val = (i==j ? 1. : 0.);
      assert(std::abs(unit1(i,j) - val) < 1.e-10);
    }
  }

  // 2x2 KSymMatrix.

  trkf::KSymMatrix<2>::type m2;
  for(unsigned int i = 0; i < m2.size1(); ++i) {
    for(unsigned int j = 0; j <= i; ++j)
      m2(i,j) = i+j+1.;
  }
  trkf::KSymMatrix<2>::type minv2(m2);
  ok = trkf::syminvert(minv2);
  assert(ok);
  trkf::ublas::matrix<double> unit2 = prod(m2, minv2);
  std::cout << m2 << std::endl;
  std::cout << minv2 << std::endl;
  std::cout << unit2 << std::endl;
  for(unsigned int i = 0; i < m2.size1(); ++i) {
    for(unsigned int j = 0; j <= i; ++j) {
      double val = (i==j ? 1. : 0.);
      assert(std::abs(unit2(i,j) - val) < 1.e-10);
    }
  }

  // 3x3 KSymMatrix.

  trkf::KSymMatrix<3>::type m3;
  for(unsigned int i = 0; i < m3.size1(); ++i) {
    for(unsigned int j = 0; j <= i; ++j) {
      m3(i,j) = i+j+1.;
      if(i==j)
	m3(i,j) += 1.;
    }
  }
  trkf::KSymMatrix<3>::type minv3(m3);
  ok = trkf::syminvert(minv3);
  assert(ok);
  trkf::ublas::matrix<double> unit3 = prod(m3, minv3);
  std::cout << m3 << std::endl;
  std::cout << minv3 << std::endl;
  std::cout << unit3 << std::endl;
  for(unsigned int i = 0; i < m3.size1(); ++i) {
    for(unsigned int j = 0; j <= i; ++j) {
      double val = (i==j ? 1. : 0.);
      assert(std::abs(unit3(i,j) - val) < 1.e-10);
    }
  }

  // 4x4 KSymMatrix.

  trkf::KSymMatrix<4>::type m4;
  for(unsigned int i = 0; i < m4.size1(); ++i) {
    for(unsigned int j = 0; j <= i; ++j) {
      m4(i,j) = i+j+1.;
      if(i==j)
	m4(i,j) += 1.;
    }
  }
  trkf::KSymMatrix<4>::type minv4(m4);
  ok = trkf::syminvert(minv4);
  assert(ok);
  trkf::ublas::matrix<double> unit4 = prod(m4, minv4);
  std::cout << m4 << std::endl;
  std::cout << minv4 << std::endl;
  std::cout << unit4 << std::endl;
  for(unsigned int i = 0; i < m4.size1(); ++i) {
    for(unsigned int j = 0; j <= i; ++j) {
      double val = (i==j ? 1. : 0.);
      assert(std::abs(unit4(i,j) - val) < 1.e-10);
    }
  }

  // 5x5 TrackError.

  trkf::TrackError m5;
  for(unsigned int i = 0; i < m5.size1(); ++i) {
    for(unsigned int j = 0; j <= i; ++j) {
      m5(i,j) = i+j+1.;
      if(i==j)
	m5(i,j) += 1.;
    }
  }
  trkf::TrackError minv5(m5);
  ok = trkf::syminvert(minv5);
  assert(ok);
  trkf::ublas::matrix<double> unit5 = prod(m5, minv5);
  std::cout << m5 << std::endl;
  std::cout << minv5 << std::endl;
  std::cout << unit5 << std::endl;
  for(unsigned int i = 0; i < m5.size1(); ++i) {
    for(unsigned int j = 0; j <= i; ++j) {
      double val = (i==j ? 1. : 0.);
      assert(std::abs(unit5(i,j) - val) < 1.e-10);
    }
  }

  // Done (success).

  std::cout << "LATest: All tests passed." << std::endl;

  return 0;
}
