/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "PecletFunction.h"
#include "SimdInterface.h"

#include <vector>
#include <cmath>

namespace {

constexpr double tolerance = 1.0e-6;

}

TEST(PecletFunction, classic_double)
{
  const double A = 5.0;
  const double hybridFactor = 1.0;
  std::vector<double> pecletNumbers = {0.0, 1.0, std::sqrt(5.0), 1e5};
  std::vector<double> pecletFactors = {0.0, 1.0/6.0, 0.5, 1.0};

  sierra::nalu::ClassicPecletFunction<double> pecFunc(A, hybridFactor);

  for (int i=0; i < 4; i++) {
    EXPECT_NEAR(pecFunc.execute(pecletNumbers[i]),
                pecletFactors[i],
                tolerance);
  }
}

TEST(PecletFunction, classic_simd)
{
  const DoubleType A = 5.0;
  const DoubleType hybridFactor = 1.0;
  std::vector<DoubleType> pecletNumbers = {0.0, 1.0, std::sqrt(5.0), 1e5};
  std::vector<double> pecletFactors = {0.0, 1.0/6.0, 0.5, 1.0};

  sierra::nalu::ClassicPecletFunction<DoubleType> pecFunc(A, hybridFactor);

  for (int i=0; i < 4; i++) {
    const DoubleType pecFac = pecFunc.execute(pecletNumbers[i]);
    for (int is=0; is < stk::simd::ndoubles; is++) {
      EXPECT_NEAR(stk::simd::get_data(pecFac, is),
                  pecletFactors[i], tolerance);
    }
  }
}

TEST(PecletFunction, tanh_double)
{
  const double c1 = 5000.0;
  const double c2 = 200.0;
  std::vector<double> pecletNumbers = {-c1 - 10.0 * c2, c1, c1 + 10.0 * c2};
  std::vector<double> pecletFactors = {0.0, 0.5, 1.0};

  sierra::nalu::TanhFunction<double> pecFunc(c1, c2);

  for (int i=0; i < 3; i++) {
    EXPECT_NEAR(pecFunc.execute(pecletNumbers[i]),
                pecletFactors[i],
                tolerance);
  }
}

TEST(PecletFunction, tanh_simd)
{
  const DoubleType c1 = 5000.0;
  const DoubleType c2 = 200.0;
  std::vector<DoubleType> pecletNumbers = {-10.0 * c2, c1, c1 + 10.0 * c2};
  std::vector<double> pecletFactors = {0.0, 0.5, 1.0};

  sierra::nalu::TanhFunction<DoubleType> pecFunc(c1, c2);

  for (int i=0; i < 3; i++) {
    const DoubleType pecFac = pecFunc.execute(pecletNumbers[i]);
    for (int is=0; is < stk::simd::ndoubles; is++) {
      EXPECT_NEAR(stk::simd::get_data(pecFac, is),
                  pecletFactors[i], tolerance);
    }
  }
}
