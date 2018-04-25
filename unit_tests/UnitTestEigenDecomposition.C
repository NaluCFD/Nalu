#include <gtest/gtest.h>
#include <limits>
#include <random>
#include <stdexcept>

#include "EigenDecomposition.h"

// NGP-based includes
#include "SimdInterface.h"
#include "KokkosInterface.h"

#include "UnitTestUtils.h"


#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <master_element/MasterElement.h>

#include <fstream>



namespace {

    std::random_device rd; 
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> rand(0.0, 1.0);
 
    double a11 = rand(rng);
    double a22 = rand(rng);
    double a33 = rand(rng);
    double a12 = rand(rng);
    double a13 = rand(rng);
    double a23 = rand(rng);

static double A3d_rand[3][3] = {
  {a11, a12, a13}, {a12, a22, a23}, {a13, a23, a33},
};

static double A3d_fixed[3][3] = {
  {0.0562500000000000,  0.1875000000000000,  0.0125000000000000},
  {0.1875000000000000,  0.0093750000000000, -0.4218750000000000},
  {0.0125000000000000, -0.4218750000000000,  0.0031250000000000},
};

static double A2d_rand[2][2] = {
  {a11, a12}, {a12, a22},
};

static double A2d_fixed[2][2] = {
  { 0.0562500000000000, -0.0187500000000000},  
  {-0.018750000000000,   0.0093750000000000}, 
};

} // namespace

// This tests whether the correct eigenvalues are obtained 
TEST(TestEigen, testeigendecomp3d)
{
  double Q_[3][3], D_[3][3];

  sierra::nalu::EigenDecomposition eigSolver;
  eigSolver.sym_diagonalize(A3d_fixed, Q_, D_);

  // Perform tests -- lambda evaluated in Mathematica
  const double tol = 1e-15;
  const double lambda_gold[3] = {
      0.056736539229635605, 0.46782517604126655, -0.45581171527090225};

  for (unsigned j = 0 ; j < 3; ++j) {
    EXPECT_NEAR(D_[j][j], lambda_gold[j], tol);
  }
}

TEST(TestEigen, testeigendecomp2d)
{
  double Q_[2][2], D_[2][2];

  sierra::nalu::EigenDecomposition eigSolver;
  eigSolver.sym_diagonalize(A2d_fixed, Q_, D_);

  // Perform tests -- lambda evaluated in Mathematica
  const double tol = 1e-15;
  const double lambda_gold[2] = {
      0.06282714486296648, 0.0027978551370335227};

  for (unsigned j = 0; j < 2; ++j) {
    EXPECT_NEAR(D_[j][j], lambda_gold[j], tol);
  }
}

// This tests that the eignevalue decomposition and reconstruction
// returns the original matrix
TEST(TestEigen, testeigendecompandreconstruct3d)
{
  double b_[3][3], Q_[3][3], D_[3][3];

  sierra::nalu::EigenDecomposition eigSolver;
  eigSolver.sym_diagonalize(A3d_rand, Q_, D_);
  eigSolver.reconstruct_matrix_from_decomposition(D_, Q_, b_);

  // Perform tests
  const double tol = 1e-15;

  // Reconstructed matrix should be within tol of original
  for (unsigned j = 0 ; j < 3; ++j) {
    for (unsigned i = 0; i < 3; ++i) {
      EXPECT_NEAR(b_[i][j], A3d_rand[i][j], tol);
    }
  }
}

// This tests that the eignevalue decomposition and reconstruction
// returns the original matrix
TEST(TestEigen, testeigendecompandreconstruct2d)
{
  double b_[2][2], Q_[2][2], D_[2][2];

  sierra::nalu::EigenDecomposition eigSolver;
  eigSolver.sym_diagonalize(A2d_rand, Q_, D_);
  eigSolver.reconstruct_matrix_from_decomposition(D_, Q_, b_);

  // Perform tests
  const double tol = 1e-15;
  
  // Reconstructed matrix should be within tol of original
  for (unsigned j = 0 ; j < 2; ++j) {
    for (unsigned i = 0; i < 2; ++i) {
      EXPECT_NEAR(b_[i][j], A2d_rand[i][j], tol);
    }
  }
}

