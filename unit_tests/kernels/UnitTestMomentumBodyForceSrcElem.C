/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernels/UnitTestKernelUtils.h"
#include "UnitTestUtils.h"

#include "MomentumBodyForceSrcElemKernel.h"

#include <random>

TEST_F(MomentumKernelHex8Mesh, body_force)
{
  std::mt19937 rng;
  rng.seed(0); // fixed seed
  std::uniform_real_distribution<double> coeff(-1, +1);

  fill_mesh_and_init_fields();

  // Setup solution options for default advection kernel
  solnOpts_.meshMotion_ = false;
  solnOpts_.meshDeformation_ = false;
  solnOpts_.externalMeshDeformation_ = false;

  std::vector<double> bfVec = { coeff(rng), coeff(rng), coeff(rng) };
  solnOpts_.srcTermParamMap_.insert({"body_force", bfVec});

  // Initialize the kernel driver
  unit_test_kernel_utils::TestKernelDriver assembleKernels(
    bulk_, partVec_, coordinates_, spatialDim_, stk::topology::HEX_8);

  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> kernel(
    new sierra::nalu::MomentumBodyForceSrcElemKernel<sierra::nalu::AlgTraitsHex8>(
      bulk_, solnOpts_, assembleKernels.dataNeededByKernels_));

  // Add to kernels to be tested
  assembleKernels.activeKernels_.push_back(kernel.get());

  EXPECT_NO_THROW(assembleKernels.execute());

  EXPECT_EQ(assembleKernels.lhs_.dimension(0), 24u);
  EXPECT_EQ(assembleKernels.lhs_.dimension(1), 24u);
  EXPECT_EQ(assembleKernels.rhs_.dimension(0), 24u);

  // Exact solution
  std::vector<double> rhsExact(24,0.0);
  for (size_t i= 0; i < 24; i += 3) {
    rhsExact[i + 0] = 0.125 * bfVec[0];
    rhsExact[i + 1] = 0.125 * bfVec[1];
    rhsExact[i + 2] = 0.125 * bfVec[2];
  }

  unit_test_kernel_utils::expect_all_near(assembleKernels.rhs_,rhsExact.data());
}
