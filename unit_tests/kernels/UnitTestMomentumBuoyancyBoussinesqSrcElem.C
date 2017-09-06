/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernels/UnitTestKernelUtils.h"
#include "UnitTestUtils.h"

#include "MomentumBuoyancyBoussinesqSrcElemKernel.h"

#include <random>

TEST_F(MomentumKernelHex8Mesh, buoyancy_boussinesq)
{
  std::mt19937 rng;
  rng.seed(0); // fixed seed
  std::uniform_real_distribution<double> ref_densities(0.8,1.3);

  fill_mesh_and_init_fields();

  // Setup solution options for default advection kernel
  solnOpts_.meshMotion_ = false;
  solnOpts_.meshDeformation_ = false;
  solnOpts_.externalMeshDeformation_ = false;
  solnOpts_.gravity_.resize(spatialDim_, 0.0);
  solnOpts_.gravity_[2] = -9.81;
  solnOpts_.referenceDensity_ = ref_densities(rng);
  solnOpts_.referenceTemperature_ = 298;
  solnOpts_.thermalExpansionCoeff_ = 1.0;

  // Initialize the kernel driver
  unit_test_kernel_utils::TestKernelDriver assembleKernels(
    bulk_, partVec_, coordinates_, spatialDim_, stk::topology::HEX_8);

  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> kernel(
    new sierra::nalu::MomentumBuoyancyBoussinesqSrcElemKernel<sierra::nalu::AlgTraitsHex8>(
      bulk_, solnOpts_, assembleKernels.dataNeededByKernels_));

  // Add to kernels to be tested
  assembleKernels.activeKernels_.push_back(kernel.get());

  assembleKernels.execute();

  EXPECT_EQ(assembleKernels.lhs_.dimension(0), 24u);
  EXPECT_EQ(assembleKernels.lhs_.dimension(1), 24u);
  EXPECT_EQ(assembleKernels.rhs_.dimension(0), 24u);

  double expFac = solnOpts_.referenceDensity_ * solnOpts_.thermalExpansionCoeff_;

  // Exact solution
  std::vector<double> rhsExact(24,0.0);
  for (size_t i=2; i < 24; i += 3)
    rhsExact[i] = 0.125 * solnOpts_.gravity_[2] * expFac * (300.0 - solnOpts_.referenceTemperature_);

  unit_test_kernel_utils::expect_all_near(assembleKernels.rhs_,rhsExact.data());
}
