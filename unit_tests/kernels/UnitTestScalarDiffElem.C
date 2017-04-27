/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernels/UnitTestKernelUtils.h"
#include "UnitTestUtils.h"

#include "ScalarDiffElemKernel.h"

namespace {
namespace hex8_golds {
namespace heatcond {
namespace cvfem_diff {
static constexpr double lhs[8][8] = {
  {0.421875, -0.046875, -0.078125, -0.046875, -0.046875, -0.078125, -0.046875, -0.078125},
  {-0.046875, 0.421875, -0.046875, -0.078125, -0.078125, -0.046875, -0.078125, -0.046875},
  {-0.078125, -0.046875, 0.421875, -0.046875, -0.046875, -0.078125, -0.046875, -0.078125},
  {-0.046875, -0.078125, -0.046875, 0.421875, -0.078125, -0.046875, -0.078125, -0.046875},
  {-0.046875, -0.078125, -0.046875, -0.078125, 0.421875, -0.046875, -0.078125, -0.046875},
  {-0.078125, -0.046875, -0.078125, -0.046875, -0.046875, 0.421875, -0.046875, -0.078125},
  {-0.046875, -0.078125, -0.046875, -0.078125, -0.078125, -0.046875, 0.421875, -0.046875},
  {-0.078125, -0.046875, -0.078125, -0.046875, -0.046875, -0.078125, -0.046875, 0.421875},
};
} // cvfem_diff
} // heatcond
} // hex8_golds
} // anonymous namespace


/// Scalar diffusion kernel applied to heat conduction equation
TEST_F(HeatCondKernelHex8Mesh, cvfem_diff)
{
  fill_mesh_and_init_fields();

  // Setup solution options for default advection kernel
  solnOpts_.meshMotion_ = false;
  solnOpts_.meshDeformation_ = false;
  solnOpts_.externalMeshDeformation_ = false;

  // Initialize the kernel driver
  unit_test_kernel_utils::TestKernelDriver assembleKernels(
    bulk_, partVec_, coordinates_, 1, stk::topology::HEX_8);

  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> kernel(
    new sierra::nalu::ScalarDiffElemKernel<sierra::nalu::AlgTraitsHex8>(
      bulk_, solnOpts_, temperature_, thermalCond_,
      assembleKernels.dataNeededByKernels_));

  // Add to kernels to be tested
  assembleKernels.activeKernels_.push_back(kernel.get());

  assembleKernels.execute();

  EXPECT_EQ(assembleKernels.lhs_.dimension(0), 8u);
  EXPECT_EQ(assembleKernels.lhs_.dimension(1), 8u);
  EXPECT_EQ(assembleKernels.rhs_.dimension(0), 8u);

  namespace gold_values = hex8_golds::heatcond::cvfem_diff;
  unit_test_kernel_utils::expect_all_near(assembleKernels.rhs_, 0.0);
  unit_test_kernel_utils::expect_all_near<8>(assembleKernels.lhs_, gold_values::lhs);
}
