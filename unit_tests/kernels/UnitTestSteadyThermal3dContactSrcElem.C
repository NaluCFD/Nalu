/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernels/UnitTestKernelUtils.h"
#include "UnitTestUtils.h"

#include "user_functions/SteadyThermal3dContactSrcElemKernel.h"

namespace {
namespace hex8_golds {
namespace steady_3d_thermal {

static constexpr double rhs[8] = {
  2.26627114475e-16, -7.55423714915e-17,
  -3.77711857458e-16, -7.55423714915e-17,
  -7.55423714915e-17, -3.77711857458e-16,
  -6.79881343424e-16, -3.77711857458e-16,
};
} // steady_3d_thermal
} // hex8_golds
} // anonymous namespace

/// Steady 3D MMS source term
TEST_F(HeatCondKernelHex8Mesh, steady_3d_thermal)
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
    new sierra::nalu::SteadyThermal3dContactSrcElemKernel<sierra::nalu::AlgTraitsHex8>(
      bulk_, solnOpts_, assembleKernels.dataNeededByKernels_));

  // Add to kernels to be tested
  assembleKernels.activeKernels_.push_back(kernel.get());

  assembleKernels.execute();

  EXPECT_EQ(assembleKernels.lhs_.dimension(0), 8u);
  EXPECT_EQ(assembleKernels.lhs_.dimension(1), 8u);
  EXPECT_EQ(assembleKernels.rhs_.dimension(0), 8u);

  namespace gold_values = ::hex8_golds::steady_3d_thermal;

  unit_test_kernel_utils::expect_all_near(assembleKernels.rhs_, gold_values::rhs);
  unit_test_kernel_utils::expect_all_near<8>(assembleKernels.lhs_, 0.0);
}
