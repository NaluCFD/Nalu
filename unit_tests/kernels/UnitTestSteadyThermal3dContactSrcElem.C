/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernels/UnitTestKernelUtils.h"
#include "kernels/UnitTestKernelGolds.h"
#include "UnitTestUtils.h"

#include "user_functions/SteadyThermal3dContactSrcElemKernel.h"

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

  namespace gold_values = unit_test_golds::hex8_golds::heatcond::steady_3d_thermal;

  unit_test_kernel_utils::expect_all_near(assembleKernels.rhs_, gold_values::rhs);
  unit_test_kernel_utils::expect_all_near<8>(assembleKernels.lhs_, 0.0);
}
