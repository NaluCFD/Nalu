/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernels/UnitTestKernelUtils.h"
#include "kernels/UnitTestKernelGolds.h"
#include "UnitTestUtils.h"

#include "ContinuityAdvElemKernel.h"

/// Continuity advection with default Solution options
TEST_F(ContinuityKernelHex8Mesh, advection_default)
{
  fill_mesh_and_init_fields();

  // Setup solution options for default advection kernel
  solnOpts_.meshMotion_ = false;
  solnOpts_.meshDeformation_ = false;
  solnOpts_.externalMeshDeformation_ = false;
  solnOpts_.cvfemShiftMdot_ = false;
  solnOpts_.cvfemShiftPoisson_ = false;
  solnOpts_.cvfemReducedSensPoisson_ = false;
  solnOpts_.mdotInterpRhoUTogether_ = true;

  // Initialize the kernel driver
  unit_test_kernel_utils::TestKernelDriver assembleKernels(
    bulk_, partVec_, coordinates_, 1, stk::topology::HEX_8);

  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> advKernel(
    new sierra::nalu::ContinuityAdvElemKernel<sierra::nalu::AlgTraitsHex8>(
      bulk_, solnOpts_, assembleKernels.dataNeededByKernels_));

  // Register the kernel for execution
  assembleKernels.activeKernels_.push_back(advKernel.get());

  // Populate LHS and RHS
  assembleKernels.execute();

  namespace gold_values = unit_test_golds::hex8_golds::continuity::advection_default;

  unit_test_kernel_utils::expect_all_near(assembleKernels.rhs_, gold_values::rhs);
  unit_test_kernel_utils::expect_all_near<8>(assembleKernels.lhs_, gold_values::lhs);
}

/// Continuity advection kernel
///
/// `reduced_sens_cvfem_poisson: true`
///
TEST_F(ContinuityKernelHex8Mesh, advection_reduced_sens_cvfem_poisson)
{
  fill_mesh_and_init_fields();

  // Setup solution options for default advection kernel
  solnOpts_.meshMotion_ = false;
  solnOpts_.meshDeformation_ = false;
  solnOpts_.externalMeshDeformation_ = false;
  solnOpts_.cvfemShiftMdot_ = false;
  solnOpts_.cvfemShiftPoisson_ = false;
  solnOpts_.cvfemReducedSensPoisson_ = true;
  solnOpts_.mdotInterpRhoUTogether_ = true;

  // Initialize the kernel driver
  unit_test_kernel_utils::TestKernelDriver assembleKernels(
    bulk_, partVec_, coordinates_, 1, stk::topology::HEX_8);

  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> advKernel(
    new sierra::nalu::ContinuityAdvElemKernel<sierra::nalu::AlgTraitsHex8>(
      bulk_, solnOpts_, assembleKernels.dataNeededByKernels_));

  // Register the kernel for execution
  assembleKernels.activeKernels_.push_back(advKernel.get());

  // Populate LHS and RHS
  assembleKernels.execute();

  namespace gold_values =
    unit_test_golds::hex8_golds::continuity::advection_reduced_sensitivities;

  unit_test_kernel_utils::expect_all_near(assembleKernels.rhs_, gold_values::rhs);
  unit_test_kernel_utils::expect_all_near<8>(assembleKernels.lhs_, gold_values::lhs);
}
