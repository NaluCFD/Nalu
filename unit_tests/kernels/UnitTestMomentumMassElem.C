/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernels/UnitTestKernelUtils.h"
#include "UnitTestUtils.h"

#include "MomentumMassElemKernel.h"

namespace {
namespace hex8_golds {
namespace momentum_time_derivative {
static constexpr double rhs[24] = {
  4.04377248484e-15, 5.57458098377e-15, 0,
  1.36621259534e-14, 7.10538948271e-15, 0,
  1.21313174545e-14, 1.67237429513e-14, 0,
  2.5129639859e-15,  1.51929344524e-14, 0,
  4.04377248484e-15, 5.57458098377e-15, 0,
  1.36621259534e-14, 7.10538948271e-15, 0,
  1.21313174545e-14, 1.67237429513e-14, 0,
  2.5129639859e-15,  1.51929344524e-14, 0,
};

static constexpr double lhs[24][24] = {
  {
    0.52734375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0,
    0.17578125, 0, 0, 0.05859375, 0, 0, 0.01953125, 0, 0, 0.05859375, 0, 0,
  },
  {
    0, 0.52734375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0,
    0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.01953125, 0, 0, 0.05859375, 0,
  },
  {
    0, 0, 0.52734375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.17578125,
    0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.01953125, 0, 0, 0.05859375,
  },
  {
    0.17578125, 0, 0, 0.52734375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0,
    0.05859375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.01953125, 0, 0,
  },
  {
    0, 0.17578125, 0, 0, 0.52734375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0,
    0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.01953125, 0,
  },
  {
    0, 0, 0.17578125, 0, 0, 0.52734375, 0, 0, 0.17578125, 0, 0, 0.05859375,
    0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.01953125,
  },
  {
    0.05859375, 0, 0, 0.17578125, 0, 0, 0.52734375, 0, 0, 0.17578125, 0, 0,
    0.01953125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0,
  },
  {
    0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.52734375, 0, 0, 0.17578125, 0,
    0, 0.01953125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0,
  },
  {
    0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.52734375, 0, 0, 0.17578125,
    0, 0, 0.01953125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.05859375,
  },
  {
    0.17578125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.52734375, 0, 0,
    0.05859375, 0, 0, 0.01953125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0,
  },
  {
    0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.52734375, 0,
    0, 0.05859375, 0, 0, 0.01953125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0,
  },
  {
    0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.52734375,
    0, 0, 0.05859375, 0, 0, 0.01953125, 0, 0, 0.05859375, 0, 0, 0.17578125,
  },
  {
    0.17578125, 0, 0, 0.05859375, 0, 0, 0.01953125, 0, 0, 0.05859375, 0, 0,
    0.52734375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0,
  },
  {
    0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.01953125, 0, 0, 0.05859375, 0,
    0, 0.52734375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0,
  },
  {
    0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.01953125, 0, 0, 0.05859375,
    0, 0, 0.52734375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.17578125,
  },
  {
    0.05859375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.01953125, 0, 0,
    0.17578125, 0, 0, 0.52734375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0,
  },
  {
    0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.01953125, 0,
    0, 0.17578125, 0, 0, 0.52734375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0,
  },
  {
    0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.01953125,
    0, 0, 0.17578125, 0, 0, 0.52734375, 0, 0, 0.17578125, 0, 0, 0.05859375,
  },
  {
    0.01953125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0,
    0.05859375, 0, 0, 0.17578125, 0, 0, 0.52734375, 0, 0, 0.17578125, 0, 0,
  },
  {
    0, 0.01953125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.05859375, 0,
    0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.52734375, 0, 0, 0.17578125, 0,
  },
  {
    0, 0, 0.01953125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.05859375,
    0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.52734375, 0, 0, 0.17578125,
  },
  {
    0.05859375, 0, 0, 0.01953125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0,
    0.17578125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.52734375, 0, 0,
  },
  {
    0, 0.05859375, 0, 0, 0.01953125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0,
    0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.52734375, 0,
  },
  {
    0, 0, 0.05859375, 0, 0, 0.01953125, 0, 0, 0.05859375, 0, 0, 0.17578125,
    0, 0, 0.17578125, 0, 0, 0.05859375, 0, 0, 0.17578125, 0, 0, 0.52734375,
  },
};

} // momentum_time_derivative

namespace momentum_time_derivative_lumped {
static constexpr double rhs[24] = {
  0, 0, 0,
  1.92367069372e-14, 3.06161699787e-15, 0,
  1.61750899393e-14, 2.22983239351e-14, 0,
  -3.06161699787e-15, 1.92367069372e-14, 0,
  0, 0, 0,
  1.92367069372e-14, 3.06161699787e-15, 0,
  1.61750899393e-14, 2.22983239351e-14, 0,
  -3.06161699787e-15, 1.92367069372e-14, 0,
};
} // momentum_time_derivative_lumped
} // hex8_golds
} // anonymous namespace

TEST_F(MomentumKernelHex8Mesh, momentum_time_derivative)
{
  fill_mesh_and_init_fields();

  // Setup solution options for default advection kernel
  solnOpts_.meshMotion_ = false;
  solnOpts_.meshDeformation_ = false;
  solnOpts_.externalMeshDeformation_ = false;

  // Initialize the kernel driver
  unit_test_kernel_utils::TestKernelDriver assembleKernels(
    bulk_, partVec_, coordinates_, spatialDim_, stk::topology::HEX_8);

  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> massKernel(
    new sierra::nalu::MomentumMassElemKernel<sierra::nalu::AlgTraitsHex8>(
      bulk_, solnOpts_, assembleKernels.dataNeededByKernels_, false));

  // Add to kernels to be tested
  assembleKernels.activeKernels_.push_back(massKernel.get());

  // Mass terms need time integration information
  sierra::nalu::TimeIntegrator timeIntegrator;
  timeIntegrator.timeStepN_ = 0.1;
  timeIntegrator.timeStepNm1_ = 0.1;
  timeIntegrator.gamma1_ = 1.0;
  timeIntegrator.gamma2_ = -1.0;
  timeIntegrator.gamma3_ = 0.0;

  // Call Kernel setup with the time integrator to setup Kernel values
  massKernel->setup(timeIntegrator);

  // Populate LHS and RHS
  assembleKernels.execute();

  EXPECT_EQ(assembleKernels.lhs_.dimension(0), 24u);
  EXPECT_EQ(assembleKernels.lhs_.dimension(1), 24u);
  EXPECT_EQ(assembleKernels.rhs_.dimension(0), 24u);

  namespace gold_values = ::hex8_golds::momentum_time_derivative;
  unit_test_kernel_utils::expect_all_near(assembleKernels.rhs_, gold_values::rhs);
  unit_test_kernel_utils::expect_all_near<24>(assembleKernels.lhs_, gold_values::lhs);
}

TEST_F(MomentumKernelHex8Mesh, momentum_time_derivative_lumped)
{
  fill_mesh_and_init_fields();

  // Setup solution options for default advection kernel
  solnOpts_.meshMotion_ = false;
  solnOpts_.meshDeformation_ = false;
  solnOpts_.externalMeshDeformation_ = false;

  // Initialize the kernel driver
  unit_test_kernel_utils::TestKernelDriver assembleKernels(
    bulk_, partVec_, coordinates_, spatialDim_, stk::topology::HEX_8);

  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> massKernel(
    new sierra::nalu::MomentumMassElemKernel<sierra::nalu::AlgTraitsHex8>(
      bulk_, solnOpts_, assembleKernels.dataNeededByKernels_, true));

  // Add to kernels to be tested
  assembleKernels.activeKernels_.push_back(massKernel.get());

  // Mass terms need time integration information
  sierra::nalu::TimeIntegrator timeIntegrator;
  timeIntegrator.timeStepN_ = 0.1;
  timeIntegrator.timeStepNm1_ = 0.1;
  timeIntegrator.gamma1_ = 1.0;
  timeIntegrator.gamma2_ = -1.0;
  timeIntegrator.gamma3_ = 0.0;

  // Call Kernel setup with the time integrator to setup Kernel values
  massKernel->setup(timeIntegrator);

  // Populate LHS and RHS
  assembleKernels.execute();

  namespace gold_values = ::hex8_golds::momentum_time_derivative_lumped;
  // Exact LHS expected
  std::vector<double> lhsExact(576, 0.0);
  for (int i=0; i<24; i++)
    lhsExact[i*24+i] = 1.25;

  EXPECT_EQ(assembleKernels.lhs_.dimension(0), 24u);
  EXPECT_EQ(assembleKernels.lhs_.dimension(1), 24u);
  EXPECT_EQ(assembleKernels.rhs_.dimension(0), 24u);
  unit_test_kernel_utils::expect_all_near(assembleKernels.rhs_, gold_values::rhs);
  unit_test_kernel_utils::expect_all_near(assembleKernels.lhs_, lhsExact.data());
}
