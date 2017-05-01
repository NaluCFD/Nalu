/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernels/UnitTestKernelUtils.h"
#include "UnitTestUtils.h"

#include "ScalarAdvDiffElemKernel.h"

namespace {
namespace hex8_golds {
namespace advection_diffusion {

static constexpr double rhs[8] = {
  -0.0000138750000000, 0.0000138750000000, -0.0000138750000000, 0.0000138750000000, 
  0.0000138750000000, -0.0000138750000000, 0.0000138750000000, -0.0000138750000000, };

static constexpr double lhs[8][8] = {
  { 0.0000078046875000, 0.0001495247342374, -0.0000014453125000, -0.0001512591092374, 
    -0.0000008671875000, 0.0000486853280791, -0.0000008671875000, -0.0000515759530791,  },
  { -0.0004520429527123, -0.0008945468429247, -0.0007528267961872, -0.0003022291559749, 
    -0.0001518372342374, -0.0003016510309749, -0.0002520985153957, -0.0001011284686583,  },
  { -0.0000014453125000, 0.0004503085777123, 0.0000078046875000, -0.0004520429527123, 
    -0.0000008671875000, 0.0001489466092374, -0.0000008671875000, -0.0001518372342374,  },
  { 0.0004503085777123, 0.0002993385309749, 0.0007510924211872, 0.0009101562179247, 
    0.0001489466092374, 0.0000993940936583, 0.0002492078903957, 0.0002999166559749,  },
  { -0.0000008671875000, -0.0000133055695312, -0.0000008671875000, 0.0000104149445312, 
    0.0000078046875000, -0.0000364479585937, -0.0000014453125000, 0.0000347135835937,  },
  { 0.0000341354585937, 0.0000702943546875, 0.0000578559726562, 0.0000228533265625, 
    0.0001058751257812, 0.0002212893140624, 0.0001770366679686, 0.0000697162296875,  },
  { -0.0000008671875000, -0.0000370260835937, -0.0000008671875000, 0.0000341354585937, 
    -0.0000014453125000, -0.0001076095007812, 0.0000078046875000, 0.0001058751257812,  },
  { -0.0000370260835937, -0.0000245877015625, -0.0000607465976562, -0.0000720287296875, 
    -0.0001076095007812, -0.0000726068546875, -0.0001787710429686, -0.0002056799390624,  },
};

} // advection_diffusion
} // hex8_golds
} // anonymous namespace

/// Scalar advection/diffusion (will use mixture fraction as scalar)
TEST_F(MixtureFractionKernelHex8Mesh, advection_diffusion)
{
  // FIXME: only test on one core
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) 
    return;

  fill_mesh_and_init_fields();

  // Setup solution options for default advection kernel
  solnOpts_.meshMotion_ = false;
  solnOpts_.meshDeformation_ = false;
  solnOpts_.externalMeshDeformation_ = false;

  // Initialize the kernel driver
  unit_test_kernel_utils::TestKernelDriver assembleKernels(
    bulk_, partVec_, coordinates_, 1, stk::topology::HEX_8);

  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> advKernel(
    new sierra::nalu::ScalarAdvDiffElemKernel<sierra::nalu::AlgTraitsHex8>(
     bulk_, solnOpts_, mixFraction_, viscosity_, assembleKernels.dataNeededByKernels_));

  // Register the kernel for execution
  assembleKernels.activeKernels_.push_back(advKernel.get());

  // Populate LHS and RHS
  assembleKernels.execute();

  EXPECT_EQ(assembleKernels.lhs_.dimension(0), 8u);
  EXPECT_EQ(assembleKernels.lhs_.dimension(1), 8u);
  EXPECT_EQ(assembleKernels.rhs_.dimension(0), 8u);

  namespace gold_values = hex8_golds::advection_diffusion;
  unit_test_kernel_utils::expect_all_near(assembleKernels.rhs_, gold_values::rhs);
  unit_test_kernel_utils::expect_all_near<8>(assembleKernels.lhs_, gold_values::lhs);
}
