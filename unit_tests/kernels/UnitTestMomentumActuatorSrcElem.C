/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernels/UnitTestKernelUtils.h"
#include "UnitTestUtils.h"
#include "UnitTestHelperObjects.h"

#include "MomentumActuatorSrcElemKernel.h"

#include <random>

TEST_F(ActuatorSourceKernelHex8Mesh, actuator_source)
{
  fill_mesh_and_init_fields();

  // Setup solution options for default advection kernel
  solnOpts_.meshMotion_ = false;
  solnOpts_.meshDeformation_ = false;
  solnOpts_.externalMeshDeformation_ = false;

  unit_test_utils::HelperObjectsNewME helperObjs(bulk_, stk::topology::HEX_8, 3, partVec_[0]);

  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> kernel(
    new sierra::nalu::MomentumActuatorSrcElemKernel<sierra::nalu::AlgTraitsHex8>(
      bulk_, solnOpts_, helperObjs.assembleElemSolverAlg->dataNeededByKernels_));

  // Add to kernels to be tested
  helperObjs.assembleElemSolverAlg->activeKernels_.push_back(kernel.get());

  helperObjs.assembleElemSolverAlg->execute();

  EXPECT_EQ(helperObjs.linsys->lhs_.dimension(0), 24u);
  EXPECT_EQ(helperObjs.linsys->lhs_.dimension(1), 24u);
  EXPECT_EQ(helperObjs.linsys->rhs_.dimension(0), 24u);

  // Exact solution
  std::vector<double> rhsExact(24,0.0);
  for (size_t i=0; i < 8; i += 1) {
    for (size_t j=0; j < 3; j+= 1) {
      rhsExact[i*3+j] = 0.125 * (j+1) ;
    }
  }

  unit_test_kernel_utils::expect_all_near(helperObjs.linsys->rhs_,rhsExact.data());
}

