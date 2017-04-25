/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "UnitTestKernelUtils.h"
#include "UnitTestUtils.h"
#include "UnitTestKernelGolds.h"

#include "Kernel.h"
#include "ElemDataRequests.h"
#include "ScratchViews.h"
#include "AlgTraits.h"
#include "KokkosInterface.h"
#include "TimeIntegrator.h"

// Computational Kernels to be tested
#include "ContinuityMassElemKernel.h"
#include "ContinuityAdvElemKernel.h"
#include "ScalarDiffElemKernel.h"
#include "user_functions/SteadyThermal3dContactSrcElemKernel.h"

#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>

namespace {

void expect_all_near(
  const sierra::nalu::SharedMemView<double*>& calcValue,
  const double* exactValue,
  const double tol = 1.0e-15)
{
  const auto length = calcValue.dimension(0);

  for (int i=0; i < length; ++i) {
    EXPECT_NEAR(calcValue[i], exactValue[i], tol);
  }
}

template<int N>
void expect_all_near(
  const sierra::nalu::SharedMemView<double**>& calcValue,
  const double (*exactValue)[N],
  const double tol = 1.0e-15)
{
  const auto dim1 = calcValue.dimension(0);
  const auto dim2 = calcValue.dimension(1);
  EXPECT_EQ(dim2, N);

  for (int i=0; i < dim1; ++i) {
    for (int j=0; j < dim2; ++j) {
      EXPECT_NEAR(calcValue(i,j), exactValue[i][j], tol);
    }
  }
}

void expect_all_near(
  const sierra::nalu::SharedMemView<double*>& calcValue,
  const double exactValue,
  const double tol = 1.0e-15)
{
  const auto length = calcValue.dimension(0);

  for (int i=0; i < length; ++i) {
    EXPECT_NEAR(calcValue[i], exactValue, tol);
  }
}

template<int N>
void expect_all_near(
  const sierra::nalu::SharedMemView<double**>& calcValue,
  const double exactValue,
  const double tol = 1.0e-15)
{
  const auto dim1 = calcValue.dimension(0);
  const auto dim2 = calcValue.dimension(1);
  EXPECT_EQ(dim2, N);

  for (int i=0; i < dim1; ++i) {
    for (int j=0; j < dim2; ++j) {
      EXPECT_NEAR(calcValue(i,j), exactValue, tol);
    }
  }
}

/** Driver class that mimics Assemble*SolverAlgorithm
 *
 * It is the caller's responsibility to populate `activeKernels_` and call the
 * `setup` method on the activated Kernels before calling the `execute` method
 * of this class.
 *
 * The execute method will assert that it is being called with a one-element
 * mesh when it loops over the buckets.
 */
class TestKernelDriver
{
public:
  TestKernelDriver(
    const stk::mesh::BulkData& bulk,
    const stk::mesh::PartVector& partVec,
    const VectorFieldType* coordinates,
    const int numDof = 1,
    const stk::topology topo = stk::topology::HEX_8)
    : bulk_(bulk),
      partVec_(partVec),
      coordinates_(coordinates),
      topo_(topo),
      numDof_(numDof)
  {}

  void execute()
  {
    const stk::mesh::MetaData& meta = bulk_.mesh_meta_data();

    stk::mesh::Selector s_locally_owned_union = (
      meta.locally_owned_part() & stk::mesh::selectUnion(partVec_));

    const auto& buckets = bulk_.get_buckets(stk::topology::ELEM_RANK,
                                            s_locally_owned_union);

    EXPECT_EQ(buckets.size(), 1u);

    const int rhsSize = topo_.num_nodes() * numDof_;
    const int lhsSize = rhsSize * rhsSize;
    const int bytes_per_team = 0;
    const int num_bytes_for_kernels = sierra::nalu::get_num_bytes_pre_req_data(
        dataNeededByKernels_, meta.spatial_dimension()) ;
    const int bytes_per_thread = (rhsSize + lhsSize) * sizeof(double) +
      num_bytes_for_kernels;

    auto team_exec = sierra::nalu::get_team_policy(
      buckets.size(), bytes_per_team, bytes_per_thread);

    Kokkos::parallel_for(team_exec, [&](const sierra::nalu::TeamHandleType& team) {
        auto& b = *buckets[team.league_rank()];
        const auto length = b.size();
        EXPECT_EQ(length, 1u);

        sierra::nalu::ScratchViews preReqData(
          team, bulk_, topo_, dataNeededByKernels_);
        rhs_ = sierra::nalu::get_shmem_view_1D(team, rhsSize);
        lhs_ = sierra::nalu::get_shmem_view_2D(team, rhsSize, rhsSize);

        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(team, length), [&](const size_t& k) {
            stk::mesh::Entity element = b[k];
            sierra::nalu::fill_pre_req_data(
              dataNeededByKernels_, bulk_, topo_, element, coordinates_, preReqData);

            for (int i=0; i < rhsSize; i++) {
              rhs_(i) = 0.0;
              for (int j=0; j < rhsSize; j++) {
                lhs_(i,j) = 0.0;
              }
            }

            for (size_t i=0; i < activeKernels_.size(); ++i)
              activeKernels_[i]->execute(lhs_, rhs_, element, preReqData);
          });

      });
  }

  /** Convenience function to dump LHS and RHS
   *
   * Used to generate the gold values as well as for debugging
   */
  void dump_lhs_and_rhs()
  {
    using unit_test_utils::nalu_out;
    const int rhsSize = rhs_.dimension(0);

    // Dump the RHS
    nalu_out() << std::endl
               << "static constexpr double rhs[" << rhsSize << "] = {"
               << std::endl << "  ";
    for (int i=0; i < rhsSize; i++) {
      nalu_out() << std::setprecision(12) << rhs_(i) << ", ";
    }
    nalu_out() << "};" << std::endl;

    // Dump the LHS
    nalu_out() << std::endl
               << "static constexpr double lhs[" << rhsSize << "]["
               << rhsSize << "] = {" << std::endl;
    for (int i=0; i < rhsSize; i++) {
      nalu_out() << "  { ";
      for (int j=0; j < rhsSize; j++) {
        nalu_out() << std::setprecision(12) << lhs_(i,j) << ", ";
      }
      std::cerr << " }," << std::endl;
    }
    nalu_out() << "};" << std::endl << std::endl;
  }

  std::vector<sierra::nalu::Kernel*> activeKernels_;
  sierra::nalu::ElemDataRequests dataNeededByKernels_;
  sierra::nalu::SharedMemView<double*> rhs_;
  sierra::nalu::SharedMemView<double**> lhs_;

private:
  const stk::mesh::BulkData& bulk_;
  const stk::mesh::PartVector& partVec_;
  const VectorFieldType* coordinates_;
  const stk::topology topo_{stk::topology::HEX_8};
  const int numDof_;
};

} // anonymous namespace

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
  TestKernelDriver assembleKernels(bulk_, partVec_, coordinates_,
                                   1, stk::topology::HEX_8);

  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> advKernel(
    new sierra::nalu::ContinuityAdvElemKernel<sierra::nalu::AlgTraitsHex8>(
      bulk_, solnOpts_, assembleKernels.dataNeededByKernels_));

  // Register the kernel for execution
  assembleKernels.activeKernels_.push_back(advKernel.get());

  // Populate LHS and RHS
  assembleKernels.execute();

  namespace gold_values = unit_test_golds::hex8_golds::continuity::advection_default;

  expect_all_near(assembleKernels.rhs_, gold_values::rhs);
  expect_all_near<8>(assembleKernels.lhs_, gold_values::lhs);
}

TEST_F(ContinuityKernelHex8Mesh, advection_reduced_sensitivities)
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
  TestKernelDriver assembleKernels(bulk_, partVec_, coordinates_,
                                   1, stk::topology::HEX_8);

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

  expect_all_near(assembleKernels.rhs_, gold_values::rhs);
  expect_all_near<8>(assembleKernels.lhs_, gold_values::lhs);
}

TEST_F(ContinuityKernelHex8Mesh, density_time_derivative)
{
  fill_mesh_and_init_fields();

  // Setup solution options for default advection kernel
  solnOpts_.meshMotion_ = false;
  solnOpts_.meshDeformation_ = false;
  solnOpts_.externalMeshDeformation_ = false;

  // Initialize the kernel driver
  TestKernelDriver assembleKernels(bulk_, partVec_, coordinates_,
                                   1, stk::topology::HEX_8);

  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> massKernel(
    new sierra::nalu::ContinuityMassElemKernel<sierra::nalu::AlgTraitsHex8>(
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

  expect_all_near(assembleKernels.rhs_,-12.5);
  expect_all_near<8>(assembleKernels.lhs_,0.0);
}

TEST_F(ContinuityKernelHex8Mesh, density_time_derivative_lumped)
{
  fill_mesh_and_init_fields();

  // Setup solution options for default advection kernel
  solnOpts_.meshMotion_ = false;
  solnOpts_.meshDeformation_ = false;
  solnOpts_.externalMeshDeformation_ = false;

  // Initialize the kernel driver
  TestKernelDriver assembleKernels(bulk_, partVec_, coordinates_,
                                   1, stk::topology::HEX_8);

  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> massKernel(
    new sierra::nalu::ContinuityMassElemKernel<sierra::nalu::AlgTraitsHex8>(
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

  expect_all_near(assembleKernels.rhs_,-12.5);
  expect_all_near<8>(assembleKernels.lhs_,0.0);
}

TEST_F(HeatCondKernelHex8Mesh, steady_3d_thermal)
{
  fill_mesh_and_init_fields();

  // Setup solution options for default advection kernel
  solnOpts_.meshMotion_ = false;
  solnOpts_.meshDeformation_ = false;
  solnOpts_.externalMeshDeformation_ = false;

  // Initialize the kernel driver
  TestKernelDriver assembleKernels(bulk_, partVec_, coordinates_,
                                   1, stk::topology::HEX_8);

  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> kernel(
    new sierra::nalu::SteadyThermal3dContactSrcElemKernel<sierra::nalu::AlgTraitsHex8>(
      bulk_, solnOpts_, assembleKernels.dataNeededByKernels_));

  // Add to kernels to be tested
  assembleKernels.activeKernels_.push_back(kernel.get());

  assembleKernels.execute();

  namespace gold_values = unit_test_golds::hex8_golds::heatcond::steady_3d_thermal;

  expect_all_near(assembleKernels.rhs_, gold_values::rhs);
  expect_all_near<8>(assembleKernels.lhs_, 0.0);
}

TEST_F(HeatCondKernelHex8Mesh, cvfem_diff)
{
  fill_mesh_and_init_fields();

  // Setup solution options for default advection kernel
  solnOpts_.meshMotion_ = false;
  solnOpts_.meshDeformation_ = false;
  solnOpts_.externalMeshDeformation_ = false;

  // Initialize the kernel driver
  TestKernelDriver assembleKernels(bulk_, partVec_, coordinates_,
                                   1, stk::topology::HEX_8);

  // Initialize the kernel
  std::unique_ptr<sierra::nalu::Kernel> kernel(
    new sierra::nalu::ScalarDiffElemKernel<sierra::nalu::AlgTraitsHex8>(
      bulk_, solnOpts_, temperature_, thermalCond_,
      assembleKernels.dataNeededByKernels_));

  // Add to kernels to be tested
  assembleKernels.activeKernels_.push_back(kernel.get());

  assembleKernels.execute();

  namespace gold_values = unit_test_golds::hex8_golds::heatcond::cvfem_diff;
  expect_all_near(assembleKernels.rhs_, 0.0);
  expect_all_near<8>(assembleKernels.lhs_, gold_values::lhs);
}
