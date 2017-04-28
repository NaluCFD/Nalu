/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef UNITTESTKERNELUTILS_H
#define UNITTESTKERNELUTILS_H

#include <gtest/gtest.h>
#include "UnitTestUtils.h"

#include "SolutionOptions.h"
#include "Kernel.h"
#include "ElemDataRequests.h"
#include "ScratchViews.h"
#include "AlgTraits.h"
#include "KokkosInterface.h"
#include "TimeIntegrator.h"

#include <gtest/gtest.h>

#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace unit_test_kernel_utils {

void velocity_test_function(
  const stk::mesh::BulkData&,
  const VectorFieldType& coordinates,
  VectorFieldType& velocity);

void pressure_test_function(
  const stk::mesh::BulkData&,
  const VectorFieldType& coordinates,
  ScalarFieldType& pressure);

void dpdx_test_function(
  const stk::mesh::BulkData&,
  const VectorFieldType& coordinates,
  VectorFieldType& dpdx);

void temperature_test_function(
  const stk::mesh::BulkData& bulk,
  const VectorFieldType& coordinates,
  ScalarFieldType& temperature);

void calc_mass_flow_rate_scs(
  const stk::mesh::BulkData&,
  const stk::topology&,
  const VectorFieldType&,
  const ScalarFieldType&,
  const VectorFieldType&,
  const GenericFieldType&);

void expect_all_near(
  const sierra::nalu::SharedMemView<double*>& calcValue,
  const double* exactValue,
  const double tol = 1.0e-15);

void expect_all_near(
  const sierra::nalu::SharedMemView<double*>& calcValue,
  const double exactValue,
  const double tol = 1.0e-15);

void expect_all_near(
  const sierra::nalu::SharedMemView<double**>& calcValue,
  const double* exactValue,
  const double tol = 1.0e-15);

template<int N>
void expect_all_near(
  const sierra::nalu::SharedMemView<double**>& calcValue,
  const double (*exactValue)[N],
  const double tol = 1.0e-15)
{
  const int dim1 = calcValue.dimension(0);
  const int dim2 = calcValue.dimension(1);
  EXPECT_EQ(dim2, N);

  for (int i=0; i < dim1; ++i) {
    for (int j=0; j < dim2; ++j) {
      EXPECT_NEAR(calcValue(i,j), exactValue[i][j], tol);
    }
  }
}

template<int N>
void expect_all_near(
  const sierra::nalu::SharedMemView<double**>& calcValue,
  const double exactValue,
  const double tol = 1.0e-15)
{
  const int dim1 = calcValue.dimension(0);
  const int dim2 = calcValue.dimension(1);
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

    // For LHS/RHS golds expect only one element in the mesh
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
        // For LHS/RHS golds expect only one element in the mesh
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
  void dump_lhs_and_rhs(double tol = 1.0e-15)
  {
    using unit_test_utils::nalu_out;
    const int rhsSize = rhs_.dimension(0);

    // Dump the RHS
    nalu_out() << std::endl
               << "static constexpr double rhs[" << rhsSize << "] = {"
               << std::endl << "  ";
    for (int i=0; i < rhsSize; i++) {
      nalu_out() << std::setprecision(12)
                 << (std::fabs(rhs_(i)) < tol ? 0.0 : rhs_(i)) << ", ";
    }
    nalu_out() << "};" << std::endl;

    // Dump the LHS
    nalu_out() << std::endl
               << "static constexpr double lhs[" << rhsSize << "]["
               << rhsSize << "] = {" << std::endl;
    for (int i=0; i < rhsSize; i++) {
      nalu_out() << "  { ";
      for (int j=0; j < rhsSize; j++) {
        nalu_out() << std::setprecision(12)
                   << (std::fabs(lhs_(i,j)) < tol ? 0.0 : lhs_(i,j)) << ", ";
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

}

#ifndef KOKKOS_HAVE_CUDA

/** Base class for all computational kernel testing setups
 *
 *  Initializes the STK mesh data structures and stores a handle to the
 *  coordinates field. Subclasses must declare the necessary computational
 *  fields necessary for the tests.
 */
class TestKernelHex8Mesh : public ::testing::Test
{
public:
  TestKernelHex8Mesh()
    : comm_(MPI_COMM_WORLD),
      spatialDim_(3),
      meta_(spatialDim_),
      bulk_(meta_, comm_),
      solnOpts_()
  {}

  virtual ~TestKernelHex8Mesh() {}

  void fill_mesh()
  {
    unit_test_utils::fill_mesh_1_elem_per_proc_hex8(bulk_);

    partVec_ = {meta_.get_part("block_1")};

    coordinates_ = static_cast<const VectorFieldType*>(meta_.coordinate_field());

    EXPECT_TRUE(coordinates_ != nullptr);
  }

  stk::ParallelMachine comm_;
  unsigned spatialDim_;
  stk::mesh::MetaData meta_;
  stk::mesh::BulkData bulk_;
  stk::mesh::PartVector partVec_;

  sierra::nalu::SolutionOptions solnOpts_;

  const VectorFieldType* coordinates_{nullptr};
};

/** Test Fixture for Low-Mach Kernels
 *
 *  This test fixture performs the following actions:
 *    - Create a HEX8 mesh with one element
 *    - Declare `velocity`, `pressure`, `density`, and `dpdx` fields
 *    - Initialize the fields with SteadyTaylorVortex solution
 *    - `density` is initialized to 1.0
 */
class LowMachKernelHex8Mesh : public TestKernelHex8Mesh
{
public:
  LowMachKernelHex8Mesh()
    : TestKernelHex8Mesh(),
      velocity_(
        &meta_.declare_field<VectorFieldType>(
          stk::topology::NODE_RANK, "velocity",2)),
      dpdx_(
        &meta_.declare_field<VectorFieldType>(
          stk::topology::NODE_RANK, "dpdx",2)),
      density_(
        &meta_.declare_field<ScalarFieldType>(
          stk::topology::NODE_RANK, "density",2)),
      pressure_(
        &meta_.declare_field<ScalarFieldType>(
          stk::topology::NODE_RANK, "pressure",2))
  {
    stk::mesh::put_field(*velocity_, meta_.universal_part(), spatialDim_);
    stk::mesh::put_field(*dpdx_, meta_.universal_part(), spatialDim_);
    stk::mesh::put_field(*density_, meta_.universal_part(), 1);
    stk::mesh::put_field(*pressure_, meta_.universal_part(), 1);
  }

  virtual ~LowMachKernelHex8Mesh() {}

  virtual void fill_mesh_and_init_fields()
  {
    fill_mesh();

    unit_test_kernel_utils::velocity_test_function(bulk_, *coordinates_, *velocity_);
    unit_test_kernel_utils::pressure_test_function(bulk_, *coordinates_, *pressure_);
    unit_test_kernel_utils::dpdx_test_function(bulk_, *coordinates_, *dpdx_);
    stk::mesh::field_fill(1.0, *density_);
  }

  VectorFieldType* velocity_{nullptr};
  VectorFieldType* dpdx_{nullptr};
  ScalarFieldType* density_{nullptr};
  ScalarFieldType* pressure_{nullptr};
};

class ContinuityKernelHex8Mesh : public LowMachKernelHex8Mesh
{
public:
  virtual ~ContinuityKernelHex8Mesh() {}
};

class MomentumKernelHex8Mesh : public LowMachKernelHex8Mesh
{
public:
  MomentumKernelHex8Mesh()
    : LowMachKernelHex8Mesh(),
      massFlowRate_(
        &meta_.declare_field<GenericFieldType>(
          stk::topology::ELEM_RANK, "mass_flow_rate_scs")),
      viscosity_(
        &meta_.declare_field<ScalarFieldType>(
          stk::topology::NODE_RANK, "viscosity"))
  {
    const auto& meSCS = sierra::nalu::get_surface_master_element(stk::topology::HEX_8);
    stk::mesh::put_field(*massFlowRate_, meta_.universal_part(), meSCS->numIntPoints_);
    stk::mesh::put_field(*viscosity_, meta_.universal_part(), 1);
  }

  virtual ~MomentumKernelHex8Mesh() {}

  virtual void fill_mesh_and_init_fields()
  {
    LowMachKernelHex8Mesh::fill_mesh_and_init_fields();
    unit_test_kernel_utils::calc_mass_flow_rate_scs(
      bulk_, stk::topology::HEX_8, *coordinates_, *density_, *velocity_, *massFlowRate_);
    stk::mesh::field_fill(0.1, *viscosity_);
  }

  GenericFieldType* massFlowRate_{nullptr};
  ScalarFieldType* viscosity_{nullptr};
};

/** Text fixture for heat conduction equation kernels
 *
 *  This test fixture performs the following actions:
 *    - Create a HEX8 mesh with one element
 *    - Declare `temperature` and `thermal_conductivity` fields
 *    - Initialize the fields with steady 3-D thermal solution
 *    - `thermal_conductivity` is initialized to 1.0
 */
class HeatCondKernelHex8Mesh : public TestKernelHex8Mesh
{
public:
  HeatCondKernelHex8Mesh()
    : TestKernelHex8Mesh(),
      temperature_(
        &meta_.declare_field<ScalarFieldType>(
          stk::topology::NODE_RANK, "temperature",2)),
      thermalCond_(
        &meta_.declare_field<ScalarFieldType>(
          stk::topology::NODE_RANK, "thermal_conductivity",2))
  {
    stk::mesh::put_field(*temperature_, meta_.universal_part(), 1);
    stk::mesh::put_field(*thermalCond_, meta_.universal_part(), 1);
  }

  void fill_mesh_and_init_fields()
  {
    fill_mesh();

    unit_test_kernel_utils::temperature_test_function(bulk_, *coordinates_, *temperature_);
    stk::mesh::field_fill(1.0, *thermalCond_);
  }

  ScalarFieldType* temperature_{nullptr};
  ScalarFieldType* thermalCond_{nullptr};
};

#endif /* KOKKOS_HAVE_CUDA */

#endif /* UNITTESTKERNELUTILS_H */
