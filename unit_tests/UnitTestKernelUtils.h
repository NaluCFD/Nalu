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

  ~TestKernelHex8Mesh() {}

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

/** Test Fixture for Continuity Kernels
 */
class ContinuityKernelHex8Mesh : public TestKernelHex8Mesh
{
public:
  ContinuityKernelHex8Mesh()
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

  void fill_mesh_and_init_fields()
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

/** Text fixture for heat conduction equation kernels
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
