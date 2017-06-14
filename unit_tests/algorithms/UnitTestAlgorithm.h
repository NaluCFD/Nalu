/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef UNITTESTALGORITHM_H
#define UNITTESTALGORITHM_H

#include <gtest/gtest.h>
#include "UnitTestRealm.h"
#include "UnitTestUtils.h"
#include "UnitTestFieldUtils.h"

#include <memory>
#include <cassert>

#ifndef KOKKOS_HAVE_CUDA

class TestAlgorithm: public ::testing::Test
{
public:
  TestAlgorithm()
    : comm_(MPI_COMM_WORLD)
  {
    YAML::Node doc = unit_test_utils::get_default_inputs();
    naluObj_.reset(new unit_test_utils::NaluTest(doc));
  }

  virtual ~TestAlgorithm() {}

  inline sierra::nalu::Realm& create_realm(
    const YAML::Node& realm_node,
    const std::string realm_type = "multi_physics")
  {
    realm_ = &naluObj_->create_realm(realm_node, realm_type);
    return *realm_;
  }

  inline sierra::nalu::Realm& create_realm(
    const std::string realm_type = "multi_physics")
  {
    const YAML::Node realm_node = unit_test_utils::get_realm_default_node();
    realm_ = &naluObj_->create_realm(realm_node, realm_type);
    return *realm_;
  }

  void fill_mesh(const std::string mesh_spec="generated:10x10x10");

  virtual void declare_fields() = 0;

  inline sierra::nalu::Realm& realm() const
  {
    assert(realm_ != nullptr);
    return *realm_;
  }

  inline stk::mesh::MetaData& meta() const
  {
    return realm().meta_data();
  }

  inline stk::mesh::BulkData& bulk() const
  {
    return realm().bulk_data();
  }

  double field_max(const ScalarFieldType & field, stk::mesh::Selector* selector = nullptr);

  double field_min(const ScalarFieldType & field, stk::mesh::Selector* selector = nullptr);

  double field_norm(const ScalarFieldType & field, stk::mesh::Selector* selector = nullptr);

  double calc_vector_norm(const std::vector<double> & vec);

  //! Reference to test Nalu instance used to hold Simulation and Realm
  std::unique_ptr<unit_test_utils::NaluTest> naluObj_;

  //! Reference to realm instance
  sierra::nalu::Realm* realm_{nullptr};

  stk::mesh::Part* meshPart_{nullptr};
  const VectorFieldType* coordinates_{nullptr};
  stk::ParallelMachine comm_;
};

class TestTurbulenceAlgorithm: public TestAlgorithm
{
public:
  TestTurbulenceAlgorithm()
    : TestAlgorithm()
  {}

  virtual ~TestTurbulenceAlgorithm() {}

  virtual void declare_fields();

  virtual void fill_mesh_and_init_fields(
    const std::string mesh_spec="generated:10x10x10");

  ScalarFieldType* density_{nullptr};
  ScalarFieldType* viscosity_{nullptr};
  ScalarFieldType* tke_{nullptr};
  ScalarFieldType* sdr_{nullptr};
  ScalarFieldType* minDistance_{nullptr};
  GenericFieldType* dudx_{nullptr};
  ScalarFieldType* tvisc_{nullptr};
  ScalarFieldType* maxLengthScale_{nullptr};
  ScalarFieldType* fOneBlend_{nullptr};
  ScalarFieldType* evisc_{nullptr};
  ScalarFieldType* dualNodalVolume_{nullptr};
  VectorFieldType* dkdx_{nullptr};
  VectorFieldType* dwdx_{nullptr};
};

#endif /* KOKKOS_HAVE_CUDA */

#endif /* UNITTESTALGORITHM_H */
