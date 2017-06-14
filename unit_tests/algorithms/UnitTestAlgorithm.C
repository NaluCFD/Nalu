/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "UnitTestAlgorithm.h"
#include "UnitTestUtils.h"
#include "UnitTestKokkosUtils.h"
#include "kernels/UnitTestKernelUtils.h"

#include "Realm.h"

void
TestAlgorithm::fill_mesh(const std::string mesh_spec)
{
  this->declare_fields();

  unit_test_utils::fill_hex8_mesh(mesh_spec, bulk());
  meshPart_ = meta().get_part("block_1");
  coordinates_ = static_cast<const VectorFieldType*>(meta().coordinate_field());
  EXPECT_TRUE(coordinates_ != nullptr);
}

double
TestAlgorithm::field_max(const ScalarFieldType & field, stk::mesh::Selector* selector)
{
  auto& meta = this->meta();
  auto& bulk = this->bulk();
  auto sel = (selector == nullptr)? meta.locally_owned_part() : *selector;

  return unit_test_utils::field_max(field, bulk, sel);
}

double
TestAlgorithm::field_min(const ScalarFieldType & field, stk::mesh::Selector* selector)
{
  auto& meta = this->meta();
  auto& bulk = this->bulk();
  auto sel = (selector == nullptr)? meta.locally_owned_part() : *selector;

  return unit_test_utils::field_min(field, bulk, sel);
}

double
TestAlgorithm::field_norm(const ScalarFieldType & field, stk::mesh::Selector* selector)
{

  auto& meta = this->meta();
  auto& bulk = this->bulk();
  auto sel = (selector == nullptr)? meta.locally_owned_part() : *selector;

  return unit_test_utils::field_norm(field, bulk, sel);
}

double
TestAlgorithm::calc_vector_norm(const std::vector<double> & vec)
{
  size_t N = vec.size();
  size_t g_N = 0;
  double norm = 0.0;
  double g_norm = 0.0;

  for (int i = 0; i < N; ++i) {
    norm += vec[i]*vec[i];
  }
  stk::all_reduce_sum(comm_, &N, &g_N, 1);
  stk::all_reduce_sum(comm_, &norm, &g_norm, 1);
  g_norm = std::sqrt(g_norm/g_N);

  return g_norm;
}

void
TestTurbulenceAlgorithm::declare_fields()
{
  auto& meta = this->meta();
  auto spatialDim = meta.spatial_dimension();

  density_ = (
    &meta.declare_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "density"));
  viscosity_ = (
    &meta.declare_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "viscosity"));
  tke_ = (
    &meta.declare_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "turbulent_ke"));
  sdr_ = (
    &meta.declare_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "specific_dissipation_rate"));
  minDistance_ = (
    &meta.declare_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "minimum_distance_to_wall"));
  dudx_ = (
    &meta.declare_field<GenericFieldType>(
      stk::topology::NODE_RANK, "dudx"));
  tvisc_ = (
    &meta.declare_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "turbulent_viscosity"));
  maxLengthScale_ = (
    &meta.declare_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "sst_max_length_scale"));
  fOneBlend_ = (
    &meta.declare_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "sst_f_one_blending"));
  evisc_ = (
     &meta.declare_field<ScalarFieldType>(
       stk::topology::NODE_RANK, "effective_viscosity"));
  evisc_ = (
     &meta.declare_field<ScalarFieldType>(
       stk::topology::NODE_RANK, "effective_viscosity"));
  dualNodalVolume_ = (
     &meta.declare_field<ScalarFieldType>(
       stk::topology::NODE_RANK, "dual_nodal_volume"));
  dkdx_ = (
     &meta.declare_field<VectorFieldType>(
       stk::topology::NODE_RANK, "dkdx"));
  dwdx_ = (
     &meta.declare_field<VectorFieldType>(
       stk::topology::NODE_RANK, "dwdx"));

  stk::mesh::put_field(*density_, meta.universal_part(), 1);
  stk::mesh::put_field(*viscosity_, meta.universal_part(), 1);
  stk::mesh::put_field(*tke_, meta.universal_part(), 1);
  stk::mesh::put_field(*sdr_, meta.universal_part(), 1);
  stk::mesh::put_field(*minDistance_, meta.universal_part(), 1);
  stk::mesh::put_field(*dudx_, meta.universal_part(), spatialDim*spatialDim);
  stk::mesh::put_field(*tvisc_, meta.universal_part(), 1);
  stk::mesh::put_field(*maxLengthScale_, meta.universal_part(), 1);
  stk::mesh::put_field(*fOneBlend_, meta.universal_part(), 1);
  stk::mesh::put_field(*evisc_, meta.universal_part(), 1);
  stk::mesh::put_field(*dualNodalVolume_, meta.universal_part(), 1);
  stk::mesh::put_field(*dkdx_, meta.universal_part(), spatialDim);
  stk::mesh::put_field(*dwdx_, meta.universal_part(), spatialDim);
}

void
TestTurbulenceAlgorithm::fill_mesh_and_init_fields(const std::string mesh_spec)
{
  fill_mesh(mesh_spec);

  auto& bulk = this->bulk();

  unit_test_kernel_utils::density_test_function(bulk, *coordinates_, *density_);
  stk::mesh::field_fill(0.2, *viscosity_);
  unit_test_kernel_utils::tke_test_function(bulk, *coordinates_, *tke_);
  unit_test_kernel_utils::sdr_test_function(bulk, *coordinates_, *sdr_);
  unit_test_kernel_utils::minimum_distance_to_wall_test_function(bulk, *coordinates_, *minDistance_);
  unit_test_kernel_utils::dudx_test_function(bulk, *coordinates_, *dudx_);
  unit_test_kernel_utils::turbulent_viscosity_test_function(bulk, *coordinates_, *tvisc_);
  stk::mesh::field_fill(0.5, *maxLengthScale_);
  unit_test_kernel_utils::sst_f_one_blending_test_function(bulk, *coordinates_, *fOneBlend_);
  stk::mesh::field_fill(0.0, *evisc_);
  stk::mesh::field_fill(0.2, *dualNodalVolume_);
  unit_test_kernel_utils::dkdx_test_function(bulk, *coordinates_, *dkdx_);
  unit_test_kernel_utils::dwdx_test_function(bulk, *coordinates_, *dwdx_);
}
