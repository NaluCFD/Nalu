/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/SteadyThermalContactSrcHOElemSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>

#include <element_promotion/ElementDescription.h>
#include <element_promotion/operators/HighOrderOperatorsQuad.h>
#include <element_promotion/operators/HighOrderGeometryQuadVolume.h>
#include <element_promotion/operators/HighOrderSourceQuad.h>

#include <BuildTemplates.h>
#include <ScratchViews.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra{
namespace nalu{

template<class AlgTraits>
SteadyThermalContactSrcHOElemSuppAlg<AlgTraits>::SteadyThermalContactSrcHOElemSuppAlg(
  Realm &realm,
  const ElementDescription& desc,
  ElemDataRequests& dataPreReqs)
  : SupplementalAlgorithm(realm),
    a_(1.0),
    k_(1.0),
    pi_(std::acos(-1.0)),
    ops_()
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  for (int j = 0; j < AlgTraits::nodes1D_; ++j) {
    for (int i = 0; i < AlgTraits::nodes1D_; ++i) {
      v_node_map_(j,i) = desc.node_map(i,j);
    }
  }

  dataPreReqs.add_cvfem_volume_me(realm_.get_volume_master_element(AlgTraits::topo_));
  dataPreReqs.add_gathered_nodal_field(*coordinates_, AlgTraits::nDim_);
}
//--------------------------------------------------------------------------
template<class AlgTraits>
void
SteadyThermalContactSrcHOElemSuppAlg<AlgTraits>::element_execute(
  SharedMemView<double **>& /*lhs*/,
  SharedMemView<double *>& rhs,
  stk::mesh::Entity element,
  ScratchViews& scratchViews)
{
  SharedMemView<double**>& v_flatCoords = scratchViews.get_scratch_view_2D(*coordinates_);

  for (int j = 0; j < AlgTraits::nodes1D_; ++j) {
    for (int i = 0; i < AlgTraits::nodes1D_; ++i) {
      int nodeId = v_node_map_(j,i);
      double x = v_flatCoords(nodeId, 0);
      double y = v_flatCoords(nodeId, 1);

      v_coords_(0,j,i) = x;
      v_coords_(1,j,i) = y;
      v_nodalSource_(j,i) = k_/4.0*(2.0*a_*pi_)*(2.0*a_*pi_)*(cos(2.0*a_*pi_*x) + cos(2.0*a_*pi_*y));
    }
  }

  Kokkos::deep_copy(v_vol_, 0.0);
  high_order_metrics::compute_volume_metric_linear(ops_, v_coords_, v_vol_);

  Kokkos::deep_copy(v_rhs_, 0.0);
  tensor_assembly::add_volumetric_source(ops_, v_vol_, v_nodalSource_, v_rhs_);

  for (int j = 0; j < AlgTraits::nodes1D_; ++j) {
    for (int i = 0; i < AlgTraits::nodes1D_; ++i) {
      rhs(v_node_map_(j,i)) += v_rhs_(j,i);
    }
  }
}

INSTANTIATE_HOQUAD_ALGORITHM(SteadyThermalContactSrcHOElemSuppAlg);

} // namespace nalu
} // namespace Sierra
