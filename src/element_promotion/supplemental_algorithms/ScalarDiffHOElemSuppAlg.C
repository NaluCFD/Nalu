/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <element_promotion/supplemental_algorithms/ScalarDiffHOElemSuppAlg.h>

#include <element_promotion/operators/HighOrderDiffusionQuad.h>
#include <element_promotion/operators/HighOrderGeometryQuadDiffusion.h>
#include <element_promotion/operators/CoefficientMatrices.h>
#include <element_promotion/operators/DirectionEnums.h>
#include <element_promotion/ElementDescription.h>

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>

// template and scratch space
#include <BuildTemplates.h>
#include <ScratchViews.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>

// topology
#include <stk_topology/topology.hpp>

// Kokkos
#include <Kokkos_Core.hpp>

namespace sierra{
namespace nalu{

template<class AlgTraits>
ScalarDiffHOElemSuppAlg<AlgTraits>::ScalarDiffHOElemSuppAlg(
  Realm &realm,
  ScalarFieldType *scalarQ,
  ScalarFieldType *diffFluxCoeff,
  const ElementDescription& desc,
  ElemDataRequests& dataPreReqs)
  : SupplementalAlgorithm(realm),
    scalarQ_(scalarQ),
    diffFluxCoeff_(diffFluxCoeff)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // map from "exodus-style" node ordering to usual tensor-product based node ordering.
  for (int j = 0; j < AlgTraits::nodes1D_; ++j) {
    for (int i = 0; i < AlgTraits::nodes1D_; ++i) {
      v_node_map_(j*AlgTraits::nodes1D_+i) = desc.node_map(i,j);
    }
  }

  // only necessary for correctly sizing scratch views
  dataPreReqs.add_cvfem_surface_me(get_surface_master_element(AlgTraits::topo_));

  dataPreReqs.add_gathered_nodal_field(*coordinates_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*scalarQ, 1);
  dataPreReqs.add_gathered_nodal_field(*diffFluxCoeff, 1);
}
//--------------------------------------------------------------------------
template<class AlgTraits> void
ScalarDiffHOElemSuppAlg<AlgTraits>::element_execute(
  SharedMemView<double **>& lhs,
  SharedMemView<double *>& rhs,
  stk::mesh::Entity element,
  ScratchViews& scratchViews)
{
  SharedMemView<double**> v_flatCoords = scratchViews.get_scratch_view_2D(*coordinates_);
  SharedMemView<double*> v_flatScalar = scratchViews.get_scratch_view_1D(*scalarQ_);
  SharedMemView<double*> v_flatDiff  = scratchViews.get_scratch_view_1D(*diffFluxCoeff_);

  // reorder fields into the ordering expected by the alg
  for (int j = 0; j < AlgTraits::nodes1D_; ++j) {
    for (int i = 0; i < AlgTraits::nodes1D_; ++i) {
      int nodeId = v_node_map_(j*AlgTraits::nodes1D_+i);
      v_scalar_(j,i) = v_flatScalar(nodeId);
      v_diff_(j,i)   = v_flatDiff(nodeId);
      for (int d = 0; d < AlgTraits::nDim_; ++d) {
        v_coords_(d,j,i) = v_flatCoords(nodeId, d);
      }
    }
  }

  Kokkos::deep_copy(v_lhs_, 0.0);
  Kokkos::deep_copy(v_rhs_, 0.0);

  /*
   * The computation of the diffusion term, split into three steps.
   * first, compute the metric J^-T Diffusivity A.  Then use it to compute the left-hand side
   * through nested loops, and the right-hand side through a sequence of calls to dgemm.
   * In principle, if we're computing the Jacobian, then we should just form
   * the rhs along the way.  Instead, we have more efficient
   * rhs-alone calculation (elemental_diffusion_action) done separately
   * mainly for experimentation
   */
  high_order_metrics::compute_diffusion_metric_linear(ops_, v_coords_, v_diff_, v_metric_);
  tensor_assembly::elemental_diffusion_jacobian(ops_, v_metric_, v_lhs_);
  tensor_assembly::elemental_diffusion_action(ops_, v_metric_, v_scalar_, v_rhs_);


  // map lhs/rhs back to the usual ordering
  for (int j = 0; j < AlgTraits::nodesPerElement_; ++j) {
    for (int i = 0; i < AlgTraits::nodesPerElement_; ++i) {
      lhs(v_node_map_(j), v_node_map_(i)) += v_lhs_(j,i);
    }
  }

  for (int j = 0; j < AlgTraits::nodes1D_; ++j) {
    for (int i = 0; i < AlgTraits::nodes1D_; ++i) {
      rhs(v_node_map_(j*AlgTraits::nodes1D_+i)) += v_rhs_(j,i);
    }
  }
}

INSTANTIATE_HOQUAD_ALGORITHM(ScalarDiffHOElemSuppAlg)

} // namespace nalu
} // namespace Sierra
