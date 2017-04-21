/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ScalarDiffHOElemSuppAlg_h
#define ScalarDiffHOElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <AlgTraits.h>

#include <element_promotion/operators/HighOrderOperatorsQuad.h>
#include <element_promotion/CVFEMTypeDefs.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

// Kokkos
#include <Kokkos_Core.hpp>

namespace sierra{
namespace nalu{

class ElemDataRequests;
class Realm;
class ScratchViews;
class MasterElement;

template<class AlgTraits>
class ScalarDiffHOElemSuppAlg final: public SupplementalAlgorithm
{
public:
  ScalarDiffHOElemSuppAlg(
    Realm &realm,
    ScalarFieldType *scalarQ,
    ScalarFieldType *diffFluxCoeff,
    const ElementDescription& desc,
    ElemDataRequests& dataPreReqs);

  void element_execute(
    SharedMemView<double **>& lhs,
    SharedMemView<double *>& rhs,
    stk::mesh::Entity element,
    ScratchViews& scratchViews) final;

private:
  ScalarFieldType *scalarQ_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *coordinates_;

  CVFEMQuadOperators<AlgTraits::polyOrder_> ops_{};
  Kokkos::View<int[AlgTraits::nodesPerElement_]> v_node_map_{"tensor_product_node_map"};
  nodal_scalar_view<AlgTraits> v_rhs_{"scalar_ho_diff_rhs"};
  nodal_scalar_view<AlgTraits> v_diff_{"scalar_ho_diffusion_coeff"};
  nodal_scalar_view<AlgTraits> v_scalar_{"scalar_ho_q"};
  nodal_vector_view<AlgTraits> v_coords_{"coords"};
  scs_tensor_view<AlgTraits> v_metric_{"metric"};
  matrix_view<AlgTraits> v_lhs_{"scalar_ho_diff_lhs"};
};

} // namespace nalu
} // namespace Sierra

#endif
