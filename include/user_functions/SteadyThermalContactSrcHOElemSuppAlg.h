/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SteadyThermalContactSrcHOElemSuppAlg_h
#define SteadyThermalContactSrcHOElemSuppAlg_h

#include <element_promotion/CVFEMTypeDefs.h>
#include <element_promotion/operators/HighOrderOperatorsQuad.h>

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <AlgTraits.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <memory>

namespace sierra{
namespace nalu{

class Realm;
class ElementDescription;
class ElemDataRequests;

template<class AlgTraits>
class SteadyThermalContactSrcHOElemSuppAlg final : public SupplementalAlgorithm
{
public:
  SteadyThermalContactSrcHOElemSuppAlg(
    Realm &realm,
    const ElementDescription& desc,
    ElemDataRequests& dataPreReqs);

  virtual ~SteadyThermalContactSrcHOElemSuppAlg() {}

  void element_execute(
    SharedMemView<double **>& lhs,
    SharedMemView<double *>& rhs,
    stk::mesh::Entity element,
    ScratchViews& scratchViews) final;

  VectorFieldType *coordinates_;

  const double a_;
  const double k_;
  const double pi_;

  // fixed scratch space
  CVFEMQuadOperators<AlgTraits::polyOrder_> ops_;
  nodal_scalar_view<AlgTraits, int> v_node_map_{"tensor_product_node_map"};
  nodal_scalar_view<AlgTraits> v_nodalSource_{"v_nodal_source"};
  nodal_scalar_view<AlgTraits> v_vol_{"v_volume_metric"};
  nodal_scalar_view<AlgTraits> v_rhs_{"v_rhs"};
  nodal_vector_view<AlgTraits> v_coords_{"v_coords"};
};

} // namespace nalu
} // namespace Sierra

#endif
