/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level NaluUnit      */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef HighOrderSourceQuad_h
#define HighOrderSourceQuad_h

#include <element_promotion/operators/HighOrderOperatorsQuad.h>
#include <element_promotion/operators/CoefficientMatrices.h>
#include <element_promotion/CVFEMTypeDefs.h>

namespace sierra {
namespace nalu {
namespace tensor_assembly {

  template <int poly_order>
  void add_volumetric_source(
    const CVFEMQuadOperators<poly_order> ops,
    const nodal_scalar_view<AlgTraitsQuad<poly_order>> volume_metric,
    const nodal_scalar_view<AlgTraitsQuad<poly_order>> nodal_source,
    nodal_scalar_view<AlgTraitsQuad<poly_order>> rhs)
  {
    constexpr int n1D = AlgTraitsQuad<poly_order>::nodes1D_;
    for (int j = 0; j < n1D; ++j) {
      for (int i = 0; i < n1D; ++i) {
        nodal_source(j,i) *= volume_metric(j,i);
      }
    }
    ops.volume_2D(nodal_source, rhs);
  }

} // namespace TensorAssembly
} // namespace naluUnit
} // namespace Sierra

#endif
