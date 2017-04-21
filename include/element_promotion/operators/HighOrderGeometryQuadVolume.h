/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level NaluUnit      */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef HighOrderGeometryQuadVolume_h
#define HighOrderGeometryQuadVolume_h

#include <element_promotion/operators/CoefficientMatrices.h>
#include <element_promotion/operators/DirectionEnums.h>
#include <element_promotion/CVFEMTypeDefs.h>
#include <AlgTraits.h>

#include <stk_util/environment/ReportHandler.hpp>

namespace sierra {
namespace nalu {
namespace high_order_metrics
{

  template <int p>
  void compute_volume_metric_linear(
    const CVFEMQuadOperators<p> ops,
    const nodal_vector_view<AlgTraitsQuad<p>>& coordinates,
    nodal_scalar_view<AlgTraitsQuad<p>>& vol)
  {
    // Computes det(J) at nodes using a linear basis for element geometry
    const double dx_x0 = coordinates(XH, p, 0) - coordinates(XH, 0, 0);
    const double dx_x1 = coordinates(XH, 0, p) - coordinates(XH, 0, 0);
    const double dx_y0 = coordinates(XH, p, p) - coordinates(XH, p, 0);
    const double dx_y1 = coordinates(XH, p, p) - coordinates(XH, 0, p);

    const double dy_x0 = coordinates(YH, p, 0) - coordinates(YH, 0, 0);
    const double dy_x1 = coordinates(YH, 0, p) - coordinates(YH, 0, 0);
    const double dy_y0 = coordinates(YH, p, p) - coordinates(YH, p, 0);
    const double dy_y1 = coordinates(YH, p, p) - coordinates(YH, 0, p);

    const auto& mat = ops.mat_;

    for (int j = 0; j < p + 1; ++j) {
      const double dx_dyh = mat.linear_nodal_interp(0,j) * dx_x1 + mat.linear_nodal_interp(1,j) * dx_y0;
      const double dy_dyh = mat.linear_nodal_interp(0,j) * dy_x1 + mat.linear_nodal_interp(1,j) * dy_y0;

      for (int i = 0; i < p + 1; ++i) {
        const double dx_dxh = mat.linear_nodal_interp(0,i) * dx_x0 + mat.linear_nodal_interp(1,i) * dx_y1;
        const double dy_dxh = mat.linear_nodal_interp(0,i) * dy_x0 + mat.linear_nodal_interp(1,i) * dy_y1;

        // times divided by 4 for missing factor of a half in the derivatives
        vol(j,i) = 0.25 * (dx_dyh * dy_dxh  - dx_dxh * dy_dyh);
      }
    }
  }

} // namespace HighOrderGeometryQuad
} // namespace naluUnit
} // namespace Sierra

#endif
