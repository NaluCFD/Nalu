/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level NaluUnit      */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef HighOrderGeometryQuadDiffusion_h
#define HighOrderGeometryQuadDiffusion_h

#include <element_promotion/operators/CoefficientMatrices.h>
#include <element_promotion/operators/HighOrderOperatorsQuad.h>
#include <element_promotion/operators/DirectionEnums.h>
#include <element_promotion/CVFEMTypeDefs.h>
#include <AlgTraits.h>

#include <stk_util/environment/ReportHandler.hpp>

namespace sierra {
namespace nalu {
namespace high_order_metrics
{
  template <int p>
  void compute_laplacian_metric_linear(
    const CoefficientMatrices<p> mat,
    const nodal_vector_view<AlgTraitsQuad<p>> coordinates,
    scs_tensor_view<AlgTraitsQuad<p>> metric)
  {
    const double dx_x0 = coordinates(XH, p, 0) - coordinates(XH, 0, 0);
    const double dx_x1 = coordinates(XH, 0, p) - coordinates(XH, 0, 0);
    const double dx_y0 = coordinates(XH, p, p) - coordinates(XH, p, 0);
    const double dx_y1 = coordinates(XH, p, p) - coordinates(XH, 0, p);

    const double dy_x0 = coordinates(YH, p, 0) - coordinates(YH, 0, 0);
    const double dy_x1 = coordinates(YH, 0, p) - coordinates(YH, 0, 0);
    const double dy_y0 = coordinates(YH, p, p) - coordinates(YH, p, 0);
    const double dy_y1 = coordinates(YH, p, p) - coordinates(YH, 0, p);

    for (int j = 0; j < p; ++j) {
      const double dx_dyh = mat.linear_scs_interp(0,j) * dx_x0 + mat.linear_scs_interp(1,j) * dx_y1;
      const double dy_dyh = mat.linear_scs_interp(0,j) * dy_x0 + mat.linear_scs_interp(1,j) * dy_y1;

      const double orth = dx_dyh * dx_dyh + dy_dyh * dy_dyh;
      for (int i = 0; i < p+1; ++i) {
        const double dx_dxh = mat.linear_nodal_interp(0,i) * dx_x1 + mat.linear_nodal_interp(1,i) * dx_y0;
        const double dy_dxh = mat.linear_nodal_interp(0,i) * dy_x1 + mat.linear_nodal_interp(1,i) * dy_y0;

        const double inv_detj = 1.0 / (dx_dyh * dy_dxh - dx_dxh * dy_dyh);
        metric(XH,XH,j,i) =  inv_detj * orth;
        metric(XH,YH,j,i) = -inv_detj * (dx_dxh * dx_dyh + dy_dxh * dy_dyh);
      }
    }

    for (int j = 0; j < p; ++j) {
      const double dx_dxh =  mat.linear_scs_interp(0,j) * dx_x1 + mat.linear_scs_interp(1,j) * dx_y0;
      const double dy_dxh =  mat.linear_scs_interp(0,j) * dy_x1 + mat.linear_scs_interp(1,j) * dy_y0;

      const double orth = dx_dxh * dx_dxh + dy_dxh * dy_dxh;
      for (int i = 0; i < p+1; ++i) {
        const double dx_dyh = mat.linear_nodal_interp(0,i) * dx_x0 + mat.linear_nodal_interp(1,i) * dx_y1;
        const double dy_dyh = mat.linear_nodal_interp(0,i) * dy_x0 + mat.linear_nodal_interp(1,i) * dy_y1;

        const double inv_detj = 1.0 / (dx_dyh * dy_dxh - dx_dxh * dy_dyh);
        metric(YH,XH,j,i) = -inv_detj * (dx_dxh * dx_dyh + dy_dxh * dy_dyh);
        metric(YH,YH,j,i) =  inv_detj * orth;
      }
    }
  }

  template <int p>
  void compute_diffusion_metric_linear(
    const CVFEMQuadOperators<p> ops,
    const nodal_vector_view<AlgTraitsQuad<p>> coordinates,
    const nodal_scalar_view<AlgTraitsQuad<p>> diffusivity,
    scs_tensor_view<AlgTraitsQuad<p>> metric)
  {
    const double dx_x0 = coordinates(XH, p, 0) - coordinates(XH, 0, 0);
    const double dx_x1 = coordinates(XH, 0, p) - coordinates(XH, 0, 0);
    const double dx_y0 = coordinates(XH, p, p) - coordinates(XH, p, 0);
    const double dx_y1 = coordinates(XH, p, p) - coordinates(XH, 0, p);

    const double dy_x0 = coordinates(YH, p, 0) - coordinates(YH, 0, 0);
    const double dy_x1 = coordinates(YH, 0, p) - coordinates(YH, 0, 0);
    const double dy_y0 = coordinates(YH, p, p) - coordinates(YH, p, 0);
    const double dy_y1 = coordinates(YH, p, p) - coordinates(YH, 0, p);

    const auto& mat = ops.mat_;
    nodal_scalar_view<AlgTraitsQuad<p>> diffIp{""};
    ops.scs_xhat_interp(diffusivity, diffIp);

    for (int j = 0; j < p; ++j) {
      const double dx_dyh = mat.linear_scs_interp(0,j) * dx_x0 + mat.linear_scs_interp(1,j) * dx_y1;
      const double dy_dyh = mat.linear_scs_interp(0,j) * dy_x0 + mat.linear_scs_interp(1,j) * dy_y1;

      const double orth = dx_dyh * dx_dyh + dy_dyh * dy_dyh;
      for (int i = 0; i < p + 1; ++i) {
        const double dx_dxh = mat.linear_nodal_interp(0,i) * dx_x1 + mat.linear_nodal_interp(1,i) * dx_y0;
        const double dy_dxh = mat.linear_nodal_interp(0,i) * dy_x1 + mat.linear_nodal_interp(1,i) * dy_y0;

        const double fac = diffIp(i,j) / (dx_dyh * dy_dxh - dx_dxh * dy_dyh);
        metric(XH,XH,j,i) =  fac * orth;
        metric(XH,YH,j,i) = -fac * (dx_dxh * dx_dyh + dy_dxh * dy_dyh);
      }
    }

    ops.scs_yhat_interp(diffusivity, diffIp);

    for (int j = 0; j < p; ++j) {
      const double dx_dxh = mat.linear_scs_interp(0,j) * dx_x1 + mat.linear_scs_interp(1,j) * dx_y0;
      const double dy_dxh = mat.linear_scs_interp(0,j) * dy_x1 + mat.linear_scs_interp(1,j) * dy_y0;

      const double orth = dx_dxh * dx_dxh + dy_dxh * dy_dxh;
      for (int i = 0; i < p + 1; ++i) {
        const double dx_dyh = mat.linear_nodal_interp(0,i) * dx_x0 + mat.linear_nodal_interp(1,i) * dx_y1;
        const double dy_dyh = mat.linear_nodal_interp(0,i) * dy_x0 + mat.linear_nodal_interp(1,i) * dy_y1;

        const double fac = diffIp(j,i) / (dx_dyh * dy_dxh - dx_dxh * dy_dyh);
        metric(YH,XH,j,i) = -fac * (dx_dxh * dx_dyh + dy_dxh * dy_dyh);
        metric(YH,YH,j,i) =  fac * orth;
      }
    }
  }

  // need a symmetric version tensor too, e.g. for the NSO

} // namespace HighOrderGeometryQuad
} // namespace naluUnit
} // namespace Sierra

#endif
