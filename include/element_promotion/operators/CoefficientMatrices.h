/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level NaluUnit      */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef CoefficientMatrices_h
#define CoefficientMatrices_h

#include <element_promotion/operators/HighOrderCoefficients.h>

namespace sierra {
namespace nalu{

template <int p>
struct CoefficientMatrices
{
  constexpr static int poly_order = p;

  CoefficientMatrices(const double* nodeLocs, const double* scsLocs)
  : scsDeriv(coefficients::scs_derivative_weights<poly_order>(nodeLocs, scsLocs)),
    scsInterp(coefficients::scs_interpolation_weights<poly_order>(nodeLocs, scsLocs)),
    nodalWeights(coefficients::nodal_integration_weights<poly_order>(nodeLocs, scsLocs)),
    nodalDeriv(coefficients::nodal_derivative_weights<poly_order>(nodeLocs)),
    linear_nodal_interp(coefficients::linear_nodal_interpolation_weights<poly_order>(nodeLocs)),
    linear_scs_interp(coefficients::linear_scs_interpolation_weights<poly_order>(scsLocs))
  {};

  CoefficientMatrices()
  : scsDeriv(coefficients::scs_derivative_weights<poly_order>()),
    scsInterp(coefficients::scs_interpolation_weights<poly_order>()),
    nodalWeights(coefficients::nodal_integration_weights<poly_order>()),
    nodalDeriv(coefficients::nodal_derivative_weights<poly_order>()),
    difference(coefficients::difference_matrix<poly_order>()),
    linear_nodal_interp(coefficients::linear_nodal_interpolation_weights<poly_order>()),
    linear_scs_interp(coefficients::linear_scs_interpolation_weights<poly_order>())
  {};

  const scs_matrix_view<p> scsDeriv;
  const scs_matrix_view<p> scsInterp;
  const nodal_matrix_view<p> nodalWeights;
  const nodal_matrix_view<p> nodalDeriv;
  const nodal_matrix_view<p> difference;
  const linear_nodal_matrix_view<p> linear_nodal_interp;
  const linear_scs_matrix_view<p> linear_scs_interp;
};

} // namespace naluUnit
} // namespace Sierra

#endif
