/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level NaluUnit      */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef HighOrderCoefficients_h
#define HighOrderCoefficients_h

#include <element_promotion/QuadratureRule.h>
#include <element_promotion/LagrangeBasis.h>
#include <element_promotion/CVFEMTypeDefs.h>
#include <stk_util/environment/ReportHandler.hpp>

#include <Teuchos_LAPACK.hpp>

namespace sierra {
namespace nalu {
namespace coefficients {
/* Computes 1D coefficient matrices (e.g. for the derivative) for CVFEM */

template <int p>
nodal_matrix_view<p>
nodal_integration_weights(
  const double* nodeLocs,
  const double* scsLocs)
{
  /*
  * Compute integration weights for CVFEM
  * This routine calculates integratation weights, \sum W_ij \phi(x_j) = \int_{x_{scs,i}}^{x_{scs,i+1}} \phi(x') dx',
  * at fixed integration point locations, with x_scs padded to include the end-points -1 and +1 by using a
  * "moment-matching" algorithm.  If "node locs" are the locations of the nodes
  */
  constexpr int nodes1D = p + 1;
  constexpr int nodesPerElement = (p + 1) * (p + 1);

  nodal_matrix_view<p> weightLHS("vandermonde matrix");
  for (int j = 0; j < nodes1D; ++j) {
    for (int i = 0; i < nodes1D; ++i) {
      weightLHS(j,i) = std::pow(nodeLocs[j], i);
    }
  }

  nodal_matrix_view<p> weights("nodal integration weighting for each scv");
  // each node has a separate RHS
  for (int i = 0; i < nodes1D; ++i) {
    weights(0,i) = (std::pow(scsLocs[0], i + 1) - std::pow(-1.0, i + 1)) / (i + 1.0);
  }

  for (int j = 1; j < nodes1D-1; ++j) {
    for (int i = 0; i < nodes1D; ++i) {
      weights(j,i) = (std::pow(scsLocs[j], i + 1) - std::pow(scsLocs[j-1], i + 1)) / (i + 1.0);
    }
  }

  for (int i = 0; i < nodes1D; ++i) {
    weights(p,i) = (std::pow(+1.0, i + 1) - std::pow(scsLocs[p-1], i + 1)) / (i + 1.0);
  }

  int info = 1;
  std::array<int, nodesPerElement> ipiv;
  Teuchos::LAPACK<int, double>().GESV(nodes1D, nodes1D,
    &weightLHS(0,0), nodes1D,
    ipiv.data(),
    &weights(0,0), nodes1D,
    &info
  );
  ThrowRequire(info == 0);

  // GESV overwrites the RHS with the solution
  return weights;
}
//--------------------------------------------------------------------------
template <int p >
scs_matrix_view<p>
scs_interpolation_weights(const double* nodeLocs, const double* scsLocs)
{
  constexpr int nodes1D = p+1;
  scs_matrix_view<p> scsInterp("subcontrol surface interpolation matrix");

  auto basis1D = Lagrange1D(nodeLocs, p);
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < nodes1D; ++i) {
      scsInterp(j,i) = basis1D.interpolation_weight(scsLocs[j], i);
    }
  }
  return scsInterp;
}
//--------------------------------------------------------------------------
template <int p>
scs_matrix_view<p>
scs_derivative_weights(
  const double* nodeLocs,
  const double* scsLocs)
{
  constexpr int nodes1D = p+1;
  scs_matrix_view<p> scsDeriv("subcontrol surface derivative matrix");

  auto basis1D = Lagrange1D(nodeLocs, p);
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < nodes1D; ++i) {
      scsDeriv(j,i) = basis1D.derivative_weight(scsLocs[j], i);
    }
  }
  return scsDeriv;
}
//--------------------------------------------------------------------------
template <int p >
nodal_matrix_view<p>
nodal_derivative_weights(const double* nodeLocs)
{
  constexpr int nodes1D = p+1;
  nodal_matrix_view<p> nodalDeriv("nodal derivative matrix");

  auto basis1D = Lagrange1D(nodeLocs, p);
  for (int j = 0; j < nodes1D; ++j) {
    for (int i = 0; i < nodes1D; ++i) {
      nodalDeriv(j,i) = basis1D.derivative_weight(nodeLocs[j],i);
    }
  }
  return nodalDeriv;
}
//--------------------------------------------------------------------------
template <int p>
linear_scs_matrix_view<p>
linear_scs_interpolation_weights(const double* scsLocs)
{
  linear_scs_matrix_view<p> linear_scs_interp("linscs");

  for (int j = 0; j < p; ++j) {
    linear_scs_interp(0,j) = 0.5*(1 - scsLocs[j]);
    linear_scs_interp(1,j) = 0.5*(1 + scsLocs[j]);
  }

  return linear_scs_interp;
}
//--------------------------------------------------------------------------
template <int p>
linear_nodal_matrix_view<p>
linear_nodal_interpolation_weights(const double* nodeLocs)
{
  linear_nodal_matrix_view<p> linear_nodal_interp("linnodal");

  for (int j = 0; j < p+1; ++j) {
    linear_nodal_interp(0,j) = 0.5*(1 - nodeLocs[j]);
    linear_nodal_interp(1,j) = 0.5*(1 + nodeLocs[j]);
  }
  return linear_nodal_interp;
}
//--------------------------------------------------------------------------
template <int p>
nodal_matrix_view<p>
nodal_integration_weights()
{
  auto nodeLocs = gauss_lobatto_legendre_rule(p+1).first;
  auto scsLocs  = gauss_legendre_rule(p).first;

  return nodal_integration_weights<p>(nodeLocs.data(), scsLocs.data());
}
//--------------------------------------------------------------------------
template <int p>
nodal_matrix_view<p>
nodal_derivative_weights()
{
  auto nodeLocs = gauss_lobatto_legendre_rule(p+1).first;
  return nodal_derivative_weights<p>(nodeLocs.data());
}
//--------------------------------------------------------------------------
template <int p>
scs_matrix_view<p>
scs_derivative_weights()
{
  auto nodeLocs = gauss_lobatto_legendre_rule(p+1).first;
  auto scsLocs  = gauss_legendre_rule(p).first;

  return scs_derivative_weights<p>(nodeLocs.data(), scsLocs.data());
}
//--------------------------------------------------------------------------
template <int p>
scs_matrix_view<p>
scs_interpolation_weights()
{
  auto nodeLocs = gauss_lobatto_legendre_rule(p+1).first;
  auto scsLocs  = gauss_legendre_rule(p).first;

  return scs_interpolation_weights<p>(nodeLocs.data(), scsLocs.data());
}
//--------------------------------------------------------------------------
template <int p>
linear_scs_matrix_view<p>
linear_scs_interpolation_weights()
{
  auto scsLocs  = gauss_legendre_rule(p).first;
  return linear_scs_interpolation_weights<p>(scsLocs.data());
}
//--------------------------------------------------------------------------
template <int p>
linear_nodal_matrix_view<p>
linear_nodal_interpolation_weights()
{
  auto nodeLocs = gauss_lobatto_legendre_rule(p+1).first;
  return linear_nodal_interpolation_weights<p>(nodeLocs.data());
}
//--------------------------------------------------------------------------
template <int p>
nodal_matrix_view<p> difference_matrix()
{
  nodal_matrix_view<p> scatt{"diff"};
  Kokkos::deep_copy(scatt, 0.0);

  scatt(0, 0)  = -1;
  scatt(p,p-1) = +1;

  for (int j = 1; j < p; ++j) {
    scatt(j,j+0) = -1;
    scatt(j,j-1) = +1;
  }
  return scatt;
}
//--------------------------------------------------------------------------
template <int p>
nodal_matrix_view<p> identity_matrix()
{
  nodal_matrix_view<p> id{ "" };
  Kokkos::deep_copy(id,0.0);
  for (int j = 0; j < p+1; ++j) {
    id(j,j) = 1.0;
  }

  return id;
}

} // namespace CoefficientMatrices
} // namespace naluUnit
} // namespace Sierra

#endif
