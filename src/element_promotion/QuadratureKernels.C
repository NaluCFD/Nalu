/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level nalu      */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include <element_promotion/QuadratureKernels.h>

#include <element_promotion/QuadratureRule.h>

#include <Teuchos_BLAS.hpp>

#include <tuple>
#include <element_promotion/ElementDescription.h>


namespace sierra {
namespace nalu {

//==========================================================================
// Class Definition
//==========================================================================
// SGLQuadratureOps - Wraps some of the BLAS tools / data used to integrate
// over surfaces/volumes in 2D/3D using the "Segmented Gauss Lobatto"
// quadrature
//===========================================================================

SGLQuadratureOps::SGLQuadratureOps(const ElementDescription& elem)
: blas_(Teuchos::BLAS<int,double>())
{

  nodes1D_ = elem.nodes1D;
  nodesPerElement_ = elem.nodesPerElement;
  numSurfaces_ = elem.dimension * elem.polyOrder;
  nodesPerFace_ = elem.nodesPerSide;

  auto scsEndLoc = pad_end_points(gauss_legendre_rule(elem.nodes1D-1).first);
  weightTensor_ = SGL_quadrature_rule(elem.nodes1D, scsEndLoc).second;
  p_weightTensor_ = weightTensor_.data();

  work2D_.resize(elem.nodes1D*elem.nodes1D);
  p_work2D_ = work2D_.data();

  // For 3D volume integrals, the weights are saved to a large N^3 x N^3 matrix
  // where N is the number of nodes in 1D
  if (elem.dimension == 3) {
    size3D_ = numSurfaces_*nodes1D_*nodes1D_;

    // forms a (re)mapped weight matrix from the weights in the summation
    // The mapping means that the result of the integration is ordered like
    // the nodes in the element itself (whatever that may be)
    weightMatrix_.resize(elem.nodesPerElement*elem.nodesPerElement);
    for (int q = 0; q < elem.nodesPerElement; ++q) {
      const auto& lmn = elem.inverseNodeMap[q];
      for (int p = 0; p < elem.nodesPerElement; ++p) {
        const auto& ijk = elem.inverseNodeMap[p];
        weightMatrix_[p + q * elem.nodesPerElement] = weightTensor_[ijk[0] * elem.nodes1D + lmn[0]]
                                                    * weightTensor_[ijk[1] * elem.nodes1D + lmn[1]]
                                                    * weightTensor_[ijk[2] * elem.nodes1D + lmn[2]];
      }
    }
    p_weightMatrix_ = weightMatrix_.data();
  }
}
//--------------------------------------------------------------------------
void
SGLQuadratureOps::volume_2D(
  const double* nodalValuesTensor,
  double* result)
{
  /* Computes
   * \bar{\phi}_{lm} = \sum_m^n \sum_l^n \sum_j^n \sum_i W_{mj} W_{li} phi_{ij}
   */

  blas_.GEMM(
    Teuchos::NO_TRANS,
    Teuchos::NO_TRANS,
    nodes1D_, nodes1D_, nodes1D_,
    1.0,
    nodalValuesTensor, nodes1D_,
    p_weightTensor_, nodes1D_,
    0.0,
    p_work2D_, nodes1D_
  );

  blas_.GEMM(
    Teuchos::TRANS,
    Teuchos::NO_TRANS,
    nodes1D_, nodes1D_, nodes1D_,
    1.0,
    p_weightTensor_, nodes1D_,
    p_work2D_, nodes1D_,
    0.0,
    result, nodes1D_
  );
}
//--------------------------------------------------------------------------
void
SGLQuadratureOps::volume_3D(
  const double*  nodalValues,
  double* result)
{
  /* Computes
   * \bar{\phi}_{lmn} = \sum_n^N \sum_m^N \sum_l^N
   *                    \sum_k^N \sum_j^N \sum_i^N W_{nk} W_{mj} W_{li} phi_{ijk}
   *
   * Or, more directly,
   * \bar{phi}_q = \sum^{N^3}_{qp} \tilde{W}_{qp} \phi_p
   *
   * TODO(rcknaus):
   * This can be computed with fewer flops by using a sequence of smaller
   * matrix-matrix multiplications.  For now, just do a big matvec.
   */

  blas_.GEMV(
    Teuchos::NO_TRANS,
    nodesPerElement_, nodesPerElement_, 1.0, p_weightMatrix_,
    nodesPerElement_, nodalValues, 1, 0.0,
    result, 1
  );
}
//--------------------------------------------------------------------------
void SGLQuadratureOps::surface_2D(
  const double*  integrand,
  double* result,
  int line_offset)
{
  /* Computes a line integral for each of the  N cv segments along a subcontrol line
   * $\bar{\phi}_l = \sum_l^N \sum_\alpha W_{l \alpha} phi_{\alpha}$
   *
   * where greek indices indicate integration point indices
   * and roman indices are nodal indices
   */

  blas_.GEMV(
    Teuchos::TRANS,
    nodes1D_, nodes1D_,
    1.0,
    p_weightTensor_, nodes1D_,
    integrand + line_offset, 1,
    0.0,
    result + line_offset, 1
  );
}
//--------------------------------------------------------------------------
void
SGLQuadratureOps::surface_3D(
  const double* integrand,
  double* result,
  int face_offset)
{
  /* Computes a surface integral for each of the N*N cv surfaces along a subcontrol plane
   * \bar{\phi}_lm = \sum_m^n \sum_l^n
   *                 \sum_\beta^n \sum_\alpha W_{m \beta} W_{l \alpha} phi_{\alpha \beta}
   *
   * where greek indices indicate integration point indices
   * and roman indices are nodal indices
   */

  blas_.GEMM(
    Teuchos::NO_TRANS,
    Teuchos::NO_TRANS,
    nodes1D_, nodes1D_, nodes1D_,
    1.0,
    integrand + face_offset, nodes1D_,
    p_weightTensor_, nodes1D_,
    0.0,
    p_work2D_, nodes1D_
  );

  blas_.GEMM(
    Teuchos::TRANS,
    Teuchos::NO_TRANS,
    nodes1D_, nodes1D_, nodes1D_,
    1.0,
    p_weightTensor_, nodes1D_,
    p_work2D_, nodes1D_,
    0.0,
    result + face_offset, nodes1D_
  );
}
//--------------------------------------------------------------------------
void SGLQuadratureOps::surfaces_2D(const double* integrand, double* result)
{
  // Computes all surfaces integrals for a 2D Quad element

  blas_.GEMM(
    Teuchos::TRANS,
    Teuchos::NO_TRANS,
    nodes1D_, numSurfaces_, nodes1D_,
    1.0,
    p_weightTensor_, nodes1D_,
    integrand, nodes1D_,
    0.0,
    result, nodes1D_
  );
}
//--------------------------------------------------------------------------
void SGLQuadratureOps::surfaces_3D(const double* integrand, double* result)
{
  // Computes all surface integrals for a 3D Hex element
  // TODO(rcknaus): optimization

  for (int face_offset = 0; face_offset < size3D_; face_offset += nodesPerFace_) {
    blas_.GEMM(
      Teuchos::NO_TRANS,
      Teuchos::NO_TRANS,
      nodes1D_, nodes1D_, nodes1D_,
      1.0,
      integrand + face_offset, nodes1D_,
      p_weightTensor_, nodes1D_,
      0.0,
      p_work2D_, nodes1D_
    );

    blas_.GEMM(
      Teuchos::TRANS,
      Teuchos::NO_TRANS,
      nodes1D_, nodes1D_, nodes1D_,
      1.0,
      p_weightTensor_, nodes1D_,
      p_work2D_, nodes1D_,
      0.0,
      result + face_offset, nodes1D_
    );
  }
}

} // namespace nalu
} // namespace Sierra
