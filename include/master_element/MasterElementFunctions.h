/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MasterElementGeneric_h
#define MasterElementGeneric_h

#include <AlgTraits.h>

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <array>
#include <type_traits>

#include <master_element/MasterElement.h>
#include <SimdInterface.h>
#include <Kokkos_Core.hpp>
#include <element_promotion/MasterElementUtils.h>
#include <stk_util/environment/ReportHandler.hpp>

namespace sierra {
namespace nalu {

  template <typename AlgTraits, typename GradViewType, typename CoordViewType, typename OutputViewType>
  void generic_grad_op_3d(GradViewType referenceGradWeights, CoordViewType coords, OutputViewType weights)
  {
    /**
    * Given the reference gradient weights evaluated at the integration points and the cooordinates,
    * this method computes
    *     \nabla_\mathbf{x}|_{\alpha}  as \left( J^{-T} \nabla_{\mathbf{x}^\hat} \right)|_\alpha,
    *
    *  This operation can be specialized for efficiency on hex topologies (as tensor-contractions)
    *  or on tets (since the gradient is independent of \alpha).  But this can work as a fallback
    */

    using ftype = typename CoordViewType::value_type;
    static_assert(std::is_same<ftype, typename GradViewType::value_type>::value,  "Incompatiable value type for views");
    static_assert(std::is_same<ftype, typename OutputViewType::value_type>::value,  "Incompatiable value type for views");
    static_assert(GradViewType::Rank == 3, "grad view assumed to be 3D");
    static_assert(CoordViewType::Rank == 2, "Coordinate view assumed to be 2D");
    static_assert(OutputViewType::Rank == 3, "Weight view assumed to be 3D");
    static_assert(AlgTraits::nDim_ == 3, "3D method");

    ThrowAssert(AlgTraits::nodesPerElement_ == referenceGradWeights.extent(1));
    ThrowAssert(AlgTraits::nDim_ == referenceGradWeights.extent(2));
    ThrowAssert(weights.extent(0) == referenceGradWeights.extent(0));
    ThrowAssert(weights.extent(1) == referenceGradWeights.extent(1));
    ThrowAssert(weights.extent(2) == referenceGradWeights.extent(2));

    for (unsigned ip = 0; ip < referenceGradWeights.extent(0); ++ip) {
      ftype jact[3][3] = {
          {0.0, 0.0, 0.0},
          {0.0, 0.0, 0.0},
          {0.0, 0.0, 0.0}
      };

      ftype refGrad[AlgTraits::nodesPerElement_][3];
      for (int n = 0; n < AlgTraits::nodesPerElement_; ++n) {
        refGrad[n][0] = referenceGradWeights(ip, n, 0);
        refGrad[n][1] = referenceGradWeights(ip, n, 1);
        refGrad[n][2] = referenceGradWeights(ip, n, 2);

        jact[0][0] += refGrad[n][0] * coords(n, 0);
        jact[0][1] += refGrad[n][1] * coords(n, 0);
        jact[0][2] += refGrad[n][2] * coords(n, 0);

        jact[1][0] += refGrad[n][0] * coords(n, 1);
        jact[1][1] += refGrad[n][1] * coords(n, 1);
        jact[1][2] += refGrad[n][2] * coords(n, 1);

        jact[2][0] += refGrad[n][0] * coords(n, 2);
        jact[2][1] += refGrad[n][1] * coords(n, 2);
        jact[2][2] += refGrad[n][2] * coords(n, 2);
      }

      ftype adjJac[3][3];
      adjJac[0][0] = jact[1][1] * jact[2][2] - jact[2][1] * jact[1][2];
      adjJac[0][1] = jact[1][2] * jact[2][0] - jact[2][2] * jact[1][0];
      adjJac[0][2] = jact[1][0] * jact[2][1] - jact[2][0] * jact[1][1];

      adjJac[1][0] = jact[0][2] * jact[2][1] - jact[2][2] * jact[0][1];
      adjJac[1][1] = jact[0][0] * jact[2][2] - jact[2][0] * jact[0][2];
      adjJac[1][2] = jact[0][1] * jact[2][0] - jact[2][1] * jact[0][0];

      adjJac[2][0] = jact[0][1] * jact[1][2] - jact[1][1] * jact[0][2];
      adjJac[2][1] = jact[0][2] * jact[1][0] - jact[1][2] * jact[0][0];
      adjJac[2][2] = jact[0][0] * jact[1][1] - jact[1][0] * jact[0][1];

     ThrowAssertMsg(
       stk::simd::are_any(
         jact[0][0] * adjJac[0][0] + jact[1][0] * adjJac[1][0] + jact[2][0] * adjJac[2][0]
         > tiny_positive_value()
       ),
       "Problem with Jacobian determinant"
     );

     const ftype inv_detj = ftype(1.0) /
         (jact[0][0] * adjJac[0][0] + jact[1][0] * adjJac[1][0] + jact[2][0] * adjJac[2][0]);

     for (int n = 0; n < AlgTraits::nodesPerElement_; ++n) {
       weights(ip, n, 0) = inv_detj *
           (adjJac[0][0] * refGrad[n][0] + adjJac[0][1] * refGrad[n][1] + adjJac[0][2] * refGrad[n][2]);

       weights(ip, n, 1) = inv_detj *
           (adjJac[1][0] * refGrad[n][0] + adjJac[1][1] * refGrad[n][1] + adjJac[1][2] * refGrad[n][2]);

       weights(ip, n, 2) = inv_detj *
           (adjJac[2][0] * refGrad[n][0] + adjJac[2][1] * refGrad[n][1] + adjJac[2][2] * refGrad[n][2]);
     }
    }
  }

  template<typename AlgTraits, typename GradViewType,  typename CoordViewType,
           typename DetjType, typename OutputViewType>
  void generic_grad_op_3d(
    GradViewType referenceGradWeights,
    CoordViewType coords,
    OutputViewType weights,
    DetjType detj)
  {
    /**
    * Given the reference gradient weights evaluated at the integration points and the cooordinates,
    * this method computes
    *     \nabla_\mathbf{x}|_{\alpha}  as \left( J^{-T} \nabla_{\mathbf{x}^\hat} \right)|_\alpha,
    *
    *  This operation can be specialized for efficiency on hex topologies (as tensor-contractions)
    *  or on tets (since the gradient is independent of \alpha).  But this can work as a fallback
    */

    using ftype = typename CoordViewType::value_type;
    static_assert(std::is_same<ftype, typename GradViewType::value_type>::value,  "Incompatiable value type for views");
    static_assert(std::is_same<ftype, typename OutputViewType::value_type>::value,  "Incompatiable value type for views");
    static_assert(GradViewType::Rank == 3, "grad view assumed to be 3D");
    static_assert(CoordViewType::Rank == 2, "Coordinate view assumed to be 2D");
    static_assert(OutputViewType::Rank == 3, "Weight view assumed to be 3D");
    static_assert(AlgTraits::nDim_ == 3, "3D method");

    ThrowAssert(AlgTraits::nodesPerElement_ == referenceGradWeights.extent(1));
    ThrowAssert(AlgTraits::nDim_ == referenceGradWeights.extent(2));
    ThrowAssert(weights.extent(0) == referenceGradWeights.extent(0));
    ThrowAssert(weights.extent(1) == referenceGradWeights.extent(1));
    ThrowAssert(weights.extent(2) == referenceGradWeights.extent(2));
    ThrowAssert(detj.extent(0) == weights.extent(0));


    for (unsigned ip = 0; ip < referenceGradWeights.extent(0); ++ip) {

      ftype jac[3][3] = {
          {0.0, 0.0, 0.0},
          {0.0, 0.0, 0.0},
          {0.0, 0.0, 0.0}
      };

     for (int n = 0; n < AlgTraits::nodesPerElement_; ++n) {
        for (int outer_dim = 0; outer_dim < 3; ++outer_dim) {
          for (int inner_dim  = 0; inner_dim < 3; ++inner_dim) {
            jac[outer_dim][inner_dim] += referenceGradWeights(ip, n, inner_dim) * coords(n, outer_dim);
          }
        }
      }

     ftype adj_jac[3][3];
     adj_jac[0][0] = jac[1][1] * jac[2][2] - jac[2][1] * jac[1][2];
     adj_jac[0][1] = jac[0][2] * jac[2][1] - jac[2][2] * jac[0][1];
     adj_jac[0][2] = jac[0][1] * jac[1][2] - jac[1][1] * jac[0][2];

     adj_jac[1][0] = jac[1][2] * jac[2][0] - jac[2][2] * jac[1][0];
     adj_jac[1][1] = jac[0][0] * jac[2][2] - jac[2][0] * jac[0][2];
     adj_jac[1][2] = jac[0][2] * jac[1][0] - jac[1][2] * jac[0][0];

     adj_jac[2][0] = jac[1][0] * jac[2][1] - jac[2][0] * jac[1][1];
     adj_jac[2][1] = jac[0][1] * jac[2][0] - jac[2][1] * jac[0][0];
     adj_jac[2][2] = jac[0][0] * jac[1][1] - jac[1][0] * jac[0][1];

     detj(ip) = jac[0][0] * adj_jac[0][0] + jac[1][0] * adj_jac[0][1] + jac[2][0] * adj_jac[0][2];

     ThrowAssertMsg(stk::simd::are_any(detj(ip) > +tiny_positive_value()),"Problem with determinant");

     const ftype inv_detj = ftype(1.0) / detj(ip);

     for (int n = 0; n < AlgTraits::nodesPerElement_; ++n) {
       const ftype dnds[3] =
           { referenceGradWeights(ip, n, 0), referenceGradWeights(ip, n, 1), referenceGradWeights(ip, n, 2) };

       weights(ip, n, 0) = inv_detj *
           (adj_jac[0][0] * dnds[0] + adj_jac[1][0] * dnds[1] + adj_jac[2][0] * dnds[2]);

       weights(ip, n, 1) = inv_detj *
           (adj_jac[0][1] * dnds[0] + adj_jac[1][1] * dnds[1] + adj_jac[2][1] * dnds[2]);

       weights(ip, n, 2) = inv_detj *
           (adj_jac[0][2] * dnds[0] + adj_jac[1][2] * dnds[1] + adj_jac[2][2] * dnds[2]);
     }
    }
  }

  template <typename AlgTraits, typename GradViewType, typename CoordViewType, typename OutputViewType>
  void generic_gij_3d(GradViewType referenceGradWeights, CoordViewType coords, OutputViewType gup, OutputViewType glo)
  {
    using ftype = typename CoordViewType::value_type;
    static_assert(std::is_same<ftype, typename GradViewType::value_type>::value,
      "Incompatiable value type for views");
    static_assert(std::is_same<ftype, typename OutputViewType::value_type>::value,
      "Incompatiable value type for views");
    static_assert(GradViewType::Rank == 3, "grad view assumed to be 3D");
    static_assert(CoordViewType::Rank == 2, "Coordinate view assumed to be 2D");
    static_assert(OutputViewType::Rank == 3, "gij view assumed to be 3D");
    static_assert(AlgTraits::nDim_ == 3, "3D method");

    for (unsigned ip = 0; ip < referenceGradWeights.extent(0); ++ip) {

      ftype jac[3][3] = {
          {0.0, 0.0, 0.0},
          {0.0, 0.0, 0.0},
          {0.0, 0.0, 0.0}
      };

     for (int n = 0; n < AlgTraits::nodesPerElement_; ++n) {
        for (int outer_dim = 0; outer_dim < 3; ++outer_dim) {
          for (int inner_dim  = 0; inner_dim < 3; ++inner_dim) {
            jac[outer_dim][inner_dim] += referenceGradWeights(ip, n, inner_dim) * coords(n, outer_dim);
          }
        }
      }

      gup(ip, 0, 0) = jac[0][0] * jac[0][0] + jac[0][1] * jac[0][1] + jac[0][2] * jac[0][2];
      gup(ip, 0, 1) = jac[0][0] * jac[1][0] + jac[0][1] * jac[1][1] + jac[0][2] * jac[1][2];
      gup(ip, 0, 2) = jac[0][0] * jac[2][0] + jac[0][1] * jac[2][1] + jac[0][2] * jac[2][2];

      gup(ip, 1, 0) = gup(ip, 0, 1);
      gup(ip, 1, 1) = jac[1][0] * jac[1][0] + jac[1][1] * jac[1][1] + jac[1][2] * jac[1][2];
      gup(ip, 1, 2) = jac[1][0] * jac[2][0] + jac[1][1] * jac[2][1] + jac[1][2] * jac[2][2];

      gup(ip, 2, 0) = gup(ip, 0, 2);
      gup(ip, 2, 1) = gup(ip, 1, 2);
      gup(ip, 2, 2) = jac[2][0] * jac[2][0] + jac[2][1] * jac[2][1] + jac[2][2] * jac[2][2];

      // the covariant is the inverse of the contravariant by definition
      // gUpper is symmetric
      const ftype inv_detj = ftype(1.0) / (
            gup(ip, 0, 0) * ( gup(ip, 1, 1) * gup(ip, 2, 2) - gup(ip, 1, 2) * gup(ip, 1, 2) )
          - gup(ip, 0, 1) * ( gup(ip, 0, 1) * gup(ip, 2, 2) - gup(ip, 1, 2) * gup(ip, 0, 2) )
          + gup(ip, 0, 2) * ( gup(ip, 0, 1) * gup(ip, 1, 2) - gup(ip, 1, 1) * gup(ip, 0, 2) )
      );

      glo(ip, 0, 0) = inv_detj * (gup(ip, 1, 1) * gup(ip, 2, 2) - gup(ip, 1, 2) * gup(ip, 1, 2));
      glo(ip, 0, 1) = inv_detj * (gup(ip, 0, 2) * gup(ip, 1, 2) - gup(ip, 0, 1) * gup(ip, 2, 2));
      glo(ip, 0, 2) = inv_detj * (gup(ip, 0, 1) * gup(ip, 1, 2) - gup(ip, 0, 2) * gup(ip, 1, 1));

      glo(ip, 1, 0) = glo(ip, 0, 1);
      glo(ip, 1, 1) = inv_detj * (gup(ip, 0, 0) * gup(ip, 2, 2) - gup(ip, 0, 2) * gup(ip, 0, 2));
      glo(ip, 1, 2) = inv_detj * (gup(ip, 0, 2) * gup(ip, 0, 1) - gup(ip, 0, 0) * gup(ip, 1, 2));


      glo(ip, 2, 0) = glo(ip, 0, 2);
      glo(ip, 2, 1) = glo(ip, 1, 2);
      glo(ip, 2, 2) = inv_detj * (gup(ip, 0, 0) * gup(ip, 1, 1) - gup(ip, 0, 1) * gup(ip, 0, 1));

    }
  }

  template <typename AlgTraits, typename GradViewType, typename CoordViewType, typename OutputViewType>
  void generic_determinant_3d(GradViewType referenceGradWeights, CoordViewType coords, OutputViewType detj)
  {
    using ftype = typename CoordViewType::value_type;
    static_assert(std::is_same<ftype, typename GradViewType::value_type>::value,  "Incompatiable value type for views");
    static_assert(std::is_same<ftype, typename OutputViewType::value_type>::value,  "Incompatiable value type for views");
    static_assert(GradViewType::Rank == 3, "grad view assumed to be 3D");
    static_assert(CoordViewType::Rank == 2, "Coordinate view assumed to be 2D");
    static_assert(OutputViewType::Rank == 1, "Weight view assumed to be 1D");
    static_assert(AlgTraits::nDim_ == 3, "3D method");

    for (int ip = 0; ip < AlgTraits::numScsIp_; ++ip) {
      ftype jac[3][3] = {
          {0.0, 0.0, 0.0},
          {0.0, 0.0, 0.0},
          {0.0, 0.0, 0.0}
      };
      for (int n = 0; n < AlgTraits::nodesPerElement_; ++n) {
        for (int outer_dim = 0; outer_dim < AlgTraits::nDim_; ++outer_dim) {
          for (int inner_dim  = 0; inner_dim < AlgTraits::nDim_; ++inner_dim) {
            jac[outer_dim][inner_dim] += referenceGradWeights(ip, n, outer_dim) * coords(n, inner_dim);
          }
        }
      }

      detj(ip) = (
            jac[0][0] * ( jac[1][1] * jac[2][2] - jac[1][2] * jac[1][2] )
          - jac[0][1] * ( jac[0][1] * jac[2][2] - jac[1][2] * jac[0][2] )
          + jac[0][2] * ( jac[0][1] * jac[1][2] - jac[1][1] * jac[0][2] )
      );
    }
  }

  template <typename ViewType>
  void quad_area_by_triangulation(
    int ics,
    const typename ViewType::value_type areacoords[4][3],
    const ViewType& area)
  {
    /**
     * Form up the area vector consistently with the triangulation used
     * in the Grandy algorithm, on each subcontrol volume hex
     *
     * "Efficient computation of volume of
     * Hexahedral Cells", Jeffrey Grandy, LLNL, UCRL-ID-128886,
     *  October 30, 1997.
     */
    using ftype = typename ViewType::value_type;

    constexpr int triangularFacetTable[4][3] = {
        {4, 0, 1},
        {4, 1, 2},
        {4, 2, 3},
        {4, 3, 0}
    };

    ftype xmid[3];
    ftype r2[3];
    for(int k=0; k < 3; ++k) {
      xmid[k] = 0.25 * (areacoords[0][k] + areacoords[1][k] + areacoords[2][k] + areacoords[3][k]);
      area(ics,k) = ftype(0.0);
      r2[k] = areacoords[0][k] - xmid[k];
    }
    constexpr int ntriangles = 4;
    ftype r1[3];
    for(int itriangle = 0; itriangle < ntriangles; ++itriangle) {
      const int q = triangularFacetTable[itriangle][2];
      for (int k = 0; k < 3; ++k) {
        r1[k] = r2[k];
        r2[k] = areacoords[q][k] - xmid[k];
      }
      area(ics,0) += r1[1]*r2[2] - r2[1]*r1[2];
      area(ics,1) += r1[2]*r2[0] - r2[2]*r1[0];
      area(ics,2) += r1[0]*r2[1] - r2[0]*r1[1];
    }

    for (int k=0; k < 3; ++k) {
      area(ics,k) *= 0.5;
    }
  }

  template <typename RealType>
  RealType hex_volume_grandy(RealType scvcoords[8][3])
  {
    /**
     * The Grandy algorithm for computing the volume of a multilinear box
     *
     * "Efficient computation of volume of
     * Hexahedral Cells", Jeffrey Grandy, LLNL, UCRL-ID-128886,
     *  October 30, 1997.
     */
    constexpr int nTri = 24;
    constexpr int dim = 3;

    constexpr int nNodes = 8;
    constexpr int nFaces = 6;
    constexpr int npv = nNodes + nFaces;

    RealType coordv[npv][dim];

    // copy coordinates
    for (int n = 0; n < nNodes; ++n) {
      for (int d = 0; d < dim; ++d) {
        coordv[n][d] = scvcoords[n][d];
      }
    }

    constexpr int nodesPerFace = 4;
    constexpr int face_nodes[nFaces][nodesPerFace] = {
        { 0, 3, 2, 1 }, { 4, 5, 6, 7 },
        { 0, 1, 5, 4 }, { 2, 3, 7, 6 },
        { 1, 2, 6, 5 }, { 0, 4, 3, 7 }
    };

    // append face midpoint coordinates
    for (int k = 0; k < nFaces; ++k) {
      const int coordIndex = k + nNodes;
      for (int d = 0; d < dim; ++d) {
        coordv[coordIndex][d] = 0.25 * (
             coordv[face_nodes[k][0]][d]
           + coordv[face_nodes[k][1]][d]
           + coordv[face_nodes[k][2]][d]
           + coordv[face_nodes[k][3]][d]
        );
      }
    }

    constexpr int triangular_facets[nTri][3] = {
        { 0,  8,  1}, { 8,  2,  1}, { 3,  2,  8},
        { 3,  8,  0}, { 6,  9,  5}, { 7,  9,  6},
        { 4,  9,  7}, { 4,  5,  9}, {10,  0,  1},
        { 5, 10,  1}, { 4, 10,  5}, { 4,  0, 10},
        { 7,  6, 11}, { 6,  2, 11}, { 2,  3, 11},
        { 3,  7, 11}, { 6, 12,  2}, { 5, 12,  6},
        { 5,  1, 12}, { 1,  2, 12}, { 0,  4, 13},
        { 4,  7, 13}, { 7,  3, 13}, { 3,  0, 13}
    };

    RealType volume = 0.0;

    for (int k = 0; k < nTri; ++k) {
      const int p = triangular_facets[k][0];
      const int q = triangular_facets[k][1];
      const int r = triangular_facets[k][2];

      RealType triFaceMid[3];
      for (int d = 0; d < 3; ++d) {
        triFaceMid[d] = coordv[p][d] + coordv[q][d] + coordv[r][d];
      }

      enum {XC = 0, YC = 1, ZC = 2};
      RealType dxv[3];

      dxv[0] = ( coordv[q][YC] - coordv[p][YC] ) * ( coordv[r][ZC] - coordv[p][ZC] )
             - ( coordv[r][YC] - coordv[p][YC] ) * ( coordv[q][ZC] - coordv[p][ZC] );

      dxv[1] = ( coordv[r][XC] - coordv[p][XC] ) * ( coordv[q][ZC] - coordv[p][ZC] )
             - ( coordv[q][XC] - coordv[p][XC] ) * ( coordv[r][ZC] - coordv[p][ZC] );

      dxv[2] = ( coordv[q][XC] - coordv[p][XC] ) * ( coordv[r][YC] - coordv[p][YC] )
             - ( coordv[r][XC] - coordv[p][XC] ) * ( coordv[q][YC] - coordv[p][YC] );

      volume += triFaceMid[0] * dxv[0] + triFaceMid[1] * dxv[1] + triFaceMid[2] * dxv[2];
    }
    volume /= RealType(18.0);
    return volume;
  }

  template <typename CoordViewType>
  void subdivide_hex_8(CoordViewType coords, typename CoordViewType::value_type coordv[27][3])
  {
    /**
     * Subdivide the coordinates of a hex8 element into 8 hexs along edge, face, and volume midpoints
     */
    constexpr int dim = 3;
    constexpr int numBaseNodes = 8;

    for (int n = 0; n < numBaseNodes; ++n) {
      for (int d = 0; d < dim; ++d) {
        coordv[n][d] = coords(n,d);
      }
    }

    // Face-by-face ordering for the subdivided hex.  This is different than what is done
    // for a multilinear Hex27 element, which has equivalent nodal locations.
    for (int d = 0; d < 3; ++d) {
      // Face 1
      coordv[8][d] = 0.5 * (coords(0,d) + coords(1,d));
      coordv[9][d] = 0.5 * (coords(1,d) + coords(2,d));
      coordv[10][d] = 0.5 * (coords(2,d) + coords(3,d));
      coordv[11][d] = 0.5 * (coords(3,d) + coords(0,d));
      coordv[12][d] = 0.25 * (coords(0,d) + coords(1,d) + coords(2,d) + coords(3,d));

      // Face 2
      coordv[13][d] = 0.5 * (coords(4,d) + coords(5,d));
      coordv[14][d] = 0.5 * (coords(5,d) + coords(6,d));
      coordv[15][d] = 0.5 * (coords(6,d) + coords(7,d));
      coordv[16][d] = 0.5 * (coords(7,d) + coords(4,d));
      coordv[17][d] = 0.25 * (coords(4,d) + coords(5,d) + coords(6,d) + coords(7,d));

      // Face 3
      coordv[18][d] = 0.5 * (coords(1,d) + coords(5,d));
      coordv[19][d] = 0.5 * (coords(0,d) + coords(4,d));
      coordv[20][d] = 0.25 * (coords(0,d) + coords(1,d) + coords(4,d) + coords(5,d));

      // Face 4
      coordv[21][d] = 0.5 * (coords(3,d) + coords(7,d));
      coordv[22][d] = 0.5 * (coords(2,d) + coords(6,d));
      coordv[23][d] = 0.25 * (coords(2,d) + coords(3,d) + coords(6,d) + coords(7,d));

      // Face 5
      coordv[24][d] = 0.25 * (coords(1,d) + coords(2,d) + coords(5,d) + coords(6,d));

      // Face 7
      coordv[25][d] = 0.25 * (coords(0,d) + coords(3,d) + coords(4,d) + coords(7,d));

      // Volume centroid
      coordv[26][d] = 0.;
      for (int nd = 0; nd < 8; ++nd) {
        coordv[26][d] += coords(nd,d);
      }
      coordv[26][d] *= 0.125;
    }
  }

} // namespace nalu
} // namespace Sierra

#endif
