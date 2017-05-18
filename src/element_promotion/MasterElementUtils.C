/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporatlion.                                   */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#include <element_promotion/MasterElementUtils.h>
#include <element_promotion/LagrangeBasis.h>

#include <NaluEnv.h>
#include <FORTRAN_Proto.h>

#include <stk_util/environment/ReportHandler.hpp>

#include <array>
#include <limits>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <element_promotion/ElementDescription.h>

namespace sierra{
namespace nalu{

  double ddot(const double* u, const double* v, int n)
  {
    double val = 0.0;
    for (int i = 0; i < n; ++i) {
      val += u[i] * v[i];
    }
    return val;
  }

  double determinant33(const double* mat)
  {
    enum { XX = 0, XY = 1, XZ = 2, YX = 3, YY = 4, YZ = 5, ZX = 6, ZY = 7, ZZ = 8 };
    return (mat[XX] * (mat[YY] * mat[ZZ] - mat[YZ] * mat[ZY]) - mat[XY] * (mat[YX] * mat[ZZ] - mat[YZ] * mat[ZX]) +
        mat[XZ] * (mat[YX] * mat[ZY] - mat[YY] * mat[ZX]));
  }

  void action_of_inverse33(
    const double* POINTER_RESTRICT A,
    const double* POINTER_RESTRICT x,
    double*  POINTER_RESTRICT b)
  {
    ThrowAssert(determinant33(A) > tiny_positive_value());
    enum { XX = 0, XY = 1, XZ = 2, YX = 3, YY = 4, YZ = 5, ZX = 6, ZY = 7, ZZ = 8 };
    enum { X_RANK1 = 0, Y_RANK1 = 1, Z_RANK1 = 2 };

    const double inv_detA = 1.0 / determinant33(A);
    b[X_RANK1] = ((A[YY] * A[ZZ] - A[YZ] * A[ZY]) * x[X_RANK1] + (A[XZ] * A[ZY] - A[XY] * A[ZZ]) * x[Y_RANK1] +
        (A[XY] * A[YZ] - A[XZ] * A[YY]) * x[Z_RANK1]) *
        inv_detA;

    b[Y_RANK1] = ((A[YZ] * A[ZX] - A[YX] * A[ZZ]) * x[X_RANK1] + (A[XX] * A[ZZ] - A[XZ] * A[ZY]) * x[Y_RANK1] +
        (A[YX] * A[ZY] - A[YY] * A[ZX]) * x[Z_RANK1]) *
        inv_detA;

    b[Z_RANK1] = ((A[YX] * A[ZY] - A[YY] * A[ZX]) * x[X_RANK1] + (A[XY] * A[XY] - A[ZX] * A[XX]) * x[Y_RANK1] +
        (A[XX] * A[YY] - A[XY] * A[YX]) * x[Z_RANK1]) *
        inv_detA;
  }

  double vecnorm_sq3(const double* x) {
    return (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
  }

  bool isoparameteric_coordinates_for_point_3d(
    sierra::nalu::LagrangeBasis& basis,
    const double* POINTER_RESTRICT elemNodalCoords,
    const double* POINTER_RESTRICT pointCoord,
    double* POINTER_RESTRICT isoParCoord,
    std::array<double,3> initialGuess,
    int maxIter,
    double tol,
    double deltaLimit)
  {
    int nNodes = basis.num_nodes();
    constexpr int dim = 3;
    std::array<double, dim> guess = initialGuess;
    std::array<double, dim> delta;
    int iter = 0;

    do {
      // interpolate coordinate at guess
      const auto& weights = basis.point_interpolation_weights(guess.data());

      std::array<double, dim> error_vec;
      error_vec[0] = pointCoord[0] - ddot(weights.data(), elemNodalCoords + 0 * nNodes, nNodes);
      error_vec[1] = pointCoord[1] - ddot(weights.data(), elemNodalCoords + 1 * nNodes, nNodes);
      error_vec[2] = pointCoord[2] - ddot(weights.data(), elemNodalCoords + 2 * nNodes, nNodes);

      // update guess along gradient of mapping from physical-to-reference coordinates

      // transpose of the jacobian of the forward mapping
      const auto& deriv = basis.point_derivative_weights(guess.data());
      std::array<double, dim * dim> jact{};
      for(int j = 0; j < nNodes; ++j) {
        jact[0] += deriv[0 + j * dim] * elemNodalCoords[j + 0 * nNodes];
        jact[1] += deriv[0 + j * dim] * elemNodalCoords[j + 1 * nNodes];
        jact[2] += deriv[0 + j * dim] * elemNodalCoords[j + 2 * nNodes];

        jact[3] += deriv[1 + j * dim] * elemNodalCoords[j + 0 * nNodes];
        jact[4] += deriv[1 + j * dim] * elemNodalCoords[j + 1 * nNodes];
        jact[5] += deriv[1 + j * dim] * elemNodalCoords[j + 2 * nNodes];

        jact[6] += deriv[2 + j * dim] * elemNodalCoords[j + 0 * nNodes];
        jact[7] += deriv[2 + j * dim] * elemNodalCoords[j + 1 * nNodes];
        jact[8] += deriv[2 + j * dim] * elemNodalCoords[j + 2 * nNodes];
      }

      // apply its inverse on the error vector
      action_of_inverse33(jact.data(), error_vec.data(), delta.data());

      // update guess.  Break if update is running away & report failure
      if (vecnorm_sq3(delta.data()) > deltaLimit) {
        iter = maxIter;
        break;
      }

      guess[0] += delta[0];
      guess[1] += delta[1];
      guess[2] += delta[2];

    } while (vecnorm_sq3(delta.data()) > tol && (++iter < maxIter));

    isoParCoord[0] = guess[0];
    isoParCoord[1] = guess[1];
    isoParCoord[2] = guess[2];

    return (iter < maxIter);
  }

  double determinant22(const double* mat)
  {
    enum { XX = 0, XY = 1, YX = 2, YY = 3 };
    return (mat[XX] * mat[YY] - mat[XY] * mat[YX]);
  }

  void action_of_inverse22(
    const double* POINTER_RESTRICT A,
    const double* POINTER_RESTRICT x,
    double*  POINTER_RESTRICT b)
  {
    ThrowAssert(determinant22(A) > tiny_positive_value());
    enum { XX = 0, XY = 1, YX = 2, YY = 3 };
    enum { X_RANK1 = 0, Y_RANK1 = 1};

    const double inv_detA = 1.0 / determinant22(A);
    b[X_RANK1] =  (A[YY] * x[X_RANK1] - A[XY] * x[Y_RANK1]) * inv_detA;
    b[Y_RANK1] = -(A[YX] * x[X_RANK1] - A[XX] * x[Y_RANK1]) * inv_detA;
  }

  double vecnorm_sq2(const double* x) {
    return (x[0] * x[0] + x[1] * x[1]);
  }

  bool isoparameteric_coordinates_for_point_2d(
    sierra::nalu::LagrangeBasis& basis,
    const double* POINTER_RESTRICT elemNodalCoords,
    const double* POINTER_RESTRICT pointCoord,
    double* POINTER_RESTRICT isoParCoord,
    std::array<double,2> initialGuess,
    int maxIter,
    double tol,
    double deltaLimit)
  {
    int nNodes = basis.num_nodes();
    constexpr int dim = 2;
    std::array<double, dim> guess = initialGuess;
    std::array<double, dim> delta;
    int iter = 0;

    do {
      // interpolate coordinate at guess
      const auto& weights = basis.point_interpolation_weights(guess.data());

      std::array<double, dim> error_vec;
      error_vec[0] = pointCoord[0] - ddot(weights.data(), elemNodalCoords + 0 * nNodes, nNodes);
      error_vec[1] = pointCoord[1] - ddot(weights.data(), elemNodalCoords + 1 * nNodes, nNodes);

      // update guess along gradient of mapping from physical-to-reference coordinates

      // transpose of the jacobian of the forward mapping
      const auto& deriv = basis.point_derivative_weights(guess.data());
      std::array<double, dim * dim> jact{};
      for(int j = 0; j < nNodes; ++j) {
        jact[0] += deriv[0 + j * dim] * elemNodalCoords[j + 0 * nNodes];
        jact[1] += deriv[0 + j * dim] * elemNodalCoords[j + 1 * nNodes];

        jact[2] += deriv[1 + j * dim] * elemNodalCoords[j + 0 * nNodes];
        jact[3] += deriv[1 + j * dim] * elemNodalCoords[j + 1 * nNodes];
      }

      // apply its inverse on the error vector
      action_of_inverse22(jact.data(), error_vec.data(), delta.data());

      // update guess.  Break if update is running away & report failure
      if (vecnorm_sq2(delta.data()) > deltaLimit) {
        iter = maxIter;
        break;
      }

      guess[0] += delta[0];
      guess[1] += delta[1];
    } while (vecnorm_sq2(delta.data()) > tol && (++iter < maxIter));

    isoParCoord[0] = guess[0];
    isoParCoord[1] = guess[1];

    return (iter < maxIter);
  }
}
}
