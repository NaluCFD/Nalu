/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporatlion.                                   */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#include <master_element/MasterElementUtils.h>

#include <element_promotion/LagrangeBasis.h>
#include <master_element/TensorOps.h>

#include <NaluEnv.h>
#include <FORTRAN_Proto.h>

#include <stk_util/environment/ReportHandler.hpp>

#include <array>
#include <limits>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace sierra{
namespace nalu{

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
        jact[1] += deriv[1 + j * dim] * elemNodalCoords[j + 0 * nNodes];
        jact[2] += deriv[2 + j * dim] * elemNodalCoords[j + 0 * nNodes];

        jact[3] += deriv[0 + j * dim] * elemNodalCoords[j + 1 * nNodes];
        jact[4] += deriv[1 + j * dim] * elemNodalCoords[j + 1 * nNodes];
        jact[5] += deriv[2 + j * dim] * elemNodalCoords[j + 1 * nNodes];

        jact[6] += deriv[0 + j * dim] * elemNodalCoords[j + 2 * nNodes];
        jact[7] += deriv[1 + j * dim] * elemNodalCoords[j + 2 * nNodes];
        jact[8] += deriv[2 + j * dim] * elemNodalCoords[j + 2 * nNodes];
      }

      // apply its inverse on the error vector
      solve33(jact.data(), error_vec.data(), delta.data());

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
        jact[1] += deriv[1 + j * dim] * elemNodalCoords[j + 0 * nNodes];
        jact[2] += deriv[0 + j * dim] * elemNodalCoords[j + 1 * nNodes];
        jact[3] += deriv[1 + j * dim] * elemNodalCoords[j + 1 * nNodes];
      }

      // apply its inverse on the error vector
      solve22(jact.data(), error_vec.data(), delta.data());

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

void general_grad_op(SharedMemView<DoubleType***>& deriv,
                 SharedMemView<DoubleType**>& coordel,
                 SharedMemView<DoubleType***>& gradop)
{
  size_t nint = deriv.dimension(0);
  size_t npe = deriv.dimension(1);

  DoubleType dx_ds1, dx_ds2, dx_ds3;
  DoubleType dy_ds1, dy_ds2, dy_ds3;
  DoubleType dz_ds1, dz_ds2, dz_ds3;
  DoubleType ds1_dx, ds1_dy, ds1_dz;
  DoubleType ds2_dx, ds2_dy, ds2_dz;
  DoubleType ds3_dx, ds3_dy, ds3_dz;

  double realmin = 2.2250738585072014e-308;

  for(size_t ki=0; ki<nint; ++ki) {
    dx_ds1 = 0.0;
    dx_ds2 = 0.0;
    dx_ds3 = 0.0;
    dy_ds1 = 0.0;
    dy_ds2 = 0.0;
    dy_ds3 = 0.0;
    dz_ds1 = 0.0;
    dz_ds2 = 0.0;
    dz_ds3 = 0.0;
//
// calculate the jacobian at the integration station -
//
    for(size_t kn=0; kn<npe; ++kn) {
      dx_ds1 = dx_ds1+deriv(ki,kn,0)*coordel(kn,0);
      dx_ds2 = dx_ds2+deriv(ki,kn,1)*coordel(kn,0);
      dx_ds3 = dx_ds3+deriv(ki,kn,2)*coordel(kn,0);

      dy_ds1 = dy_ds1+deriv(ki,kn,0)*coordel(kn,1);
      dy_ds2 = dy_ds2+deriv(ki,kn,1)*coordel(kn,1);
      dy_ds3 = dy_ds3+deriv(ki,kn,2)*coordel(kn,1);

      dz_ds1 = dz_ds1+deriv(ki,kn,0)*coordel(kn,2);
      dz_ds2 = dz_ds2+deriv(ki,kn,1)*coordel(kn,2);
      dz_ds3 = dz_ds3+deriv(ki,kn,2)*coordel(kn,2);
    }

//
// calculate the determinate of the jacobian at the integration station
//
    DoubleType det_j = dx_ds1*( dy_ds2*dz_ds3 - dz_ds2*dy_ds3 )
                     + dy_ds1*( dz_ds2*dx_ds3 - dx_ds2*dz_ds3 )
                     + dz_ds1*( dx_ds2*dy_ds3 - dy_ds2*dx_ds3 );
//
// protect against a negative or small value for the determinate of the
// jacobian. The value of real_min represents
// the smallest Real value (based upon the precision set for this
// compilation) which the machine can represent -
//
    DoubleType test = stk::math::if_then_else(det_j > 1.e+6*realmin, det_j, 1.0);
    DoubleType denom = 1.0/test;
//
// compute the gradient operators at the integration station -
//
    ds1_dx = denom*(dy_ds2*dz_ds3 - dz_ds2*dy_ds3);
    ds2_dx = denom*(dz_ds1*dy_ds3 - dy_ds1*dz_ds3);
    ds3_dx = denom*(dy_ds1*dz_ds2 - dz_ds1*dy_ds2);

    ds1_dy = denom*(dz_ds2*dx_ds3 - dx_ds2*dz_ds3);
    ds2_dy = denom*(dx_ds1*dz_ds3 - dz_ds1*dx_ds3);
    ds3_dy = denom*(dz_ds1*dx_ds2 - dx_ds1*dz_ds2);

    ds1_dz = denom*(dx_ds2*dy_ds3 - dy_ds2*dx_ds3);
    ds2_dz = denom*(dy_ds1*dx_ds3 - dx_ds1*dy_ds3);
    ds3_dz = denom*(dx_ds1*dy_ds2 - dy_ds1*dx_ds2);

    for(size_t kn=0; kn<npe; ++kn) {
      gradop(ki,kn,0) =
          deriv(ki,kn,0)*ds1_dx
        + deriv(ki,kn,1)*ds2_dx
        + deriv(ki,kn,2)*ds3_dx;

      gradop(ki,kn,1) =
          deriv(ki,kn,0)*ds1_dy
        + deriv(ki,kn,1)*ds2_dy
        + deriv(ki,kn,2)*ds3_dy;

      gradop(ki,kn,2) =
          deriv(ki,kn,0)*ds1_dz
        + deriv(ki,kn,1)*ds2_dz
        + deriv(ki,kn,2)*ds3_dz;
    }
  }
}

void threeD_gij(int npe, int nint,
                const SharedMemView<DoubleType***>& deriv,
                const SharedMemView<DoubleType**>& cordel,
                SharedMemView<DoubleType***>& gupperij,
                SharedMemView<DoubleType***>& glowerij)
{
  DoubleType dx_ds[3][3], ds_dx[3][3];
  DoubleType det_j = 0.0, denom = 0.0;
  const double realmin = 2.2250738585072014e-308;

  for(int ki=0; ki<nint; ++ki) {
    dx_ds[0][0] = 0.0;
    dx_ds[0][1] = 0.0;
    dx_ds[0][2] = 0.0;
    dx_ds[1][0] = 0.0;
    dx_ds[1][1] = 0.0;
    dx_ds[1][2] = 0.0;
    dx_ds[2][0] = 0.0;
    dx_ds[2][1] = 0.0;
    dx_ds[2][2] = 0.0;
// 
// calculate the jacobian at the integration station -
    for(int kn=0; kn<npe; ++kn) {
       dx_ds[0][0] += deriv(ki,kn,0)*cordel(kn,0);
       dx_ds[1][0] += deriv(ki,kn,1)*cordel(kn,0);
       dx_ds[2][0] += deriv(ki,kn,2)*cordel(kn,0);
//
       dx_ds[0][1] += deriv(ki,kn,0)*cordel(kn,1);
       dx_ds[1][1] += deriv(ki,kn,1)*cordel(kn,1);
       dx_ds[2][1] += deriv(ki,kn,2)*cordel(kn,1);
//
       dx_ds[0][2] += deriv(ki,kn,0)*cordel(kn,2);
       dx_ds[1][2] += deriv(ki,kn,1)*cordel(kn,2);
       dx_ds[2][2] += deriv(ki,kn,2)*cordel(kn,2);
    }

// calculate the determinate of the Jacobian at the integration station -
    det_j= dx_ds[0][0]*(dx_ds[1][1]*dx_ds[2][2]-dx_ds[1][2]*dx_ds[2][1])
         + dx_ds[0][1]*(dx_ds[1][2]*dx_ds[2][0]-dx_ds[1][0]*dx_ds[2][2])
         + dx_ds[0][2]*(dx_ds[1][0]*dx_ds[2][1]-dx_ds[1][1]*dx_ds[2][0]);

// clip
    denom = stk::math::if_then_else( det_j <= 1.e+6*realmin, 1.0, 1.0/det_j );

// calculate the inverse Jacobian
    ds_dx[0][0]= denom*(dx_ds[1][1]*dx_ds[2][2]-dx_ds[1][2]*dx_ds[2][1]);
    ds_dx[1][0]= denom*(dx_ds[1][2]*dx_ds[2][0]-dx_ds[1][0]*dx_ds[2][2]);
    ds_dx[2][0]= denom*(dx_ds[1][0]*dx_ds[2][1]-dx_ds[1][1]*dx_ds[2][0]);
//
    ds_dx[0][1]= denom*(dx_ds[0][2]*dx_ds[2][1]-dx_ds[0][1]*dx_ds[2][2]);
    ds_dx[1][1]= denom*(dx_ds[0][0]*dx_ds[2][2]-dx_ds[0][2]*dx_ds[2][0]);
    ds_dx[2][1]= denom*(dx_ds[0][1]*dx_ds[2][0]-dx_ds[0][0]*dx_ds[2][1]);
//
    ds_dx[0][2]= denom*(dx_ds[0][1]*dx_ds[1][2]-dx_ds[0][2]*dx_ds[1][1]);
    ds_dx[1][2]= denom*(dx_ds[0][2]*dx_ds[1][0]-dx_ds[0][0]*dx_ds[1][2]);
    ds_dx[2][2]= denom*(dx_ds[0][0]*dx_ds[1][1]-dx_ds[0][1]*dx_ds[1][0]);
//
    for(int i=0; i<3; ++i) {
       for(int j=0; j<3; ++j) {
          gupperij(ki,j,i) =
              dx_ds[0][i]*dx_ds[0][j]+dx_ds[1][i]*dx_ds[1][j]
              +dx_ds[2][i]*dx_ds[2][j];
         glowerij(ki,j,i) =
              ds_dx[i][0]*ds_dx[j][0]+ds_dx[i][1]*ds_dx[j][1]
              +ds_dx[i][2]*ds_dx[j][2];
       }
    }
  }
}

}
}
