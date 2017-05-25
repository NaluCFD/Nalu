/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef MasterElementUtils_h
#define MasterElementUtils_h

#include <master_element/MasterElement.h>
#include <element_promotion/TensorProductQuadratureRule.h>
#include <element_promotion/LagrangeBasis.h>

#include <element_promotion/ElementDescription.h>
#include <element_promotion/HexNElementDescription.h>
#include <element_promotion/QuadNElementDescription.h>

#include <vector>
#include <array>
#include <limits>

#ifdef __INTEL_COMPILER
#define POINTER_RESTRICT restrict
#else
#define POINTER_RESTRICT __restrict__
#endif

namespace sierra{
namespace nalu{

  inline constexpr double tiny_positive_value() {
    return (1.0e6*std::numeric_limits<double>::min());
  }

  double ddot(const double* u, const double* v, int n);

  // computes b = A^-1 x;
  void action_of_inverse22(
    const double* POINTER_RESTRICT A,
    const double* POINTER_RESTRICT x,
    double*  POINTER_RESTRICT b);

  // computes b = A^-1 x;
  void action_of_inverse33(
    const double* POINTER_RESTRICT A,
    const double* POINTER_RESTRICT x,
    double*  POINTER_RESTRICT b);

  bool isoparameteric_coordinates_for_point_3d(
    sierra::nalu::LagrangeBasis& basis,
    const double* POINTER_RESTRICT elemNodalCoords,
    const double* POINTER_RESTRICT pointCoord,
    double* POINTER_RESTRICT isoParCoord,
    std::array<double,3> initialGuess,
    int maxIter,
    double tolerance,
    double deltaLimit = 1.0e4
  );

  bool isoparameteric_coordinates_for_point_2d(
    sierra::nalu::LagrangeBasis& basis,
    const double* POINTER_RESTRICT elemNodalCoords,
    const double* POINTER_RESTRICT pointCoord,
    double* POINTER_RESTRICT isoParCoord,
    std::array<double,2> initialGuess,
    int maxIter,
    double tolerance,
    double deltaLimit = 1.0e4
  );

}
}

#endif
