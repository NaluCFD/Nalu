/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef HermitePolynomialInterpolation_h
#define HermitePolynomialInterpolation_h

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

namespace sierra{
namespace nalu{

class HermitePolynomialInterpolation
{
public:

  HermitePolynomialInterpolation() {}
  virtual ~HermitePolynomialInterpolation() {}
  virtual void do_hermite( 
    const double *elemNodalQ, 
    const double *elemNodalCoords, 
    const double *elemNodalDqdx, 
    const double *haloCoord,
    double &interpResult) = 0;
};

class HermitePolynomialInterpolationFourPoint : public HermitePolynomialInterpolation
{
public:

  HermitePolynomialInterpolationFourPoint();
  virtual ~HermitePolynomialInterpolationFourPoint() {}
  void do_hermite( 
    const double *elemNodalQ, 
    const double *elemNodalCoords, 
    const double *elemNodalDqdx, 
    const double *haloCoord,
    double &interpResult);

  int numberPoints_;
  Teuchos::SerialDenseMatrix<int, double> A_;
  Teuchos::SerialDenseVector<int, double> b_;
  Teuchos::LAPACK<int, double> lapack_;

};

class HermitePolynomialInterpolationEightPoint : public HermitePolynomialInterpolation
{
public:

  HermitePolynomialInterpolationEightPoint();
  virtual ~HermitePolynomialInterpolationEightPoint() {}
  void do_hermite(
    const double *elemNodalQ,
    const double *elemNodalCoords,
    const double *elemNodalDqdx,
    const double *haloCoord,
    double &interpResult);

  int numberPoints_;
  Teuchos::SerialDenseMatrix<int, double> A_;
  Teuchos::SerialDenseVector<int, double> b_;
  Teuchos::LAPACK<int, double> lapack_;

};

} // namespace nalu
} // namespace Sierra

#endif
