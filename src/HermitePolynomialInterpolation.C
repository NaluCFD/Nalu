/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <HermitePolynomialInterpolation.h>

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// HermitePolynomialInterpolationFourPoint - Hermite for a four node sample
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
HermitePolynomialInterpolationFourPoint::HermitePolynomialInterpolationFourPoint() :
  HermitePolynomialInterpolation(),
  numberPoints_(4),
  A_(12,12),
  b_(12)
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- do_hermite ------------------------------------------------------
//--------------------------------------------------------------------------
void
HermitePolynomialInterpolationFourPoint::do_hermite(
  const double *elemNodalQ, 
  const double *elemNodalCoords, 
  const double *elemNodalDqdx, 
  const double *haloCoord,
  double &interpResult)
{

  // load the matrix
  for ( int ni = 0 ; ni < numberPoints_; ++ni ) {

    const double xC = elemNodalCoords[0*numberPoints_+ni];
    const double yC = elemNodalCoords[1*numberPoints_+ni];

    // contributions from F(x,y)
    A_(ni,0)  = 1.0;
    A_(ni,1)  = xC;
    A_(ni,2)  = yC;
    A_(ni,3)  = xC*yC;
    A_(ni,4)  = xC*xC;
    A_(ni,5)  = yC*yC;
    A_(ni,6)  = yC*xC*xC;
    A_(ni,7)  = xC*yC*yC;
    A_(ni,8)  = xC*xC*xC;
    A_(ni,9)  = yC*yC*yC;
    A_(ni,10) = yC*xC*xC*xC;
    A_(ni,11) = xC*yC*yC*yC;
    b_(ni) = elemNodalQ[ni];

    // contributions from Fx(xj,yj)
    const int offSetFx = ni+numberPoints_;
    A_(offSetFx,0)  = 0.0;
    A_(offSetFx,1)  = 1.0;
    A_(offSetFx,2)  = 0.0;
    A_(offSetFx,3)  = yC;
    A_(offSetFx,4)  = 2.0*xC;
    A_(offSetFx,5)  = 0.0;
    A_(offSetFx,6)  = 2.0*xC*yC;
    A_(offSetFx,7)  = yC*yC;
    A_(offSetFx,8)  = 3.0*xC*xC;
    A_(offSetFx,9)  = 0.0;
    A_(offSetFx,10) = 3.0*yC*xC*xC;
    A_(offSetFx,11) = yC*yC*yC; 
    b_(offSetFx) = elemNodalDqdx[0*numberPoints_+ni];

    // contributions from Fy(xj,yj)
    const int offSetFy = ni+2*numberPoints_;
    A_(offSetFy,0)  = 0.0;
    A_(offSetFy,1)  = 0.0;
    A_(offSetFy,2)  = 1.0;
    A_(offSetFy,3)  = xC;
    A_(offSetFy,4)  = 0.0;
    A_(offSetFy,5)  = 2.0*yC;
    A_(offSetFy,6)  = xC*xC;
    A_(offSetFy,7)  = 2.0*xC*yC;
    A_(offSetFy,8)  = 0.0;
    A_(offSetFy,9)  = 3.0*yC*yC;
    A_(offSetFy,10) = xC*xC*xC;
    A_(offSetFy,11) = 3.0*xC*yC*yC; 
    b_(offSetFy) = elemNodalDqdx[1*numberPoints_+ni];
    
  }

  // Perform an LU factorization of this matrix. 
  int ipiv[12], info;
  char TRANS = 'N';
  lapack_.GETRF( 12, 12, A_.values(), A_.stride(), ipiv, &info ); 
  
  // Solve the linear system; solution "x" saved off in "b"
  lapack_.GETRS( TRANS, 12, 1, A_.values(), A_.stride(), 
		ipiv, b_.values(), b_.stride(), &info );  
  
  // extract the value
  const double xH = haloCoord[0];
  const double yH = haloCoord[1];

  interpResult = b_(0) + b_(1)*xH + b_(2)*yH + b_(3)*xH*yH 
    + b_(4)*xH*xH + b_(5)*yH*yH + b_(6)*yH*xH*xH + b_(7)*xH*yH*yH 
    + b_(8)*xH*xH*xH + b_(9)*yH*yH*yH + b_(10)*yH*xH*xH*xH + b_(11)*xH*yH*yH*yH;

}

//==========================================================================
// Class Definition
//==========================================================================
// HermitePolynomialInterpolationEightPoint - Hermite for a seight node sample
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
HermitePolynomialInterpolationEightPoint::HermitePolynomialInterpolationEightPoint() :
  HermitePolynomialInterpolation(),
  numberPoints_(8),
  A_(32,32),
  b_(32)
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- do_hermite ------------------------------------------------------
//--------------------------------------------------------------------------
void
HermitePolynomialInterpolationEightPoint::do_hermite(
  const double *elemNodalQ,
  const double *elemNodalCoords,
  const double *elemNodalDqdx,
  const double *haloCoord,
  double &interpResult)
{

  // load the matrix
  for ( int ni = 0 ; ni < numberPoints_; ++ni ) {

    const double xC = elemNodalCoords[0*numberPoints_+ni];
    const double yC = elemNodalCoords[1*numberPoints_+ni];
    const double zC = elemNodalCoords[2*numberPoints_+ni];

    // save square and cube
    const double xCsq = xC*xC;
    const double yCsq = yC*yC;
    const double zCsq = zC*zC;
    const double xCcb = xC*xC*xC;
    const double yCcb = yC*yC*yC;
    const double zCcb = zC*zC*zC;

    // contributions from F(x,y)
    A_(ni,0)  = 1.0;
    A_(ni,1)  = xC;
    A_(ni,2)  = yC;
    A_(ni,3)  = zC;
    A_(ni,4)  = xC*yC;
    A_(ni,5)  = xC*zC;
    A_(ni,6)  = yC*zC;
    A_(ni,7)  = xC*yC*zC;
    A_(ni,8)  = xCsq;
    A_(ni,9)  = yCsq;
    A_(ni,10) = zCsq;
    A_(ni,11) = yC*xCsq;
    A_(ni,12) = zC*xCsq;
    A_(ni,13) = yC*zC*xCsq;
    A_(ni,14) = xC*yCsq;
    A_(ni,15) = zC*yCsq;
    A_(ni,16) = xC*zC*yCsq;
    A_(ni,17) = xC*zCsq;
    A_(ni,18) = yC*zCsq;
    A_(ni,19) = xC*yC*zCsq;
    A_(ni,20) = xCcb;
    A_(ni,21) = yCcb;
    A_(ni,22) = zCcb;
    A_(ni,23) = yC*xCcb;
    A_(ni,24) = zC*xCcb;
    A_(ni,25) = yC*zC*xCcb;
    A_(ni,26) = xC*yCcb;
    A_(ni,27) = zC*yCcb;
    A_(ni,28) = xC*zC*yCcb;
    A_(ni,29) = xC*zCcb;
    A_(ni,30) = yC*zCcb;
    A_(ni,31) = xC*yC*zCcb;
    b_(ni) = elemNodalQ[ni];

    // contributions from Fx(xj,yj,zj)
    const int offSetFx = ni+numberPoints_;
    A_(offSetFx,0)  = 0.0;
    A_(offSetFx,1)  = 1.0;
    A_(offSetFx,2)  = 0.0;
    A_(offSetFx,3)  = 0.0;
    A_(offSetFx,4)  = yC;
    A_(offSetFx,5)  = zC;
    A_(offSetFx,6)  = 0.0;
    A_(offSetFx,7)  = yC*zC;
    A_(offSetFx,8)  = 2.0*xC;
    A_(offSetFx,9)  = 0.0;
    A_(offSetFx,10) = 0.0;
    A_(offSetFx,11) = yC*2.0*xC;
    A_(offSetFx,12) = zC*2.0*xC;
    A_(offSetFx,13) = yC*zC*2.0*xC;
    A_(offSetFx,14) = yCsq;
    A_(offSetFx,15) = 0.0;
    A_(offSetFx,16) = zC*yCsq;
    A_(offSetFx,17) = zCsq;
    A_(offSetFx,18) = 0.0;
    A_(offSetFx,19) = yC*zCsq;
    A_(offSetFx,20) = 3.0*xCsq;
    A_(offSetFx,21) = 0.0;
    A_(offSetFx,22) = 0.0;
    A_(offSetFx,23) = yC*3.0*xCsq;
    A_(offSetFx,24) = zC*3.0*xCsq;
    A_(offSetFx,25) = yC*zC*3.0*xCsq;
    A_(offSetFx,26) = yCcb;
    A_(offSetFx,27) = 0.0;
    A_(offSetFx,28) = zC*yCcb;
    A_(offSetFx,29) = zCcb;
    A_(offSetFx,30) = 0.0;
    A_(offSetFx,31) = yC*zCcb;
    b_(offSetFx) = elemNodalDqdx[0*numberPoints_+ni];

    // contributions from Fy(xj,yj,zj)
    const int offSetFy = ni+2*numberPoints_;
    A_(offSetFy,0)  = 0.0;
    A_(offSetFy,1)  = 0.0;
    A_(offSetFy,2)  = 1.0;
    A_(offSetFy,3)  = 0.0;
    A_(offSetFy,4)  = xC;
    A_(offSetFy,5)  = 0.0;
    A_(offSetFy,6)  = zC;
    A_(offSetFy,7)  = xC*zC;
    A_(offSetFy,8)  = 0.0;
    A_(offSetFy,9)  = 2.0*yC;
    A_(offSetFy,10) = 0.0;
    A_(offSetFy,11) = xCsq;
    A_(offSetFy,12) = 0.0;
    A_(offSetFy,13) = zC*xCsq;
    A_(offSetFy,14) = xC*2.0*yC;
    A_(offSetFy,15) = zC*2.0*yC;
    A_(offSetFy,16) = xC*zC*2.0*yC;
    A_(offSetFy,17) = 0.0;
    A_(offSetFy,18) = zCsq;
    A_(offSetFy,19) = xC*zCsq;
    A_(offSetFy,20) = 0.0;
    A_(offSetFy,21) = 3.0*yCsq;
    A_(offSetFy,22) = 0.0;
    A_(offSetFy,23) = xCcb;
    A_(offSetFy,24) = 0.0;
    A_(offSetFy,25) = zC*xCcb;
    A_(offSetFy,26) = xC*3.0*yCsq;
    A_(offSetFy,27) = zC*3.0*yCsq;
    A_(offSetFy,28) = xC*zC*3.0*yCsq;
    A_(offSetFy,29) = 0.0;
    A_(offSetFy,30) = zCcb;
    A_(offSetFy,31) = xC*zCcb;
    b_(offSetFy) = elemNodalDqdx[1*numberPoints_+ni];

    // contributions from Fz(xj,yj,zj)
    const int offSetFz = ni+3*numberPoints_;
    A_(offSetFz,0)  = 0.0;
    A_(offSetFz,1)  = 0.0;
    A_(offSetFz,2)  = 0.0;
    A_(offSetFz,3)  = 1.0;
    A_(offSetFz,4)  = 0.0;
    A_(offSetFz,5)  = xC;
    A_(offSetFz,6)  = yC;
    A_(offSetFz,7)  = xC*yC;
    A_(offSetFz,8)  = 0.0;
    A_(offSetFz,9)  = 0.0;
    A_(offSetFz,10) = 2.0*zC;
    A_(offSetFz,11) = 0.0;
    A_(offSetFz,12) = xCsq;
    A_(offSetFz,13) = yC*xCsq;
    A_(offSetFz,14) = 0.0;
    A_(offSetFz,15) = yCsq;
    A_(offSetFz,16) = xC*yCsq;
    A_(offSetFz,17) = xC*2.0*zC;
    A_(offSetFz,18) = yC*2.0*zC;
    A_(offSetFz,19) = xC*yC*2.0*zC;
    A_(offSetFz,20) = 0.0;
    A_(offSetFz,21) = 0.0;
    A_(offSetFz,22) = 3.0*zCsq;
    A_(offSetFz,23) = 0.0;
    A_(offSetFz,24) = xCcb;
    A_(offSetFz,25) = yC*xCcb;
    A_(offSetFz,26) = 0.0;
    A_(offSetFz,27) = yCcb;
    A_(offSetFz,28) = xC*yCcb;
    A_(offSetFz,29) = xC*3.0*zCsq;
    A_(offSetFz,30) = yC*3.0*zCsq;
    A_(offSetFz,31) = xC*yC*3.0*zCsq;
    b_(offSetFz) = elemNodalDqdx[2*numberPoints_+ni];

  }

  // Perform an LU factorization of this matrix.
  int ipiv[32], info;
  char TRANS = 'N';
  lapack_.GETRF( 32, 32, A_.values(), A_.stride(), ipiv, &info );

  // Solve the linear system; solution "x" saved off in "b"
  lapack_.GETRS( TRANS, 32, 1, A_.values(), A_.stride(),
    ipiv, b_.values(), b_.stride(), &info );

  // extract the value
  const double xH = haloCoord[0];
  const double yH = haloCoord[1];
  const double zH = haloCoord[2];
  const double xHsq = xH*xH;
  const double yHsq = yH*yH;
  const double zHsq = zH*zH;
  const double xHcb = xH*xH*xH;
  const double yHcb = yH*yH*yH;
  const double zHcb = zH*zH*zH;

  interpResult = b_(0) + b_(1)*xH + b_(2)*yH + b_(3)*zH + b_(4)*xH*yH
    + b_(5)*xH*zH + b_(6)*yH*zH + b_(7)*xH*yH*zH + b_(8)*xHsq + b_(9)*yHsq + b_(10)*zHsq
    + b_(11)*yH*xHsq + b_(12)*zH*xHsq + b_(13)*yH*zH*xHsq + b_(14)*xH*yHsq + b_(15)*zH*yHsq + b_(16)*xH*zH*yHsq
    + b_(17)*xH*zHsq + b_(18)*yH*zHsq + b_(19)*xH*yH*zHsq + b_(20)*xHcb + b_(21)*yHcb + b_(22)*zHcb
    + b_(23)*yH*xHcb + b_(24)*zH*xHcb + b_(25)*yH*zH*xHcb + b_(26)*xH*yHcb + b_(27)*zH*yHcb + b_(28)*xH*zH*yHcb
    + b_(29)*xH*zHcb + b_(30)*yH*zHcb + b_(31)*xH*yH*zHcb;
}



} // namespace nalu
} // namespace Sierra
