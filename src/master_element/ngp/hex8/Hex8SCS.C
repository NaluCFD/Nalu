/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "master_element/ngp/hex8/Hex8SCS.h"
#include "master_element/MasterElement.h"
#include "AlgTraits.h"
#include "NaluEnv.h"

#include <cmath>
#include <limits>

namespace sierra{
namespace nalu{


//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Hex8SCS::Hex8SCS()
  : HexSCS()
{
  // nothing to do as it is all in HexSCS
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Hex8SCS::~Hex8SCS()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- hex8_derivative -------------------------------------------------
//--------------------------------------------------------------------------
void Hex8SCS::hex8_derivative(
  const int npts,
  SharedMemView<DoubleType**> &intgLoc,
  SharedMemView<DoubleType***> &deriv)
{
  const double half = 0.50;
  const double one4th = 0.25;
  for (int  j = 0; j < npts; ++j) {

    const double s1 = intgLoc_(j,0);
    const double s2 = intgLoc_(j,1);
    const double s3 = intgLoc_(j,2);
    const double s1s2 = s1*s2;
    const double s2s3 = s2*s3;
    const double s1s3 = s1*s3;

    // shape function derivative in the s1 direction -
    deriv(j,0,0) = half*( s3 + s2 ) - s2s3 - one4th;
    deriv(j,1,0) = half*(-s3 - s2 ) + s2s3 + one4th;
    deriv(j,2,0) = half*(-s3 + s2 ) - s2s3 + one4th;
    deriv(j,3,0) = half*(+s3 - s2 ) + s2s3 - one4th;
    deriv(j,4,0) = half*(-s3 + s2 ) + s2s3 - one4th;
    deriv(j,5,0) = half*(+s3 - s2 ) - s2s3 + one4th;
    deriv(j,6,0) = half*(+s3 + s2 ) + s2s3 + one4th;
    deriv(j,7,0) = half*(-s3 - s2 ) - s2s3 - one4th;

    // shape function derivative in the s2 direction -
    deriv(j,0,1) = half*( s3 + s1 ) - s1s3 - one4th;
    deriv(j,1,1) = half*( s3 - s1 ) + s1s3 - one4th;
    deriv(j,2,1) = half*(-s3 + s1 ) - s1s3 + one4th;
    deriv(j,3,1) = half*(-s3 - s1 ) + s1s3 + one4th;
    deriv(j,4,1) = half*(-s3 + s1 ) + s1s3 - one4th;
    deriv(j,5,1) = half*(-s3 - s1 ) - s1s3 - one4th;
    deriv(j,6,1) = half*( s3 + s1 ) + s1s3 + one4th;
    deriv(j,7,1) = half*( s3 - s1 ) - s1s3 + one4th;
    
    // shape function derivative in the s3 direction -
    deriv(j,0,2) = half*( s2 + s1 ) - s1s2 - one4th;
    deriv(j,1,2) = half*( s2 - s1 ) + s1s2 - one4th;
    deriv(j,2,2) = half*(-s2 - s1 ) - s1s2 - one4th;
    deriv(j,3,2) = half*(-s2 + s1 ) + s1s2 - one4th;
    deriv(j,4,2) = half*(-s2 - s1 ) + s1s2 + one4th;
    deriv(j,5,2) = half*(-s2 + s1 ) - s1s2 + one4th;
    deriv(j,6,2) = half*( s2 + s1 ) + s1s2 + one4th;
    deriv(j,7,2) = half*( s2 - s1 ) - s1s2 + one4th;
  }
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Hex8SCS::shape_fcn(SharedMemView<DoubleType**> &shpfc)
{
  hex_shape_fcn(&numIntPoints_, &intgLoc_, shpfc);
}

//--------------------------------------------------------------------------
//-------- hex8_shape_fcn --------------------------------------------------
//--------------------------------------------------------------------------
void
Hex8SCS::hex8_shape_fcn(
  const int  &npts,
  const double *isoParCoord, 
  SharedMemView<DoubleType**> &shape_fcn)
{
  const double one8th = 0.125;
  const double half = 0.50;
  for ( int j = 0; j < numIntPoints_; ++j ) {
    
    const double s1 = isoParCoord(j,0);
    const double s2 = isoParCoord(j,1);
    const double s3 = isoParCoord(j,2);
    
    shape_fcn(j,0) = one8th + one4th*(-s1 - s2 - s3)
      + half*( s2*s3 + s3*s1 + s1*s2 ) - s1*s2*s3;
    shape_fcn(j,1) = one8th + one4th*( s1 - s2 - s3)
      + half*( s2*s3 - s3*s1 - s1*s2 ) + s1*s2*s3;
    shape_fcn(j,2) = one8th + one4th*( s1 + s2 - s3)
      + half*(-s2*s3 - s3*s1 + s1*s2 ) - s1*s2*s3;
    shape_fcn(j,3) = one8th + one4th*(-s1 + s2 - s3)
      + half*(-s2*s3 + s3*s1 - s1*s2 ) + s1*s2*s3;
    shape_fcn(j,4) = one8th + one4th*(-s1 - s2 + s3)
      + half*(-s2*s3 - s3*s1 + s1*s2 ) + s1*s2*s3;
    shape_fcn(j,5) = one8th + one4th*( s1 - s2 + s3)
      + half*(-s2*s3 + s3*s1 - s1*s2 ) - s1*s2*s3;
    shape_fcn(j,6) = one8th + one4th*( s1 + s2 + s3)
      + half*( s2*s3 + s3*s1 + s1*s2 ) + s1*s2*s3;
    shape_fcn(j,7) = one8th + one4th*(-s1 + s2 + s3)
      + half*( s2*s3 - s3*s1 - s1*s2 ) - s1*s2*s3;    
  }
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Hex8SCS::shifted_shape_fcn(SharedMemView<DoubleType**> &shpfc)
{
  hex_shape_fcn(&numIntPoints_, &intgLocShift_[0], shpfc);
}

}
}
