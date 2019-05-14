/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "master_element/MasterElement.h"
#include "master_element/Tri6FEM.h"

namespace sierra{
namespace nalu{


//-------- tri6_derivative -------------------------------------------------
template <typename DerivType>
void tri6_derivative(
  const int npts,
  const double *intgLoc,
  DerivType &deriv)
{
  for (int ip = 0; ip < npts; ++ip) {
    DoubleType s = intgLoc[ip*2+0];
    DoubleType t = intgLoc[ip*2+1];
    
    /* deriv(ip,ic,j) */
    deriv(ip,0,0) = -1.0*(1.0-2.0*s-2.0*t) - 2.0*(1.0-s-t);
    deriv(ip,1,0) = (2.0*s-1.0) + 2.0*s;
    deriv(ip,2,0) = 0.0;
    deriv(ip,3,0) = 4.0*(1.0-s-t) - 4.0*s;
    deriv(ip,4,0) = 4.0*t;
    deriv(ip,5,0) = -4.0*t;

    deriv(ip,0,1) = -1.0*(1.0-2.0*s-2.0*t) - 2.0*(1.0-s-t);
    deriv(ip,1,1) = 0.0;
    deriv(ip,2,1) = (2.0*t-1.0) + 2.0*t;
    deriv(ip,3,1) = -4.0*s;
    deriv(ip,4,1) = 4.0*s;
    deriv(ip,5,1) = 4.0*(1.0-s-t) - 4.0*t;
  }
}
  
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Tri6FEM::Tri6FEM()  
  : MasterElement()
{
  nDim_ = 3;
  nodesPerElement_ = 6;
  numIntPoints_ = 7; // quadratic, degree 5

  // weights
  weights_.resize(numIntPoints_);
  weights_[0] = 0.1125000000000000;
  weights_[1] = 0.0629695902724135;
  weights_[2] = 0.0629695902724135;
  weights_[3] = 0.0629695902724135;
  weights_[4] = 0.066197076394253;
  weights_[5] = 0.066197076394253;
  weights_[6] = 0.066197076394253;

  // standard integration location
  intgLoc_.resize(numIntPoints_*(nDim_-1));    
  intgLoc_[0]  = 0.333333333333333; intgLoc_[1]  = 0.333333333333333;
  intgLoc_[2]  = 0.101286507323456; intgLoc_[3]  = 0.101286507323456;
  intgLoc_[4]  = 0.797426985353087; intgLoc_[5]  = 0.101286507323456;
  intgLoc_[6]  = 0.101286507323456; intgLoc_[7]  = 0.797426985353087;
  intgLoc_[8]  = 0.470142064105115; intgLoc_[9]  = 0.470142064105115;
  intgLoc_[10] = 0.059715871789770; intgLoc_[11] = 0.470142064105115;
  intgLoc_[12] = 0.470142064105115; intgLoc_[13] = 0.059715871789770;
  
  // shifted; n/a
} 

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Tri6FEM::~Tri6FEM()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- determinant_fem -------------------------------------------------
//--------------------------------------------------------------------------
void Tri6FEM::determinant_fem(
  SharedMemView<DoubleType**>&coords,
  SharedMemView<DoubleType***>&deriv,
  SharedMemView<DoubleType*>&det_j)
{
  tri6_derivative(numIntPoints_, &intgLoc_[0], deriv);

  for (int ip = 0; ip < numIntPoints_; ++ip) {
    
    NALU_ALIGNED DoubleType dxds[3]= {0.0, 0.0, 0.0};
    NALU_ALIGNED DoubleType dxdt[3]= {0.0, 0.0, 0.0};
   
    for ( int n = 0; n < nodesPerElement_; ++n ) {
      const DoubleType derivS = deriv(ip,n,0);
      dxds[0] += derivS*coords(n,0);
      dxds[1] += derivS*coords(n,1);
      dxds[2] += derivS*coords(n,2);

      const DoubleType derivT = deriv(ip,n,1);
      dxdt[0] += derivT*coords(n,0);
      dxdt[1] += derivT*coords(n,1);
      dxdt[2] += derivT*coords(n,2);
    }
    
    DoubleType detXY =  dxds[0]*dxdt[1] - dxdt[0]*dxds[1];
    DoubleType detYZ =  dxds[1]*dxdt[2] - dxdt[1]*dxds[2];
    DoubleType detXZ = -dxds[0]*dxdt[2] + dxdt[0]*dxds[2];
    
    det_j(ip) = stk::math::sqrt(detXY*detXY + detYZ*detYZ + detXZ*detXZ);
  }
}

//--------------------------------------------------------------------------
//-------- normal ----------------------------------------------------------
//--------------------------------------------------------------------------
void Tri6FEM::normal_fem(
  SharedMemView<DoubleType**>&coords,
  SharedMemView<DoubleType***>&deriv,
  SharedMemView<DoubleType**>&normal)
{
  tri6_derivative(numIntPoints_, &intgLoc_[0], deriv);  

  for (int ip = 0; ip < numIntPoints_; ++ip) {
    
    NALU_ALIGNED DoubleType dxds[3]= {0.0, 0.0, 0.0};
    NALU_ALIGNED DoubleType dxdt[3]= {0.0, 0.0, 0.0};
   
    for ( int n = 0; n < nodesPerElement_; ++n ) {
      const DoubleType derivS = deriv(ip,n,0);
      dxds[0] += derivS*coords(n,0);
      dxds[1] += derivS*coords(n,1);
      dxds[2] += derivS*coords(n,2);

      const DoubleType derivT = deriv(ip,n,1);
      dxdt[0] += derivT*coords(n,0);
      dxdt[1] += derivT*coords(n,1);
      dxdt[2] += derivT*coords(n,2);
    }  
    DoubleType detXY =  dxds[0]*dxdt[1] - dxdt[0]*dxds[1];
    DoubleType detYZ =  dxds[1]*dxdt[2] - dxdt[1]*dxds[2];
    DoubleType detXZ = -dxds[0]*dxdt[2] + dxdt[0]*dxds[2];
    
    DoubleType det = stk::math::sqrt(detXY*detXY + detYZ*detYZ + detXZ*detXZ);
    
    normal(ip,0) = detYZ/det;
    normal(ip,1) = detXZ/det;
    normal(ip,2) = detXY/det;
  }
}
  
//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Tri6FEM::shape_fcn(SharedMemView<DoubleType**> &shpfc)
{
  tri6_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Tri6FEM::shape_fcn(double *shpfc)
{
  tri6_shape_fcn(numIntPoints_,&intgLoc_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- tri6_shape_fcn --------------------------------------------------
//--------------------------------------------------------------------------
void
Tri6FEM::tri6_shape_fcn(
  const int  &npts,
  const double *isoParCoord,
  SharedMemView<DoubleType**> &shpfc)
{
  for ( int ip = 0; ip < npts; ++ip ) {
    const int rowIpc = 2*ip;

    const DoubleType s = isoParCoord[rowIpc+0];
    const DoubleType t = isoParCoord[rowIpc+1];
    
    shpfc(ip,0) = (1.0-s-t)*(1.0-2.0*s-2.0*t);
    shpfc(ip,1) = s*(2.0*s-1.0);
    shpfc(ip,2) = t*(2.0*t-1.0); 
    shpfc(ip,3) = 4.0*s*(1-s-t);;
    shpfc(ip,4) = 4.0*s*t; 
    shpfc(ip,5) = 4.0*t*(1.0-s-t);
  }
}

//--------------------------------------------------------------------------
//-------- tri6_shape_fcn --------------------------------------------------
//--------------------------------------------------------------------------
void
Tri6FEM::tri6_shape_fcn(
  const int  &npts,
  const double *isoParCoord,
  double *shpfc)
{
  const int npe = nodesPerElement_;
  for ( int ip = 0; ip < npts; ++ip ) {
    const int rowIpc = 2*ip;
    const int rowSfc = npe*ip;

    const double s = isoParCoord[rowIpc+0];
    const double t = isoParCoord[rowIpc+1];
    
    shpfc[rowSfc  ] = (1.0-s-t)*(1.0-2.0*s-2.0*t);
    shpfc[rowSfc+1] = s*(2.0*s-1.0);
    shpfc[rowSfc+2] = t*(2.0*t-1.0); 
    shpfc[rowSfc+3] = 4.0*s*(1-s-t);;
    shpfc[rowSfc+4] = 4.0*s*t; 
    shpfc[rowSfc+5] = 4.0*t*(1.0-s-t);
  }
}

} // namespace nalu
} // namespace sierra
