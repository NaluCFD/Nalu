/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <master_element/Hex8FEM.h>
#include <master_element/MasterElement.h>

#include <FORTRAN_Proto.h>

#include <cmath>
#include <iostream>

namespace sierra{
namespace nalu{

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Hex8FEM::Hex8FEM()
  : MasterElement()
{
  nDim_ = 3;
  nodesPerElement_ = 8;
  numIntPoints_ = 8;

  // weights; -1:1
  weights_.resize(numIntPoints_);
  for ( int k = 0; k < numIntPoints_; ++k )
    weights_[k] = 1.0;

  // standard integration location +/ sqrt(3)/3
  intgLoc_.resize(24);
  double c[2] = {std::sqrt(3.0)/3.0,-std::sqrt(3.0)/3.0};
  
  // tensor product
  int l = 0;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k, ++l) {
        intgLoc_[l*3] = c[i];
        intgLoc_[l*3+1] = c[j];
        intgLoc_[l*3+2] = c[k];
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Hex8FEM::~Hex8FEM()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- determinant ------------------ n/a ------------------------------
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
//-------- grad_op ---------------------------------------------------------
//--------------------------------------------------------------------------
void Hex8FEM::grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  hex8_fem_derivative(numIntPoints_, &intgLoc_[0], deriv);
  
  SIERRA_FORTRAN(hex_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );
 
}

//--------------------------------------------------------------------------
//-------- face_grad_op ----------------------------------------------------
//--------------------------------------------------------------------------
void Hex8FEM::face_grad_op(
  const int nelem,
  const int face_ordinal,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  int lerr = 0;
  int npf = 4;
  int ndim = 3;
  
  const int nface = 1;
  double dpsi[24];
  double grad[24];
  
  for ( int n=0; n<nelem; n++ ) {
    
    for ( int k=0; k<npf; k++ ) {

      const int row = 12*face_ordinal + k*ndim;
      hex8_fem_derivative
        ( nface,
          &intgExpFace_[row], dpsi );
      
      SIERRA_FORTRAN(hex_gradient_operator)
        ( &nface,
          &nodesPerElement_,
          &nface,
          dpsi,
          &coords[24*n], grad, &det_j[npf*n+k], error, &lerr );
      
      if ( lerr )
        std::cout << "sorry, issue with face_grad_op.." << std::endl;
      
      for ( int j=0; j<24; j++) {
        gradop[k*nelem*24+n*24+j] = grad[j];
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- general_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Hex8FEM::general_shape_fcn(
  const int numIp,
  const double *isoParCoord,
  double *shpfc) {

  hex8_fem_shape_fcn(numIp,isoParCoord,shpfc);
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Hex8FEM::shape_fcn(double *shpfc)
{
  hex8_fem_shape_fcn(numIntPoints_,&intgLoc_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- hex8_fem_shape_fcn ----------------------------------------------
//--------------------------------------------------------------------------
void
Hex8FEM::hex8_fem_shape_fcn(
  const int  &numIp,
  const double *isoParCoord, 
  double *shpfc)
{
  // -1:1 isoparametric range
  const double npe = nodesPerElement_;
  for ( int ip = 0; ip < numIp; ++ip ) {
    
    const int rowIpc = 3*ip;
    const int rowSfc = npe*ip;
    
    const double s1 = isoParCoord[rowIpc];
    const double s2 = isoParCoord[rowIpc+1];
    const double s3 = isoParCoord[rowIpc+2];
    shpfc[rowSfc  ] = 0.125*(1.0-s1)*(1.0-s2)*(1.0-s3);
    shpfc[rowSfc+1] = 0.125*(1.0+s1)*(1.0-s2)*(1.0-s3);
    shpfc[rowSfc+2] = 0.125*(1.0+s1)*(1.0+s2)*(1.0-s3);
    shpfc[rowSfc+3] = 0.125*(1.0-s1)*(1.0+s2)*(1.0-s3);
    shpfc[rowSfc+4] = 0.125*(1.0-s1)*(1.0-s2)*(1.0+s3);
    shpfc[rowSfc+5] = 0.125*(1.0+s1)*(1.0-s2)*(1.0+s3);
    shpfc[rowSfc+6] = 0.125*(1.0+s1)*(1.0+s2)*(1.0+s3);
    shpfc[rowSfc+7] = 0.125*(1.0-s1)*(1.0+s2)*(1.0+s3);
  }
}

//--------------------------------------------------------------------------
//-------- hex8_fem_derivative ---------------------------------------------
//--------------------------------------------------------------------------
void
Hex8FEM::hex8_fem_derivative(
  const int npt, const double* par_coord,
  double* deriv)
{
  for (int i = 0; i < npt; ++i) {
    deriv[i*nodesPerElement_*3    ] = -0.125*(1.0-par_coord[i*3+1])*(1.0-par_coord[i*3+2]);
    deriv[i*nodesPerElement_*3 + 3] =  0.125*(1.0-par_coord[i*3+1])*(1.0-par_coord[i*3+2]);
    deriv[i*nodesPerElement_*3 + 6] =  0.125*(1.0+par_coord[i*3+1])*(1.0-par_coord[i*3+2]);
    deriv[i*nodesPerElement_*3 + 9] = -0.125*(1.0+par_coord[i*3+1])*(1.0-par_coord[i*3+2]);
    deriv[i*nodesPerElement_*3 + 12  ] = -0.125*(1.0-par_coord[i*3+1])*(1.0+par_coord[i*3+2]);
    deriv[i*nodesPerElement_*3 + 15] =  0.125*(1.0-par_coord[i*3+1])*(1.0+par_coord[i*3+2]);
    deriv[i*nodesPerElement_*3 + 18] =  0.125*(1.0+par_coord[i*3+1])*(1.0+par_coord[i*3+2]);
    deriv[i*nodesPerElement_*3 + 21] = -0.125*(1.0+par_coord[i*3+1])*(1.0+par_coord[i*3+2]);

    deriv[i*nodesPerElement_*3 + 1] = -0.125*(1.0-par_coord[i*3])*(1.0-par_coord[i*3+2]);
    deriv[i*nodesPerElement_*3 + 4] = -0.125*(1.0+par_coord[i*3])*(1.0-par_coord[i*3+2]);
    deriv[i*nodesPerElement_*3 + 7] =  0.125*(1.0+par_coord[i*3])*(1.0-par_coord[i*3+2]);
    deriv[i*nodesPerElement_*3 + 10] =  0.125*(1.0-par_coord[i*3])*(1.0-par_coord[i*3+2]);
    deriv[i*nodesPerElement_*3 + 13] = -0.125*(1.0-par_coord[i*3])*(1.0+par_coord[i*3+2]);
    deriv[i*nodesPerElement_*3 + 16] = -0.125*(1.0+par_coord[i*3])*(1.0+par_coord[i*3+2]);
    deriv[i*nodesPerElement_*3 + 19] =  0.125*(1.0+par_coord[i*3])*(1.0+par_coord[i*3+2]);
    deriv[i*nodesPerElement_*3 + 22] =  0.125*(1.0-par_coord[i*3])*(1.0+par_coord[i*3+2]);

    deriv[i*nodesPerElement_*3 + 2] = -0.125*(1.0-par_coord[i*3])*(1.0-par_coord[i*3+1]);
    deriv[i*nodesPerElement_*3 + 5] = -0.125*(1.0+par_coord[i*3])*(1.0-par_coord[i*3+1]);
    deriv[i*nodesPerElement_*3 + 8] =  -0.125*(1.0+par_coord[i*3])*(1.0+par_coord[i*3+1]);
    deriv[i*nodesPerElement_*3 + 11] =  -0.125*(1.0-par_coord[i*3])*(1.0+par_coord[i*3+1]);
    deriv[i*nodesPerElement_*3 + 14] = 0.125*(1.0-par_coord[i*3])*(1.0-par_coord[i*3+1]);
    deriv[i*nodesPerElement_*3 + 17] = 0.125*(1.0+par_coord[i*3])*(1.0-par_coord[i*3+1]);
    deriv[i*nodesPerElement_*3 + 20] =  0.125*(1.0+par_coord[i*3])*(1.0+par_coord[i*3+1]);
    deriv[i*nodesPerElement_*3 + 23] =  0.125*(1.0-par_coord[i*3])*(1.0+par_coord[i*3+1]);
  }                            
}

}
}
