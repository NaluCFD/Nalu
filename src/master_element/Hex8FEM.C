/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <master_element/Hex8FEM.h>
#include <master_element/MasterElement.h>
#include <master_element/MasterElementFunctions.h>
#include <master_element/TensorOps.h>

#include <FORTRAN_Proto.h>
#include <NaluEnv.h>

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
  const double gIP = std::sqrt(3.0)/3.0;
  intgLoc_[0]  = -gIP; intgLoc_[1]  = -gIP; intgLoc_[2]  = -gIP; 
  intgLoc_[3]  = +gIP; intgLoc_[4]  = -gIP; intgLoc_[5]  = -gIP; 
  intgLoc_[6]  = +gIP; intgLoc_[7]  = +gIP; intgLoc_[8]  = -gIP;
  intgLoc_[9]  = -gIP; intgLoc_[10] = +gIP; intgLoc_[11] = -gIP;
  intgLoc_[12] = -gIP; intgLoc_[13] = -gIP; intgLoc_[14] = +gIP;
  intgLoc_[15] = +gIP; intgLoc_[16] = -gIP; intgLoc_[17] = +gIP;
  intgLoc_[18] = +gIP; intgLoc_[19] = +gIP; intgLoc_[20] = +gIP;
  intgLoc_[21] = -gIP; intgLoc_[22] = +gIP; intgLoc_[23] = +gIP;

  // shifted to nodes (Gauss Lobatto)
  intgLocShift_.resize(24);
  const double glIP = 1.0;
  intgLocShift_[0]  = -glIP; intgLocShift_[1]  = -glIP; intgLocShift_[2]  = -glIP; 
  intgLocShift_[3]  = +glIP; intgLocShift_[4]  = -glIP; intgLocShift_[5]  = -glIP; 
  intgLocShift_[6]  = +glIP; intgLocShift_[7]  = +glIP; intgLocShift_[8]  = -glIP;
  intgLocShift_[9]  = -glIP; intgLocShift_[10] = +glIP; intgLocShift_[11] = -glIP;
  intgLocShift_[12] = -glIP; intgLocShift_[13] = -glIP; intgLocShift_[14] = +glIP;
  intgLocShift_[15] = +glIP; intgLocShift_[16] = -glIP; intgLocShift_[17] = +glIP;
  intgLocShift_[18] = +glIP; intgLocShift_[19] = +glIP; intgLocShift_[20] = +glIP;
  intgLocShift_[21] = -glIP; intgLocShift_[22] = +glIP; intgLocShift_[23] = +glIP;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Hex8FEM::~Hex8FEM()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- determinant_fem -------------------------------------------------
//--------------------------------------------------------------------------
void Hex8FEM::determinant_fem(
  SharedMemView<DoubleType**>&coords,
  SharedMemView<DoubleType***>&deriv,
  SharedMemView<DoubleType*>&det_j)
{
  hex8_fem_derivative(numIntPoints_, &intgLoc_[0], deriv);
  generic_determinant_3d<AlgTraitsHex8>(deriv, coords, det_j);
}

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
//-------- grad_op ---------------------------------------------------------
//--------------------------------------------------------------------------
void Hex8FEM::grad_op_fem(
  SharedMemView<DoubleType**>&coords,
  SharedMemView<DoubleType***>&gradop,
  SharedMemView<DoubleType***>&deriv,
  SharedMemView<DoubleType*>&det_j)
{
  hex8_fem_derivative(numIntPoints_, &intgLoc_[0], deriv);
  generic_grad_op_fem<AlgTraitsHex8>(deriv, coords, gradop, det_j);
}

//--------------------------------------------------------------------------
//-------- shifted_grad_op -------------------------------------------------
//--------------------------------------------------------------------------
void Hex8FEM::shifted_grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  hex8_fem_derivative(numIntPoints_, &intgLocShift_[0], deriv);
  
  SIERRA_FORTRAN(hex_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr ); 
}


//--------------------------------------------------------------------------
//-------- shifted_grad_op -------------------------------------------------
//--------------------------------------------------------------------------
void Hex8FEM::shifted_grad_op_fem(
  SharedMemView<DoubleType**>&coords,
  SharedMemView<DoubleType***>&gradop,
  SharedMemView<DoubleType***>&deriv,
  SharedMemView<DoubleType*>&det_j)
{
  hex8_fem_derivative(numIntPoints_, &intgLocShift_[0], deriv);
  generic_grad_op_fem<AlgTraitsHex8>(deriv, coords, gradop, det_j);
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
        NaluEnv::self().naluOutput() << "sorry, issue with face_grad_op.." << std::endl;
      
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
  double *shpfc)
{
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
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Hex8FEM::shifted_shape_fcn(double *shpfc)
{
  hex8_fem_shape_fcn(numIntPoints_,&intgLocShift_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Hex8FEM::shape_fcn(SharedMemView<DoubleType**> &shpfc)
{
  hex8_fem_shape_fcn(numIntPoints_,&intgLoc_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Hex8FEM::shifted_shape_fcn(SharedMemView<DoubleType**> &shpfc)
{
  hex8_fem_shape_fcn(numIntPoints_,&intgLocShift_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- gij -------------------------------------------------------------
//--------------------------------------------------------------------------
void Hex8FEM::gij(
  const double *coords,
  double *gupperij,
  double *glowerij,
  double *deriv)
{
  hex8_fem_derivative(numIntPoints_, &intgLoc_[0], deriv);
  SIERRA_FORTRAN(threed_gij)
    ( &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gupperij, glowerij);
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
  const int npe = nodesPerElement_;
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
//-------- hex8_fem_shape_fcn ----------------------------------------------
//--------------------------------------------------------------------------
void
Hex8FEM::hex8_fem_shape_fcn(
  const int  &numIp,
  const double *isoParCoord,
  SharedMemView<DoubleType**> shpfc)
{
  // -1:1 isoparametric range
  for ( int ip = 0; ip < numIp; ++ip ) {
    const int rowIpc = 3*ip;
    const DoubleType s1 = isoParCoord[rowIpc+0];
    const DoubleType s2 = isoParCoord[rowIpc+1];
    const DoubleType s3 = isoParCoord[rowIpc+2];
    shpfc(ip,0) = 0.125*(1.0-s1)*(1.0-s2)*(1.0-s3);
    shpfc(ip,1) = 0.125*(1.0+s1)*(1.0-s2)*(1.0-s3);
    shpfc(ip,2) = 0.125*(1.0+s1)*(1.0+s2)*(1.0-s3);
    shpfc(ip,3) = 0.125*(1.0-s1)*(1.0+s2)*(1.0-s3);
    shpfc(ip,4) = 0.125*(1.0-s1)*(1.0-s2)*(1.0+s3);
    shpfc(ip,5) = 0.125*(1.0+s1)*(1.0-s2)*(1.0+s3);
    shpfc(ip,6) = 0.125*(1.0+s1)*(1.0+s2)*(1.0+s3);
    shpfc(ip,7) = 0.125*(1.0-s1)*(1.0+s2)*(1.0+s3);
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
    deriv[i*nodesPerElement_*3 + 12] = -0.125*(1.0-par_coord[i*3+1])*(1.0+par_coord[i*3+2]);
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

void
Hex8FEM::hex8_fem_derivative(
  const int npt, const double* par_coord,
  SharedMemView<DoubleType***> deriv)
{
  for (int ip = 0; ip < npt; ++ip) {
    DoubleType s1 = par_coord[ip*3+0];
    DoubleType s2 = par_coord[ip*3+1];
    DoubleType s3 = par_coord[ip*3+2];

    deriv(ip,0,0) = -0.125*(1.0-s2)*(1.0-s3);
    deriv(ip,1,0) =  0.125*(1.0-s2)*(1.0-s3);
    deriv(ip,2,0) =  0.125*(1.0+s2)*(1.0-s3);
    deriv(ip,3,0) = -0.125*(1.0+s2)*(1.0-s3);
    deriv(ip,4,0) = -0.125*(1.0-s2)*(1.0+s3);
    deriv(ip,5,0) =  0.125*(1.0-s2)*(1.0+s3);
    deriv(ip,6,0) =  0.125*(1.0+s2)*(1.0+s3);
    deriv(ip,7,0) = -0.125*(1.0+s2)*(1.0+s3);

    deriv(ip,0,1) = -0.125*(1.0-s1)*(1.0-s3);
    deriv(ip,1,1) = -0.125*(1.0+s1)*(1.0-s3);
    deriv(ip,2,1) =  0.125*(1.0+s1)*(1.0-s3);
    deriv(ip,3,1) =  0.125*(1.0-s1)*(1.0-s3);
    deriv(ip,4,1) = -0.125*(1.0-s1)*(1.0+s3);
    deriv(ip,5,1) = -0.125*(1.0+s1)*(1.0+s3);
    deriv(ip,6,1) =  0.125*(1.0+s1)*(1.0+s3);
    deriv(ip,7,1) =  0.125*(1.0-s1)*(1.0+s3);

    deriv(ip,0,2) = -0.125*(1.0-s1)*(1.0-s2);
    deriv(ip,1,2) = -0.125*(1.0+s1)*(1.0-s2);
    deriv(ip,2,2) = -0.125*(1.0+s1)*(1.0+s2);
    deriv(ip,3,2) = -0.125*(1.0-s1)*(1.0+s2);
    deriv(ip,4,2) =  0.125*(1.0-s1)*(1.0-s2);
    deriv(ip,5,2) =  0.125*(1.0+s1)*(1.0-s2);
    deriv(ip,6,2) =  0.125*(1.0+s1)*(1.0+s2);
    deriv(ip,7,2) =  0.125*(1.0-s1)*(1.0+s2);
  }
}

}
}
