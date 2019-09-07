/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <master_element/Tet10FEM.h>
#include <master_element/MasterElement.h>
#include <master_element/MasterElementFunctions.h>
#include <master_element/TensorOps.h>

#include "FORTRAN_Proto.h"
#include "AlgTraits.h"
#include "NaluEnv.h"

#include <cmath>
#include <iostream>

namespace sierra{
namespace nalu{

//-------- tet10_deriv -------------------------------------------------
template <typename DerivType>
void tet10_deriv(
  const int npts,
  const double *intgLoc,
  DerivType &deriv)
{
  for (int ip = 0; ip < npts; ++ip) {
    DoubleType s1 = intgLoc[ip*3+0];
    DoubleType s2 = intgLoc[ip*3+1];
    DoubleType s3 = intgLoc[ip*3+2];
    
    /* deriv(ip,ic,j) */
    deriv(ip,0,0) = 4.0*(s1+s2+s3)-3.0;
    deriv(ip,1,0) = 4.0*s1-1.0;
    deriv(ip,2,0) = 0.0;
    deriv(ip,3,0) = 0.0;
    deriv(ip,4,0) = 4.0*(1.0-2.0*s1-s2-s3);
    deriv(ip,5,0) = 4.0*s2;
    deriv(ip,6,0) =-4.0*s2;
    deriv(ip,7,0) =-4.0*s3;
    deriv(ip,8,0) = 4.0*s3;
    deriv(ip,9,0) = 0.0;

    deriv(ip,0,1) = 4.0*(s1+s2+s3)-3.0;
    deriv(ip,1,1) = 0.0;
    deriv(ip,2,1) = 4.0*s2-1.0;
    deriv(ip,3,1) = 0.0;
    deriv(ip,4,1) =-4.0*s1;
    deriv(ip,5,1) = 4.0*s1;
    deriv(ip,6,1) = 4.0*(1.0-s1-2.0*s2-s3);
    deriv(ip,7,1) =-4.0*s3;
    deriv(ip,8,1) = 0.0;
    deriv(ip,9,1) = 4.0*s3;

    deriv(ip,0,2) = 4.0*(s1+s2+s3)-3.0;
    deriv(ip,1,2) = 0.0;
    deriv(ip,2,2) = 0.0;
    deriv(ip,3,2) = 4.0*s3-1.0;
    deriv(ip,4,2) =-4.0*s1;
    deriv(ip,5,2) = 0.0;
    deriv(ip,6,2) =-4.0*s2;
    deriv(ip,7,2) = 4.0*(1.0-s1-s2-2.0*s3);
    deriv(ip,8,2) = 4.0*s1;
    deriv(ip,9,2) = 4.0*s2;
  }
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Tet10FEM::Tet10FEM()
  : MasterElement(),
    tri6FEM_(new Tri6FEM())
{
  nDim_ = 3;
  nodesPerElement_ = 10;
  
  if ( AlgTraitsTet10::numGp_ == 15 ) {
    // Moderate-degree tetrahedral quadrature formulas; Keast 1986 
    numIntPoints_ = 15;
    weights_.resize(15);
    intgLoc_.resize(45);
    intgLocShift_.resize(45);
    
    // weights
    weights_[0]  = 0.1817020685825351;
    weights_[1]  = 0.0361607142857143;
    weights_[2]  = 0.0361607142857143;
    weights_[3]  = 0.0361607142857143;
    weights_[4]  = 0.0361607142857143;
    weights_[5]  = 0.0698714945161738;
    weights_[6]  = 0.0698714945161738;
    weights_[7]  = 0.0698714945161738;
    weights_[8]  = 0.0698714945161738;
    weights_[9]  = 0.0656948493683187;
    weights_[10] = 0.0656948493683187;
    weights_[11] = 0.0656948493683187;
    weights_[12] = 0.0656948493683187;
    weights_[13] = 0.0656948493683187;
    weights_[14] = 0.0656948493683187;
    
    // standard integration location
    intgLoc_[0]  = 0.2500000000000000;  intgLoc_[1]  = 0.2500000000000000;  intgLoc_[2]  = 0.2500000000000000;
    intgLoc_[3]  = 0.0000000000000000;  intgLoc_[4]  = 0.3333333333333333;  intgLoc_[5]  = 0.3333333333333333;
    intgLoc_[6]  = 0.3333333333333333;  intgLoc_[7]  = 0.3333333333333333;  intgLoc_[8]  = 0.3333333333333333;
    intgLoc_[9]  = 0.3333333333333333;  intgLoc_[10] = 0.3333333333333333;  intgLoc_[11] = 0.0000000000000000;
    intgLoc_[12] = 0.3333333333333333;  intgLoc_[13] = 0.0000000000000000;  intgLoc_[14] = 0.3333333333333333;
    intgLoc_[15] = 0.7272727272727273;  intgLoc_[16] = 0.0909090909090909;  intgLoc_[17] = 0.0909090909090909;
    intgLoc_[18] = 0.0909090909090909;  intgLoc_[19] = 0.0909090909090909;  intgLoc_[20] = 0.0909090909090909;
    intgLoc_[21] = 0.0909090909090909;  intgLoc_[22] = 0.0909090909090909;  intgLoc_[23] = 0.7272727272727273;
    intgLoc_[24] = 0.0909090909090909;  intgLoc_[25] = 0.7272727272727273;  intgLoc_[26] = 0.0909090909090909;
    intgLoc_[27] = 0.4334498464263357;  intgLoc_[28] = 0.0665501535736643;  intgLoc_[29] = 0.0665501535736643;
    intgLoc_[30] = 0.0665501535736643;  intgLoc_[31] = 0.4334498464263357;  intgLoc_[32] = 0.0665501535736643;
    intgLoc_[33] = 0.0665501535736643;  intgLoc_[34] = 0.0665501535736643;  intgLoc_[35] = 0.4334498464263357;
    intgLoc_[36] = 0.0665501535736643;  intgLoc_[37] = 0.4334498464263357;  intgLoc_[38] = 0.4334498464263357;
    intgLoc_[39] = 0.4334498464263357;  intgLoc_[40] = 0.0665501535736643;  intgLoc_[41] = 0.4334498464263357;
    intgLoc_[42] = 0.4334498464263357;  intgLoc_[43] = 0.4334498464263357;  intgLoc_[44] = 0.0665501535736643;

    // shifted to nodes (Gauss Lobatto); n/a
  }
  else {
    // 16-pt rule 
    numIntPoints_ = 16;
    weights_.resize(16);
    intgLoc_.resize(48);    
    intgLocShift_.resize(48);

    // weights
    const double w1 = 8.395632350020469e-03;
    const double w2 = 1.109034477221540e-02;
    weights_[0]  = w1; weights_[1]  = w1; weights_[2]  = w1; weights_[3]  = w1;
    weights_[4]  = w2; weights_[5]  = w2; weights_[6]  = w2; weights_[7]  = w2;
    weights_[8]  = w2; weights_[9]  = w2; weights_[10] = w2; weights_[11] = w2;
    weights_[12] = w2; weights_[13] = w2; weights_[14] = w2; weights_[15] = w2;

    // standard integration location
    intgLoc_[0]  = 0.7716429020672371;  intgLoc_[1]  = 0.07611903264425430; intgLoc_[2]  = 0.07611903264425430;  
    intgLoc_[3]  = 0.07611903264425430; intgLoc_[4]  = 0.7716429020672371;  intgLoc_[5]  = 0.07611903264425430;
    intgLoc_[6]  = 0.07611903264425430; intgLoc_[7]  = 0.07611903264425430; intgLoc_[8]  = 0.7716429020672371;
    intgLoc_[9]  = 0.07611903264425430; intgLoc_[10] = 0.07611903264425430; intgLoc_[11] = 0.07611903264425430;
    intgLoc_[12] = 0.4042339134672644;  intgLoc_[13] = 0.07183164526766925; intgLoc_[14] = 0.1197005277978019;
    intgLoc_[15] = 0.4042339134672644;  intgLoc_[16] = 0.1197005277978019;  intgLoc_[17] = 0.07183164526766925;
    intgLoc_[18] = 0.1197005277978019;  intgLoc_[19] = 0.07183164526766925; intgLoc_[20] = 0.4042339134672644;
    intgLoc_[21] = 0.1197005277978019;  intgLoc_[22] = 0.4042339134672644;  intgLoc_[23] = 0.07183164526766925;
    intgLoc_[24] = 0.07183164526766925; intgLoc_[25] = 0.1197005277978019;  intgLoc_[26] = 0.4042339134672644;
    intgLoc_[27] = 0.07183164526766925; intgLoc_[28] = 0.4042339134672644;  intgLoc_[29] = 0.1197005277978019;
    intgLoc_[30] = 0.4042339134672644;  intgLoc_[31] = 0.07183164526766925; intgLoc_[32] = 0.4042339134672644;
    intgLoc_[33] = 0.07183164526766925; intgLoc_[34] = 0.4042339134672644;  intgLoc_[35] = 0.4042339134672644;
    intgLoc_[36] = 0.4042339134672644;  intgLoc_[37] = 0.4042339134672644;  intgLoc_[38] = 0.07183164526766925;
    intgLoc_[39] = 0.1197005277978019;  intgLoc_[40] = 0.4042339134672644;  intgLoc_[41] = 0.4042339134672644;
    intgLoc_[42] = 0.4042339134672644;  intgLoc_[43] = 0.1197005277978019;  intgLoc_[44] = 0.4042339134672644;
    intgLoc_[45] = 0.4042339134672644;  intgLoc_[46] = 0.4042339134672644;  intgLoc_[47] = 0.1197005277978019;

    // shifted to nodes (Gauss Lobatto); n/a
  }

  // exposed face; use helper functions
  const int numTri6Ip = tri6FEM_->numIntPoints_;
  const int numExposedFace = 4;
  
  // deal with exposed ips
  intgExpFace_.resize(nDim_*numTri6Ip*numExposedFace); // 3*7*4
  for ( int i = 0; i < numExposedFace; ++i ) 
    sidePcoords_to_elemPcoords(i,numTri6Ip, &tri6FEM_->intgLoc_[0], &intgExpFace_[i*nDim_*numTri6Ip]);

  // mapping from a side ordinal to the node ordinals on that side
  sideNodeOrdinals_ = {0, 1, 3, 4, 8, 7,
                       1, 2, 3, 5, 9, 8,
                       0, 3, 2, 7, 9, 6,
                       0, 2, 1, 6, 5, 4};
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Tet10FEM::~Tet10FEM()
{
  delete tri6FEM_;
}

//--------------------------------------------------------------------------
//-------- determinant_fem -------------------------------------------------
//--------------------------------------------------------------------------
void Tet10FEM::determinant_fem(
  SharedMemView<DoubleType**>&coords,
  SharedMemView<DoubleType***>&deriv,
  SharedMemView<DoubleType*>&det_j)
{
  tet10_deriv(numIntPoints_, &intgLoc_[0], deriv);
  generic_determinant_3d<AlgTraitsTet10>(deriv, coords, det_j);
}

//--------------------------------------------------------------------------
//-------- grad_op_fem -----------------------------------------------------
//--------------------------------------------------------------------------
void Tet10FEM::grad_op_fem(
  SharedMemView<DoubleType**>&coords,
  SharedMemView<DoubleType***>&gradop,
  SharedMemView<DoubleType***>&deriv,
  SharedMemView<DoubleType*>&det_j)
{
  tet10_deriv(numIntPoints_, &intgLoc_[0], deriv);
  generic_grad_op_fem<AlgTraitsTet10>(deriv, coords, gradop, det_j);
}

//--------------------------------------------------------------------------
//-------- shifted_grad_op_fem ---------------------------------------------
//--------------------------------------------------------------------------
void Tet10FEM::shifted_grad_op_fem(
  SharedMemView<DoubleType**>&coords,
  SharedMemView<DoubleType***>&gradop,
  SharedMemView<DoubleType***>&deriv,
  SharedMemView<DoubleType*>&det_j)
{
  tet10_deriv(numIntPoints_, &intgLocShift_[0], deriv);
  generic_grad_op_fem<AlgTraitsTet10>(deriv, coords, gradop, det_j);
}

//--------------------------------------------------------------------------
//-------- face_grad_op_fem ------------------------------------------------
//--------------------------------------------------------------------------
void Tet10FEM::face_grad_op_fem(
  int face_ordinal,
  SharedMemView<DoubleType**>& coords,
  SharedMemView<DoubleType***>& gradop,
  SharedMemView<DoubleType*>&det_j)
{
  // reminder for possible shifting
  const bool shifted = false;
  using traits = AlgTraitsTri6Tet10;

  const std::vector<double> &exp_face = shifted ? intgExpFaceShift_ : intgExpFace_;

  constexpr int derivSize = traits::numFaceIp_ * traits::nodesPerElement_ * traits::nDim_;
  DoubleType psi[derivSize];
  SharedMemView<DoubleType***> deriv(psi, traits::numFaceIp_, traits::nodesPerElement_, traits::nDim_);

  const int offset = traits::numFaceIp_ * traits::nDim_ * face_ordinal;
  tet10_deriv(traits::numFaceIp_, &exp_face[offset], deriv);
  generic_grad_op_fem<AlgTraitsTet10>(deriv, coords, gradop, det_j);
}

//--------------------------------------------------------------------------
//-------- gij -------------------------------------------------------------
//--------------------------------------------------------------------------
void Tet10FEM::gij(
  SharedMemView<DoubleType**>& coords,
  SharedMemView<DoubleType***>& gupper,
  SharedMemView<DoubleType***>& glower,
  SharedMemView<DoubleType***>& deriv)
{
  tet10_deriv(numIntPoints_, &intgLoc_[0], deriv);
  generic_gij_3d<AlgTraitsTet10>(deriv, coords, gupper, glower);
}

//--------------------------------------------------------------------------
//-------- side_node_ordinals ----------------------------------------------
//--------------------------------------------------------------------------
const int *
Tet10FEM::side_node_ordinals(int ordinal)
{
  // define face_ordinal->node_ordinal mappings for each face (ordinal);
  return &sideNodeOrdinals_[ordinal*6];
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Tet10FEM::shape_fcn(SharedMemView<DoubleType**> &shpfc)
{
  tet10_fem_shape_fcn(numIntPoints_,&intgLoc_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Tet10FEM::shifted_shape_fcn(SharedMemView<DoubleType**> &shpfc)
{
  tet10_fem_shape_fcn(numIntPoints_,&intgLocShift_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Tet10FEM::shape_fcn(double *shpfc)
{
  tet10_fem_shape_fcn(numIntPoints_,&intgLoc_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Tet10FEM::shifted_shape_fcn(double *shpfc)
{
  tet10_fem_shape_fcn(numIntPoints_,&intgLocShift_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- general_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Tet10FEM::general_shape_fcn(
  const int numIp,
  const double *isoParCoord,
  double *shpfc)
{
  tet10_fem_shape_fcn(numIp, isoParCoord, shpfc);
}

//--------------------------------------------------------------------------
//-------- gij -------------------------------------------------------------
//--------------------------------------------------------------------------
void Tet10FEM::gij(
  const double *coords,
  double *gupperij,
  double *glowerij,
  double *deriv)
{
  tet10_derivative(numIntPoints_, &intgLoc_[0], deriv);
  SIERRA_FORTRAN(threed_gij)
    ( &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gupperij, glowerij);
}

//--------------------------------------------------------------------------
//-------- isInElement -----------------------------------------------------
//--------------------------------------------------------------------------
double
Tet10FEM::isInElement(
  const double * elem_nodal_coor,
  const double * point_coor,
  double * par_coor )
{
  // Define these as compile time constants so the
  // compiler can unroll the loops.
  const int NCOEFF = 10;
  const int NCOORD = 3;
  const int MAX_ITER = 10;
  const double isInElemConverged = 1.0e-16;

  // Translate element so that (x,y) coordinates of the first node are (0,0)
  // (xp,yp) is the point that we're searching for (xi,eta), translate it too
  double xn[NCOORD][NCOEFF];
  double xp[NCOORD];
  double s_new[NCOORD];
  double s_cur[NCOORD];
  for(int i = 0; i < NCOORD; ++i) {
    const double x_ref = elem_nodal_coor[i * NCOEFF];
    xp[i] = point_coor[i] - x_ref;
    s_new[i] = 0.0;
    s_cur[i] = 0.0;
    for(int j = 0; j < NCOEFF; ++j) {
      xn[i][j] = elem_nodal_coor[i * NCOEFF + j] - x_ref;
    }
  }

  // Newton-Raphson iteration for (xi,eta)
  double jac[NCOORD][NCOORD];
  double jac_inv[NCOORD][NCOORD];
  double f[NCOORD];
  double shapefcn[NCOEFF];       // (npe,npts);
  double deriv[NCOORD * NCOEFF]; // (ncoord,npe,npts);

  int iter = 0;
  bool converged = false;

  do {
    for(int i = 0; i < NCOORD; ++i) {
      s_cur[i] = s_new[i];
    }

    // Evaluate the shape function and derivatives here.
    tet10_derivative(1, s_cur, deriv);
    general_shape_fcn(1, s_cur, shapefcn);

    for(int i = 0; i < NCOORD; ++i) {
      // Residual assembly:
      f[i] = -xp[i];
      for(int k = 0; k < NCOEFF; ++k) {
        f[i] += xn[i][k] * shapefcn[k];
      }

      // Jacobian assembly:
      for(int j = 0; j < NCOORD; ++j) {
        jac[j][i] = 0.0;
        for(int k = 0; k < NCOEFF; ++k) {
          jac[j][i] += deriv[k * NCOORD + i] * xn[j][k];
        }
      }
    }

    // invert the jacobian
    const double det_inv = 1.0 / (+jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) -
                                  jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0]) +
                                  jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]));
    
    jac_inv[0][0] = (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) * det_inv;
    jac_inv[0][1] = (jac[0][2] * jac[2][1] - jac[0][1] * jac[2][2]) * det_inv;
    jac_inv[0][2] = (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]) * det_inv;
    
    jac_inv[1][0] = (jac[1][2] * jac[2][0] - jac[1][0] * jac[2][2]) * det_inv;
    jac_inv[1][1] = (jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0]) * det_inv;
    jac_inv[1][2] = (jac[0][2] * jac[1][0] - jac[0][0] * jac[1][2]) * det_inv;
    
    jac_inv[2][0] = (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]) * det_inv;
    jac_inv[2][1] = (jac[0][1] * jac[2][0] - jac[0][0] * jac[2][1]) * det_inv;
    jac_inv[2][2] = (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]) * det_inv;
    
    // Compute the Newton correction and updated solution:
    double d_norm = 0.0;
    for(int i = 0; i < NCOORD; ++i) {
      double delta_i = 0.0;
      for(int j = 0; j < NCOORD; ++j) {
        delta_i -= jac_inv[i][j] * f[j];
      }
      s_new[i] = s_cur[i] + delta_i;
      d_norm += delta_i * delta_i;
    }

    // tookusa: The vector norm is actually the square of the L2 norm
    converged = d_norm < isInElemConverged;
  } while(!converged && ++iter < MAX_ITER);

  for(int i = 0; i < NCOORD; ++i) {
    par_coor[i] = std::numeric_limits<double>::max();
  }
  double dist = std::numeric_limits<double>::max();

  if(iter < MAX_ITER) {
    for(int i = 0; i < NCOORD; ++i) {
      par_coor[i] = s_new[i];
    }
    dist = parametric_distance(par_coor);
  }
  return dist;
}

//--------------------------------------------------------------------------
//-------- interpolatePoint ------------------------------------------------
//--------------------------------------------------------------------------
void
Tet10FEM::interpolatePoint(
  const int  &ncomp_field,
  const double *par_coord,
  const double *field,
  double *result )
{
  // 'field' is a flat array of dimension (10,ncomp_field) (Fortran ordering);
  double s1 = par_coord[0];
  double s2 = par_coord[1];
  double s3 = par_coord[2];

  for ( int i = 0; i < ncomp_field; i++ )
  {
    const int b = 10 * i;
    result[i] = (1.0-2.0*s1-2.0*s2-2.0*s3)*(1.0-s1-s2-s3)*field[b+0] 
      + s1*(2.0*s1-1.0)*field[b+1] 
      + s2*(2.0*s2-1.0)*field[b+2] 
      + s3*(2.0*s3-1.0)*field[b+3] 
      + 4.0*s1*(1.0-s1-s2-s3)*field[b+4] 
      + 4.0*s1*s2*field[b+5] 
      + 4.0*s2*(1.0-s1-s2-s3)*field[b+6] 
      + 4.0*s3*(1.0-s1-s2-s3)*field[b+7] 
      + 4.0*s1*s3*field[b+8] 
      + 4.0*s3*s2*field[b+9];
  }
}

//--------------------------------------------------------------------------
//-------- sidePcoords_to_elemPcoords --------------------------------------
//--------------------------------------------------------------------------
void
Tet10FEM::sidePcoords_to_elemPcoords(
  const int & side_ordinal,
  const int & npoints,
  const double *side_pcoords,
  double *elem_pcoords)
{
  switch(side_ordinal) {
  case 0:
    for(int i = 0; i < npoints; i++) {
      elem_pcoords[i*3+0] = side_pcoords[2*i+0];
      elem_pcoords[i*3+1] = 0.0;
      elem_pcoords[i*3+2] = side_pcoords[2*i+1];
    }
    break;
  case 1:
    for(int i = 0; i < npoints; i++) {
      elem_pcoords[i*3+0] = 1.0 - side_pcoords[2*i+0] - side_pcoords[2*i+1];
      elem_pcoords[i*3+1] = side_pcoords[2*i+0];
      elem_pcoords[i*3+2] = side_pcoords[2*i+1];
    }
    break;
  case 2:
    for(int i = 0; i < npoints; i++) {
      elem_pcoords[i*3+0] = 0.0;
      elem_pcoords[i*3+1] = side_pcoords[2*i+1];
      elem_pcoords[i*3+2] = side_pcoords[2*i+0];
    }
    break;
  case 3:
    for(int i = 0; i < npoints; i++) {
      elem_pcoords[i*3+0] = side_pcoords[2*i+1];
      elem_pcoords[i*3+1] = side_pcoords[2*i+0];
      elem_pcoords[i*3+2] = 0.0;
    }
    break;
  default:
    throw std::runtime_error("Tet10::sideMap invalid ordinal");
  }
}

//--------------------------------------------------------------------------
//-------- tet10_fem_shape_fcn ---------------------------------------------
//--------------------------------------------------------------------------
void
Tet10FEM::tet10_fem_shape_fcn(
  const int  &numIp,
  const double *isoParCoord,
  SharedMemView<DoubleType**> shpfc)
{
  for ( int ip = 0; ip < numIp; ++ip ) {
    const int rowIpc = 3*ip;
    const DoubleType s1 = isoParCoord[rowIpc+0];
    const DoubleType s2 = isoParCoord[rowIpc+1];
    const DoubleType s3 = isoParCoord[rowIpc+2];

    shpfc(ip,0) = (1.0-2.0*s1-2.0*s2-2.0*s3)*(1.0-s1-s2-s3);
    shpfc(ip,1) = s1*(2.0*s1-1.0);
    shpfc(ip,2) = s2*(2.0*s2-1.0);
    shpfc(ip,3) = s3*(2.0*s3-1.0);
    shpfc(ip,4) = 4.0*s1*(1.0-s1-s2-s3);
    shpfc(ip,5) = 4.0*s1*s2;
    shpfc(ip,6) = 4.0*s2*(1.0-s1-s2-s3);
    shpfc(ip,7) = 4.0*s3*(1.0-s1-s2-s3);
    shpfc(ip,8) = 4.0*s1*s3;
    shpfc(ip,9) = 4.0*s2*s3;
  }
}

//--------------------------------------------------------------------------
//-------- tet10_fem_shape_fcn ---------------------------------------------
//--------------------------------------------------------------------------
void
Tet10FEM::tet10_fem_shape_fcn(
  const int  &numIp,
  const double *isoParCoord, 
  double *shpfc)
{
  const int npe = nodesPerElement_;
  for ( int ip = 0; ip < numIp; ++ip ) {
    const int rowIpc = 3*ip;
    const int rowSfc = npe*ip;
   
    const double s1 = isoParCoord[rowIpc+0];
    const double s2 = isoParCoord[rowIpc+1];
    const double s3 = isoParCoord[rowIpc+2];

    shpfc[rowSfc  ] = (1.0-2.0*s1-2.0*s2-2.0*s3)*(1.0-s1-s2-s3);
    shpfc[rowSfc+1] = s1*(2.0*s1-1.0);
    shpfc[rowSfc+2] = s2*(2.0*s2-1.0);
    shpfc[rowSfc+3] = s3*(2.0*s3-1.0);
    shpfc[rowSfc+4] = 4.0*s1*(1.0-s1-s2-s3);
    shpfc[rowSfc+5] = 4.0*s1*s2;
    shpfc[rowSfc+6] = 4.0*s2*(1.0-s1-s2-s3);
    shpfc[rowSfc+7] = 4.0*s3*(1.0-s1-s2-s3);
    shpfc[rowSfc+8] = 4.0*s1*s3;
    shpfc[rowSfc+9] = 4.0*s2*s3;
  }  
}

//--------------------------------------------------------------------------
//-------- tet10_derivative ------------------------------------------------
//--------------------------------------------------------------------------
void
Tet10FEM::tet10_derivative(
  const int  &npts,
  const double *intgLoc, 
  double *deriv)
{
  for ( int ip = 0; ip < npts; ++ip) {
    const int k = ip*3;
    const int p = 30*ip;
    
    const double s1 = intgLoc[k+0];
    const double s2 = intgLoc[k+1];
    const double s3 = intgLoc[k+2];

    // node 0
    deriv[0+3*0+p] = 4.0*(s1+s2+s3)-3.0;
    deriv[1+3*0+p] = 4.0*(s1+s2+s3)-3.0;
    deriv[2+3*0+p] = 4.0*(s1+s2+s3)-3.0;

    // node 1
    deriv[0+3*1+p] = 4.0*s1-1.0;
    deriv[1+3*1+p] = 0.0;
    deriv[2+3*1+p] = 0.0;

    // node 2
    deriv[0+3*2+p] = 0.0;
    deriv[1+3*2+p] = 4.0*s2-1.0;
    deriv[2+3*2+p] = 0.0;

    // node 3
    deriv[0+3*3+p] = 0.0;
    deriv[1+3*3+p] = 0.0;
    deriv[2+3*3+p] = 4.0*s3-1.0;

    // node 4
    deriv[0+3*4+p] = 4.0*(1.0-2.0*s1-s2-s3);
    deriv[1+3*4+p] =-4.0*s1;
    deriv[2+3*4+p] =-4.0*s1;

    // node 5
    deriv[0+3*5+p] = 4.0*s2;
    deriv[1+3*5+p] = 4.0*s1;
    deriv[2+3*5+p] = 0.0;

    // node 6
    deriv[0+3*6+p] =-4.0*s2;
    deriv[1+3*6+p] = 4.0*(1.0-s1-2.0*s2-s3);
    deriv[2+3*6+p] =-4.0*s2;

    // node 7
    deriv[0+3*7+p] =-4.0*s3;
    deriv[1+3*7+p] =-4.0*s3;
    deriv[2+3*7+p] = 4.0*(1.0-s1-s2-2.0*s3);

    // node 8
    deriv[0+3*8+p] = 4.0*s3;
    deriv[1+3*8+p] = 0.0;
    deriv[2+3*8+p] = 4.0*s1;

    // node 9
    deriv[0+3*9+p] = 0.0;
    deriv[1+3*9+p] = 4.0*s3;
    deriv[2+3*9+p] = 4.0*s2;
  }
}

//--------------------------------------------------------------------------
//-------- parametric_distance ---------------------------------------------
//--------------------------------------------------------------------------
double
Tet10FEM::parametric_distance(const double *x)
{
  const double X = x[0] - 1. / 4.;
  const double Y = x[1] - 1. / 4.;
  const double Z = x[2] - 1. / 4.;
  const double dist0 = -4 * X;
  const double dist1 = -4 * Y;
  const double dist2 = -4 * Z;
  const double dist3 = 4 * (X + Y + Z);
  const double dist = std::max(std::max(dist0, dist1), std::max(dist2, dist3));
  return dist;
}

}
}
