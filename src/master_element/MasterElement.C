/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <master_element/MasterElement.h>
#include <NaluEnv.h>
#include <FORTRAN_Proto.h>

#include <stk_topology/topology.hpp>

#include <iostream>

#include <cmath>
#include <limits>
#include <array>
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_RCP.hpp"

namespace sierra{
namespace nalu{


//--------------------------------------------------------------------------
//-------- factory for surface master elements -----------------------------
//--------------------------------------------------------------------------
MasterElement*
MasterElement::create_surface_master_element(stk::topology topo)
{
  switch ( topo.value() ) {

    case stk::topology::HEX_8:
      return new HexSCS();

    case stk::topology::HEX_27:
      return new Hex27SCS();

    case stk::topology::TET_4:
      return new TetSCS();

    case stk::topology::PYRAMID_5:
      return new PyrSCS();

    case stk::topology::WEDGE_6:
      return new WedSCS();

    case stk::topology::QUAD_4:
      return new Quad3DSCS();

    case stk::topology::QUAD_9:
      return new Quad93DSCS();

    case stk::topology::TRI_3:
      return new Tri3DSCS();

    case stk::topology::QUAD_4_2D:
      return new Quad2DSCS();

    case stk::topology::QUAD_9_2D:
      return new Quad92DSCS();

    case stk::topology::TRI_3_2D:
      return new Tri2DSCS();

    case stk::topology::LINE_2:
      return new Edge2DSCS();

    case stk::topology::LINE_3:
      return new Edge32DSCS();

    case stk::topology::SHELL_QUAD_4:
      NaluEnv::self().naluOutputP0() << "SHELL_QUAD_4 only supported for io surface transfer applications" << std::endl;
      return new Quad3DSCS();

    case stk::topology::SHELL_TRI_3:
      NaluEnv::self().naluOutputP0() << "SHELL_TRI_3 only supported for io surface transfer applications" << std::endl;
      return new Tri3DSCS();

    default:
      NaluEnv::self().naluOutputP0() << "sorry, we only support hex8, tet4, pyr5, wed6,"
                                        " quad2d, quad3d, tri2d, tri3d and edge2d surface elements" << std::endl;
      NaluEnv::self().naluOutputP0() << "your type is " << topo.value() << std::endl;
      break;

  }
  return nullptr;
}

//--------------------------------------------------------------------------
//-------- factory for volume master elements ------------------------------
//--------------------------------------------------------------------------
MasterElement*
MasterElement::create_volume_master_element(stk::topology topo)
{
  switch ( topo.value() ) {

    case stk::topology::HEX_8:
      return new HexSCV();

    case stk::topology::HEX_27:
      return new Hex27SCV();

    case stk::topology::TET_4:
      return new TetSCV();

    case stk::topology::PYRAMID_5:
      return new PyrSCV();

    case stk::topology::WEDGE_6:
      return  new WedSCV();

    case stk::topology::QUAD_4_2D:
      return new Quad2DSCV();

    case stk::topology::QUAD_9_2D:
      return new Quad92DSCV();

    case stk::topology::TRI_3_2D:
      return new Tri2DSCV();

    default:
      NaluEnv::self().naluOutputP0() << "sorry, we only support hex8, tet4, wed6, "
                                        " pyr5, quad4, and tri3 volume elements" << std::endl;
      NaluEnv::self().naluOutputP0() << "your type is " << topo.value() << std::endl;
      break;
  }
  return nullptr;
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MasterElement::MasterElement()
  : nDim_(0),
    nodesPerElement_(0),
    numIntPoints_(0),
    scaleToStandardIsoFac_(1.0)
{
  // nothing else
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
MasterElement::~MasterElement()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- isoparametric_mapping -------------------------------------------
//--------------------------------------------------------------------------
double
MasterElement::isoparametric_mapping( 
  const double b,
  const double a,
  const double xi) const
{
  return xi*(b-a)/2.0 +(a+b)/2.0;
}

//--------------------------------------------------------------------------
//-------- within_tolerance ------------------------------------------------
//--------------------------------------------------------------------------
bool 
MasterElement::within_tolerance( const double & val, const double & tol )
{
  return (std::abs(val)<tol);
}

//--------------------------------------------------------------------------
//-------- vector_norm_sq --------------------------------------------------
//--------------------------------------------------------------------------
double 
MasterElement::vector_norm_sq( const double * vect, int len )
{
  double norm_sq = 0.0;
  for (int i=0; i<len; i++) {
    norm_sq += vect[i]*vect[i];
  }
  return norm_sq;
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
HexSCV::HexSCV()
  : MasterElement()
{
  nDim_ = 3;
  nodesPerElement_ = 8;
  numIntPoints_ = 8;

  // define ip node mappings
  ipNodeMap_.resize(8);
  ipNodeMap_[0] = 0; ipNodeMap_[1] = 1; ipNodeMap_[2] = 2; ipNodeMap_[3] = 3;
  ipNodeMap_[4] = 4; ipNodeMap_[5] = 5; ipNodeMap_[6] = 6; ipNodeMap_[7] = 7;

  // standard integration location
  intgLoc_.resize(24);    
  intgLoc_[0]  = -0.25; intgLoc_[1]  = -0.25; intgLoc_[2]  = -0.25; 
  intgLoc_[3]  = +0.25; intgLoc_[4]  = -0.25; intgLoc_[5]  = -0.25; 
  intgLoc_[6]  = +0.25; intgLoc_[7]  = +0.25; intgLoc_[8]  = -0.25;
  intgLoc_[9]  = -0.25; intgLoc_[10] = +0.25; intgLoc_[11] = -0.25;
  intgLoc_[12] = -0.25; intgLoc_[13] = -0.25; intgLoc_[14] = +0.25;
  intgLoc_[15] = +0.25; intgLoc_[16] = -0.25; intgLoc_[17] = +0.25;
  intgLoc_[18] = +0.25; intgLoc_[19] = +0.25; intgLoc_[20] = +0.25;
  intgLoc_[21] = -0.25; intgLoc_[22] = +0.25; intgLoc_[23] = +0.25;

  // shifted integration location
  intgLocShift_.resize(24);    
  intgLocShift_[0]  = -0.5; intgLocShift_[1]  = -0.5; intgLocShift_[2]  = -0.5; 
  intgLocShift_[3]  = +0.5; intgLocShift_[4]  = -0.5; intgLocShift_[5]  = -0.5; 
  intgLocShift_[6]  = +0.5; intgLocShift_[7]  = +0.5; intgLocShift_[8]  = -0.5;
  intgLocShift_[9]  = -0.5; intgLocShift_[10] = +0.5; intgLocShift_[11] = -0.5;
  intgLocShift_[12] = -0.5; intgLocShift_[13] = -0.5; intgLocShift_[14] = +0.5;
  intgLocShift_[15] = +0.5; intgLocShift_[16] = -0.5; intgLocShift_[17] = +0.5;
  intgLocShift_[18] = +0.5; intgLocShift_[19] = +0.5; intgLocShift_[20] = +0.5;
  intgLocShift_[21] = -0.5; intgLocShift_[22] = +0.5; intgLocShift_[23] = +0.5;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
HexSCV::~HexSCV()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
HexSCV::ipNodeMap(
  int /*ordinal*/)
{
  // define scv->node mappings
  return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void HexSCV::determinant(
  const int nelem,
  const double *coords,
  double *volume,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(hex_scv_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords,
      volume, error, &lerr );
}

//--------------------------------------------------------------------------
//-------- grad_op ---------------------------------------------------------
//--------------------------------------------------------------------------
void HexSCV::grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(hex_derivative)
    ( &numIntPoints_,
      &intgLoc_[0], deriv );
  
  SIERRA_FORTRAN(hex_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );

  if ( lerr )
    std::cout << "sorry, negative HexSCV volume.." << std::endl;
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
HexSCV::shape_fcn(double *shpfc)
{
  SIERRA_FORTRAN(hex_shape_fcn)
    (&numIntPoints_,&intgLoc_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
HexSCV::shifted_shape_fcn(double *shpfc)
{
  SIERRA_FORTRAN(hex_shape_fcn)
    (&numIntPoints_,&intgLocShift_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
HexSCS::HexSCS()
  : MasterElement()
{
  nDim_ = 3;
  nodesPerElement_ = 8;
  numIntPoints_ = 12;
  scaleToStandardIsoFac_ = 2.0;

  // define L/R mappings
  lrscv_.resize(24);
  lrscv_[0]  = 0; lrscv_[1]  = 1;
  lrscv_[2]  = 1; lrscv_[3]  = 2;
  lrscv_[4]  = 2; lrscv_[5]  = 3;
  lrscv_[6]  = 0; lrscv_[7]  = 3;
  lrscv_[8]  = 4; lrscv_[9]  = 5;
  lrscv_[10] = 5; lrscv_[11] = 6;
  lrscv_[12] = 6; lrscv_[13] = 7;
  lrscv_[14] = 4; lrscv_[15] = 7;
  lrscv_[16] = 0; lrscv_[17] = 4;
  lrscv_[18] = 1; lrscv_[19] = 5;
  lrscv_[20] = 2; lrscv_[21] = 6;
  lrscv_[22] = 3; lrscv_[23] = 7;

  // define opposing node
  oppNode_.resize(24);
  // face 0
  oppNode_[0] = 3; oppNode_[1] = 2; oppNode_[2] = 6; oppNode_[3] = 7;
  // face 1
  oppNode_[4] = 0; oppNode_[5] = 3; oppNode_[6] = 7; oppNode_[7] = 4;
  // face 2
  oppNode_[8] = 1; oppNode_[9] = 0; oppNode_[10] = 4; oppNode_[11] = 5;
  // face 3
  oppNode_[12] = 1; oppNode_[13] = 5; oppNode_[14] = 6; oppNode_[15] = 2;
  // face 4
  oppNode_[16] = 4; oppNode_[17] = 7; oppNode_[18] = 6; oppNode_[19] = 5;
  // face 5
  oppNode_[20] = 0; oppNode_[21] = 1; oppNode_[22] = 2; oppNode_[23] = 3;

  // define opposing face
  oppFace_.resize(24);
  // face 0
  oppFace_[0]  = 3;  oppFace_[1] = 1;  oppFace_[2] = 5;   oppFace_[3] = 7;
  // face 1
  oppFace_[4]  = 0;  oppFace_[5] = 2;  oppFace_[6] = 6;   oppFace_[7] = 4;
  // face 2
  oppFace_[8]  = 1;  oppFace_[9] = 3;  oppFace_[10] = 7;  oppFace_[11] = 5;
  // face 3
  oppFace_[12] = 0; oppFace_[13] = 4;  oppFace_[14] = 6;  oppFace_[15] = 2;
  // face 4
  oppFace_[16] = 8; oppFace_[17] = 11; oppFace_[18] = 10; oppFace_[19] = 9;
  // face 5
  oppFace_[20] = 8; oppFace_[21] = 9;  oppFace_[22] = 10; oppFace_[23] = 11;

  // standard integration location
  intgLoc_.resize(36);    
  intgLoc_[0]  =  0.00; intgLoc_[1]  = -0.25; intgLoc_[2]  = -0.25; // surf 1    1->2
  intgLoc_[3]  =  0.25; intgLoc_[4]  =  0.00; intgLoc_[5]  = -0.25; // surf 2    2->3
  intgLoc_[6]  =  0.00; intgLoc_[7]  =  0.25; intgLoc_[8]  = -0.25; // surf 3    3->4
  intgLoc_[9]  = -0.25; intgLoc_[10] =  0.00; intgLoc_[11] = -0.25; // surf 4    1->4
  intgLoc_[12] =  0.00; intgLoc_[13] = -0.25; intgLoc_[14] =  0.25; // surf 5    5->6
  intgLoc_[15] =  0.25; intgLoc_[16] =  0.00; intgLoc_[17] =  0.25; // surf 6    6->7
  intgLoc_[18] =  0.00; intgLoc_[19] =  0.25; intgLoc_[20] =  0.25; // surf 7    7->8
  intgLoc_[21] = -0.25; intgLoc_[22] =  0.00; intgLoc_[23] =  0.25; // surf 8    5->8
  intgLoc_[24] = -0.25; intgLoc_[25] = -0.25; intgLoc_[26] =  0.00; // surf 9    1->5
  intgLoc_[27] =  0.25; intgLoc_[28] = -0.25; intgLoc_[29] =  0.00; // surf 10   2->6
  intgLoc_[30] =  0.25; intgLoc_[31] =  0.25; intgLoc_[32] =  0.00; // surf 11   3->7
  intgLoc_[33] = -0.25; intgLoc_[34] =  0.25; intgLoc_[35] =  0.00; // surf 12   4->8  

  // shifted
  intgLocShift_.resize(36);
  intgLocShift_[0]  =  0.00; intgLocShift_[1]  = -0.50; intgLocShift_[2]  = -0.50; // surf 1    1->2
  intgLocShift_[3]  =  0.50; intgLocShift_[4]  =  0.00; intgLocShift_[5]  = -0.50; // surf 2    2->3
  intgLocShift_[6]  =  0.00; intgLocShift_[7]  =  0.50; intgLocShift_[8]  = -0.50; // surf 3    3->4
  intgLocShift_[9]  = -0.50; intgLocShift_[10] =  0.00; intgLocShift_[11] = -0.50; // surf 4    1->4
  intgLocShift_[12] =  0.00; intgLocShift_[13] = -0.50; intgLocShift_[14] =  0.50; // surf 5    5->6
  intgLocShift_[15] =  0.50; intgLocShift_[16] =  0.00; intgLocShift_[17] =  0.50; // surf 6    6->7
  intgLocShift_[18] =  0.00; intgLocShift_[19] =  0.50; intgLocShift_[20] =  0.50; // surf 7    7->8
  intgLocShift_[21] = -0.50; intgLocShift_[22] =  0.00; intgLocShift_[23] =  0.50; // surf 8    5->8
  intgLocShift_[24] = -0.50; intgLocShift_[25] = -0.50; intgLocShift_[26] =  0.00; // surf 9    1->5
  intgLocShift_[27] =  0.50; intgLocShift_[28] = -0.50; intgLocShift_[29] =  0.00; // surf 10   2->6
  intgLocShift_[30] =  0.50; intgLocShift_[31] =  0.50; intgLocShift_[32] =  0.00; // surf 11   3->7
  intgLocShift_[33] = -0.50; intgLocShift_[34] =  0.50; intgLocShift_[35] =  0.00; // surf 12   4->8  

  // exposed face
  intgExpFace_.resize(72);
  // face 0; scs 0, 1, 2, 3
  intgExpFace_[0]  = -0.25; intgExpFace_[1]  = -0.50; intgExpFace_[2]  = -0.25;
  intgExpFace_[3]  =  0.25; intgExpFace_[4]  = -0.50; intgExpFace_[5]  = -0.25;
  intgExpFace_[6]  =  0.25; intgExpFace_[7]  = -0.50; intgExpFace_[8]  =  0.25;
  intgExpFace_[9]  = -0.25; intgExpFace_[10] = -0.50; intgExpFace_[11] =  0.25;
  // face 1; scs 0, 1, 2, 3
  intgExpFace_[12] =  0.50; intgExpFace_[13] = -0.25; intgExpFace_[14] = -0.25;
  intgExpFace_[15] =  0.50; intgExpFace_[16] =  0.25; intgExpFace_[17] = -0.25;
  intgExpFace_[18] =  0.50; intgExpFace_[19] =  0.25; intgExpFace_[20] =  0.25;
  intgExpFace_[21] =  0.50; intgExpFace_[22] = -0.25; intgExpFace_[23] =  0.25;
  // face 2; scs 0, 1, 2, 3
  intgExpFace_[24] =  0.25; intgExpFace_[25] =  0.50; intgExpFace_[26] = -0.25;
  intgExpFace_[27] = -0.25; intgExpFace_[28] =  0.50; intgExpFace_[29] = -0.25;
  intgExpFace_[30] = -0.25; intgExpFace_[31] =  0.50; intgExpFace_[32] =  0.25;
  intgExpFace_[33] =  0.25; intgExpFace_[34] =  0.50; intgExpFace_[35] =  0.25;
  // face 3; scs 0, 1, 2, 3
  intgExpFace_[36] = -0.50; intgExpFace_[37] = -0.25; intgExpFace_[38] = -0.25;
  intgExpFace_[39] = -0.50; intgExpFace_[40] = -0.25; intgExpFace_[41] =  0.25;
  intgExpFace_[42] = -0.50; intgExpFace_[43] =  0.25; intgExpFace_[44] =  0.25;
  intgExpFace_[45] = -0.50; intgExpFace_[46] =  0.25; intgExpFace_[47] = -0.25;
  // face 4; scs 0, 1, 2, 3
  intgExpFace_[48] = -0.25; intgExpFace_[49] = -0.25; intgExpFace_[50] = -0.50;
  intgExpFace_[51] = -0.25; intgExpFace_[52] =  0.25; intgExpFace_[53] = -0.50;
  intgExpFace_[54] =  0.25; intgExpFace_[55] =  0.25; intgExpFace_[56] = -0.50;
  intgExpFace_[57] =  0.25; intgExpFace_[58] = -0.25; intgExpFace_[59] = -0.50;
  // face 5; scs 0, 1, 2, 3
  intgExpFace_[60] = -0.25; intgExpFace_[61] = -0.25; intgExpFace_[62] =  0.50;
  intgExpFace_[63] =  0.25; intgExpFace_[64] = -0.25; intgExpFace_[65] =  0.50;
  intgExpFace_[66] =  0.25; intgExpFace_[67] =  0.25; intgExpFace_[68] =  0.50;
  intgExpFace_[69] = -0.25; intgExpFace_[70] =  0.25; intgExpFace_[71] =  0.50;  

  // boundary integration point ip node mapping (ip on an ordinal to local node number)
  ipNodeMap_.resize(24); // 4 ips * 6 faces
  // face 0;
  ipNodeMap_[0] = 0;  ipNodeMap_[1] = 1;  ipNodeMap_[2] = 5;  ipNodeMap_[3] = 4; 
  // face 1; 
  ipNodeMap_[4] = 1;  ipNodeMap_[5] = 2;  ipNodeMap_[6] = 6;  ipNodeMap_[7] = 5; 
  // face 2;
  ipNodeMap_[8] = 2;  ipNodeMap_[9] = 3;  ipNodeMap_[10] = 7; ipNodeMap_[11] = 6; 
  // face 3;
  ipNodeMap_[12] = 0; ipNodeMap_[13] = 4; ipNodeMap_[14] = 7; ipNodeMap_[15] = 3; 
  // face 4;
  ipNodeMap_[16] = 0; ipNodeMap_[17] = 3; ipNodeMap_[18] = 2; ipNodeMap_[19] = 1; 
  // face 5;
  ipNodeMap_[20] = 4; ipNodeMap_[21] = 5; ipNodeMap_[22] = 6; ipNodeMap_[23] = 7; 

  // nodes for collocation calculations
  nodeLoc_.resize(24);
  // node 0
  nodeLoc_[0] = -0.5; nodeLoc_[1] = -0.5; nodeLoc_[2] = -0.5;
  // node 1
  nodeLoc_[3] =  0.5; nodeLoc_[4] = -0.5; nodeLoc_[5] = -0.5;
  // node 2
  nodeLoc_[6] =  0.5; nodeLoc_[7] =  0.5; nodeLoc_[8] = -0.5;
  // node 3
  nodeLoc_[9] = -0.5; nodeLoc_[10] =  0.5; nodeLoc_[11] = -0.5;
  // node 4
  nodeLoc_[12] = -0.5; nodeLoc_[13] = -0.5; nodeLoc_[14] =  0.5;
  // node 5
  nodeLoc_[15] =  0.5; nodeLoc_[16] = -0.5; nodeLoc_[17] =  0.5;
  // node 6
  nodeLoc_[18] =  0.5; nodeLoc_[19] =  0.5; nodeLoc_[20] =  0.5;
  // node 7
  nodeLoc_[21] = -0.5; nodeLoc_[22] =  0.5; nodeLoc_[23] =  0.5;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
HexSCS::~HexSCS()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
HexSCS::ipNodeMap(
  int ordinal)
{
  // define ip->node mappings for each face (ordinal); 
  return &ipNodeMap_[ordinal*4];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void HexSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  SIERRA_FORTRAN(hex_scs_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords, areav );

  // all is always well; no error checking
  *error = 0;
}

//--------------------------------------------------------------------------
//-------- grad_op ---------------------------------------------------------
//--------------------------------------------------------------------------
void HexSCS::grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(hex_derivative)
    ( &numIntPoints_,
      &intgLoc_[0], deriv );
  
  SIERRA_FORTRAN(hex_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );

  if ( lerr )
    std::cout << "sorry, negative HexSCS volume.." << std::endl;
 }

//--------------------------------------------------------------------------
//-------- shifted_grad_op -------------------------------------------------
//--------------------------------------------------------------------------
void HexSCS::shifted_grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(hex_derivative)
    ( &numIntPoints_,
      &intgLocShift_[0], deriv );
  
  SIERRA_FORTRAN(hex_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );

  if ( lerr )
    std::cout << "sorry, negative HexSCS volume.." << std::endl;
}

//--------------------------------------------------------------------------
//-------- face_grad_op ----------------------------------------------------
//--------------------------------------------------------------------------
void HexSCS::face_grad_op(
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
      SIERRA_FORTRAN(hex_derivative)
        ( &nface,
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
//-------- guij ------------------------------------------------------------
//--------------------------------------------------------------------------
void HexSCS::gij(
  const double *coords,
  double *gupperij,
  double *glowerij,
  double *deriv)
{
  SIERRA_FORTRAN(threed_gij)
    ( &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gupperij, glowerij);
}

//--------------------------------------------------------------------------
//-------- adjacentNodes ---------------------------------------------------
//--------------------------------------------------------------------------
const int *
HexSCS::adjacentNodes()
{
  // define L/R mappings
  return &lrscv_[0];
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
HexSCS::shape_fcn(double *shpfc)
{
  SIERRA_FORTRAN(hex_shape_fcn)
    (&numIntPoints_,&intgLoc_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
HexSCS::shifted_shape_fcn(double *shpfc)
{
  SIERRA_FORTRAN(hex_shape_fcn)
    (&numIntPoints_,&intgLocShift_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- opposingNodes --------------------------------------------------
//--------------------------------------------------------------------------
int
HexSCS::opposingNodes(
  const int ordinal,
  const int node)
{
  return oppNode_[ordinal*4+node];
}

//--------------------------------------------------------------------------
//-------- opposingFace --------------------------------------------------
//--------------------------------------------------------------------------
int
HexSCS::opposingFace(
  const int ordinal,
  const int node)
{
  return oppFace_[ordinal*4+node];
}

//--------------------------------------------------------------------------
//-------- isInElement -----------------------------------------------------
//--------------------------------------------------------------------------
double
HexSCS::isInElement(
    const double * elem_nodal_coor,     // (8,3)
    const double * point_coor,          // (3)
	  double * par_coor ) 
{
  const int maxNonlinearIter = 20;
  const double isInElemConverged = 1.0e-16;
  // Translate element so that (x,y,z) coordinates of the first node are (0,0,0)

  double x[] = {0.,
	      0.125*(elem_nodal_coor[1] - elem_nodal_coor[0]),
	      0.125*(elem_nodal_coor[2] - elem_nodal_coor[0]),
	      0.125*(elem_nodal_coor[3] - elem_nodal_coor[0]),
	      0.125*(elem_nodal_coor[4] - elem_nodal_coor[0]),
	      0.125*(elem_nodal_coor[5] - elem_nodal_coor[0]),
	      0.125*(elem_nodal_coor[6] - elem_nodal_coor[0]),
	      0.125*(elem_nodal_coor[7] - elem_nodal_coor[0]) };
  double y[] = {0.,
	      0.125*(elem_nodal_coor[9 ] - elem_nodal_coor[8]),
	      0.125*(elem_nodal_coor[10] - elem_nodal_coor[8]),
	      0.125*(elem_nodal_coor[11] - elem_nodal_coor[8]),
	      0.125*(elem_nodal_coor[12] - elem_nodal_coor[8]),
	      0.125*(elem_nodal_coor[13] - elem_nodal_coor[8]),
	      0.125*(elem_nodal_coor[14] - elem_nodal_coor[8]),
	      0.125*(elem_nodal_coor[15] - elem_nodal_coor[8]) };
  double z[] = {0.,
	      0.125*(elem_nodal_coor[17] - elem_nodal_coor[16]),
	      0.125*(elem_nodal_coor[18] - elem_nodal_coor[16]),
	      0.125*(elem_nodal_coor[19] - elem_nodal_coor[16]),
	      0.125*(elem_nodal_coor[20] - elem_nodal_coor[16]),
	      0.125*(elem_nodal_coor[21] - elem_nodal_coor[16]),
	      0.125*(elem_nodal_coor[22] - elem_nodal_coor[16]),
	      0.125*(elem_nodal_coor[23] - elem_nodal_coor[16]) };

  // (xp,yp,zp) is the point at which we're searching for (xi,eta,zeta)
  // (must translate this also)

  double xp = point_coor[0] - elem_nodal_coor[0];
  double yp = point_coor[1] - elem_nodal_coor[8];
  double zp = point_coor[2] - elem_nodal_coor[16];

  // Newton-Raphson iteration for (xi,eta,zeta)
  double j[9];
  double f[3];
  double shapefct[8];
  double xinew = 0.5;     // initial guess
  double etanew = 0.5;
  double zetanew = 0.5;
  double xicur = 0.5;
  double etacur = 0.5;
  double zetacur = 0.5;
  double xidiff[] = { 1.0, 1.0, 1.0 };
  int i = 0;

  do
  {
    j[0]=
      -(1.0-etacur)*(1.0-zetacur)*x[1]
      -(1.0+etacur)*(1.0-zetacur)*x[2]
      +(1.0+etacur)*(1.0-zetacur)*x[3]
      +(1.0-etacur)*(1.0+zetacur)*x[4]
      -(1.0-etacur)*(1.0+zetacur)*x[5]
      -(1.0+etacur)*(1.0+zetacur)*x[6]
      +(1.0+etacur)*(1.0+zetacur)*x[7];

    j[1]=
       (1.0+xicur)*(1.0-zetacur)*x[1]
      -(1.0+xicur)*(1.0-zetacur)*x[2]
      -(1.0-xicur)*(1.0-zetacur)*x[3]
      +(1.0-xicur)*(1.0+zetacur)*x[4]
      +(1.0+xicur)*(1.0+zetacur)*x[5]
      -(1.0+xicur)*(1.0+zetacur)*x[6]
      -(1.0-xicur)*(1.0+zetacur)*x[7];

    j[2]=
       (1.0-etacur)*(1.0+xicur)*x[1]
      +(1.0+etacur)*(1.0+xicur)*x[2]
      +(1.0+etacur)*(1.0-xicur)*x[3]
      -(1.0-etacur)*(1.0-xicur)*x[4]
      -(1.0-etacur)*(1.0+xicur)*x[5]
      -(1.0+etacur)*(1.0+xicur)*x[6]
      -(1.0+etacur)*(1.0-xicur)*x[7];

    j[3]=
      -(1.0-etacur)*(1.0-zetacur)*y[1]
      -(1.0+etacur)*(1.0-zetacur)*y[2]
      +(1.0+etacur)*(1.0-zetacur)*y[3]
      +(1.0-etacur)*(1.0+zetacur)*y[4]
      -(1.0-etacur)*(1.0+zetacur)*y[5]
      -(1.0+etacur)*(1.0+zetacur)*y[6]
      +(1.0+etacur)*(1.0+zetacur)*y[7];

    j[4]=
       (1.0+xicur)*(1.0-zetacur)*y[1]
      -(1.0+xicur)*(1.0-zetacur)*y[2]
      -(1.0-xicur)*(1.0-zetacur)*y[3]
      +(1.0-xicur)*(1.0+zetacur)*y[4]
      +(1.0+xicur)*(1.0+zetacur)*y[5]
      -(1.0+xicur)*(1.0+zetacur)*y[6]
      -(1.0-xicur)*(1.0+zetacur)*y[7];

    j[5]=
       (1.0-etacur)*(1.0+xicur)*y[1]
      +(1.0+etacur)*(1.0+xicur)*y[2]
      +(1.0+etacur)*(1.0-xicur)*y[3]
      -(1.0-etacur)*(1.0-xicur)*y[4]
      -(1.0-etacur)*(1.0+xicur)*y[5]
      -(1.0+etacur)*(1.0+xicur)*y[6]
      -(1.0+etacur)*(1.0-xicur)*y[7];

    j[6]=
      -(1.0-etacur)*(1.0-zetacur)*z[1]
      -(1.0+etacur)*(1.0-zetacur)*z[2]
      +(1.0+etacur)*(1.0-zetacur)*z[3]
      +(1.0-etacur)*(1.0+zetacur)*z[4]
      -(1.0-etacur)*(1.0+zetacur)*z[5]
      -(1.0+etacur)*(1.0+zetacur)*z[6]
      +(1.0+etacur)*(1.0+zetacur)*z[7];

    j[7]=
       (1.0+xicur)*(1.0-zetacur)*z[1]
      -(1.0+xicur)*(1.0-zetacur)*z[2]
      -(1.0-xicur)*(1.0-zetacur)*z[3]
      +(1.0-xicur)*(1.0+zetacur)*z[4]
      +(1.0+xicur)*(1.0+zetacur)*z[5]
      -(1.0+xicur)*(1.0+zetacur)*z[6]
      -(1.0-xicur)*(1.0+zetacur)*z[7];

    j[8]=
       (1.0-etacur)*(1.0+xicur)*z[1]
      +(1.0+etacur)*(1.0+xicur)*z[2]
      +(1.0+etacur)*(1.0-xicur)*z[3]
      -(1.0-etacur)*(1.0-xicur)*z[4]
      -(1.0-etacur)*(1.0+xicur)*z[5]
      -(1.0+etacur)*(1.0+xicur)*z[6]
      -(1.0+etacur)*(1.0-xicur)*z[7];

    double jdet=-(j[2]*j[4]*j[6])+j[1]*j[5]*j[6]+j[2]*j[3]*j[7]-
      j[0]*j[5]*j[7]-j[1]*j[3]*j[8]+j[0]*j[4]*j[8];

    if (!jdet) {
      i = maxNonlinearIter;
      break;
    }
    shapefct[0]=(1.0-etacur)*(1.0-xicur)*(1.0-zetacur);

    shapefct[1]=(1.0-etacur)*(1.0+xicur)*(1.0-zetacur);

    shapefct[2]=(1.0+etacur)*(1.0+xicur)*(1.0-zetacur);

    shapefct[3]=(1.0+etacur)*(1.0-xicur)*(1.0-zetacur);

    shapefct[4]=(1.0-etacur)*(1.0-xicur)*(1.0+zetacur);

    shapefct[5]=(1.0-etacur)*(1.0+xicur)*(1.0+zetacur);

    shapefct[6]=(1.0+etacur)*(1.0+xicur)*(1.0+zetacur);

    shapefct[7]=(1.0+etacur)*(1.0-xicur)*(1.0+zetacur);

    f[0]=xp-shapefct[1]*x[1]-shapefct[2]*x[2]-shapefct[3]*x[3]-shapefct[4]*x[4]-\
      shapefct[5]*x[5]-shapefct[6]*x[6]-shapefct[7]*x[7];

    f[1]=yp-shapefct[1]*y[1]-shapefct[2]*y[2]-shapefct[3]*y[3]-shapefct[4]*y[4]-\
      shapefct[5]*y[5]-shapefct[6]*y[6]-shapefct[7]*y[7];

    f[2]=zp-shapefct[1]*z[1]-shapefct[2]*z[2]-shapefct[3]*z[3]-shapefct[4]*z[4]-\
      shapefct[5]*z[5]-shapefct[6]*z[6]-shapefct[7]*z[7];

    xinew = (jdet*xicur+f[2]*(j[2]*j[4]-j[1]*j[5])-f[1]*j[2]*j[7]+f[0]*j[5]*j[7]+
	     f[1]*j[1]*j[8]-f[0]*j[4]*j[8])/jdet;

    etanew = (etacur*jdet+f[2]*(-(j[2]*j[3])+j[0]*j[5])+f[1]*j[2]*j[6]-f[0]*j[5]*j[6]-
	      f[1]*j[0]*j[8]+f[0]*j[3]*j[8])/jdet;

    zetanew = (jdet*zetacur+f[2]*(j[1]*j[3]-j[0]*j[4])-f[1]*j[1]*j[6]+
	       f[0]*j[4]*j[6]+f[1]*j[0]*j[7]-f[0]*j[3]*j[7])/jdet;

    xidiff[0] = xinew - xicur;
    xidiff[1] = etanew - etacur;
    xidiff[2] = zetanew - zetacur;
    xicur = xinew;
    etacur = etanew;
    zetacur = zetanew;

  }
  while ( !within_tolerance( vector_norm_sq(xidiff,3), isInElemConverged) && ++i < maxNonlinearIter);

  par_coor[0] = par_coor[1] = par_coor[2] = std::numeric_limits<double>::max();
  double dist = std::numeric_limits<double>::max();

  if (i <maxNonlinearIter) {
    par_coor[0] = xinew;
    par_coor[1] = etanew;
    par_coor[2] = zetanew;

    std::vector<double> xtmp(3);
    xtmp[0] = par_coor[0];
    xtmp[1] = par_coor[1];
    xtmp[2] = par_coor[2];
    dist = parametric_distance(xtmp);
  }
  return dist;
}

//--------------------------------------------------------------------------
//-------- interpolatePoint ------------------------------------------------
//--------------------------------------------------------------------------
void
HexSCS::interpolatePoint(
    const int  & ncomp_field,
    const double * par_coord,           // (3)
    const double * field,               // (8,ncomp_field)
	  double * result ) // (ncomp_field)
{
  // 'field' is a flat array of dimension (8,ncomp_field) (Fortran ordering);
  double xi   = par_coord[0];
  double eta  = par_coord[1];
  double zeta = par_coord[2];

  for ( int i = 0; i < ncomp_field; i++ )
  {
    // Base 'field array' index for ith component
    int b = 8*i;

    result[i] = 0.125000000000000*(1.00000000000000-eta )*
				  (1.00000000000000-xi  )*
				  (1.00000000000000-zeta)*field[b+0]+
				  0.125000000000000*(1.00000000000000-eta )*
				  (1.00000000000000+xi  )*
				  (1.00000000000000-zeta)*field[b+1]+
				  0.125000000000000*(1.00000000000000+eta )*
				  (1.00000000000000+xi  )*
				  (1.00000000000000-zeta)*field[b+2]+
				  0.125000000000000*(1.00000000000000+eta )*
				  (1.00000000000000-xi  )*
				  (1.00000000000000-zeta)*field[b+3]+
				  0.125000000000000*(1.00000000000000-eta )*
				  (1.00000000000000-xi  )*
				  (1.00000000000000+zeta)*field[b+4]+
				  0.125000000000000*(1.00000000000000-eta )*
				  (1.00000000000000+xi  )*
				  (1.00000000000000+zeta)*field[b+5]+
				  0.125000000000000*(1.00000000000000+eta )*
				  (1.00000000000000+xi  )*
				  (1.00000000000000+zeta)*field[b+6]+
				  0.125000000000000*(1.00000000000000+eta )*
				  (1.00000000000000-xi  )*
				  (1.00000000000000+zeta)*field[b+7];
  }
}

//--------------------------------------------------------------------------
//-------- general_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
HexSCS::general_shape_fcn(
  const int numIp,
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
    shpfc[rowSfc+1] = 0.125*(1.0-s1)*(1.0+s2)*(1.0-s3);
    shpfc[rowSfc+2] = 0.125*(1.0+s1)*(1.0+s2)*(1.0-s3);
    shpfc[rowSfc+3] = 0.125*(1.0+s1)*(1.0-s2)*(1.0-s3);
    shpfc[rowSfc+4] = 0.125*(1.0-s1)*(1.0-s2)*(1.0+s3);
    shpfc[rowSfc+5] = 0.125*(1.0-s1)*(1.0+s2)*(1.0+s3);
    shpfc[rowSfc+6] = 0.125*(1.0+s1)*(1.0+s2)*(1.0+s3);
    shpfc[rowSfc+7] = 0.125*(1.0+s1)*(1.0-s2)*(1.0+s3);
 
  }
}

//--------------------------------------------------------------------------
//-------- general_face_grad_op --------------------------------------------
//--------------------------------------------------------------------------
void 
HexSCS::general_face_grad_op(
  const int face_ordinal,
  const double *isoParCoord,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  int lerr = 0;
  const int nface = 1;

  double dpsi[24];
  double grad[24];

  SIERRA_FORTRAN(hex_derivative)
    ( &nface, &isoParCoord[0], dpsi );
      
  SIERRA_FORTRAN(hex_gradient_operator)
    ( &nface,
      &nodesPerElement_,
      &nface,
      dpsi,
      &coords[0], grad, &det_j[0], error, &lerr );
  
  if ( lerr )
    std::cout << "HexSCS::general_face_grad_op: issue.." << std::endl;
  
  for ( int j=0; j<24; j++) {
    gradop[j] = grad[j];
  }
}

//--------------------------------------------------------------------------
//-------- sidePcoords_to_elemPcoords --------------------------------------
//--------------------------------------------------------------------------
void 
HexSCS::sidePcoords_to_elemPcoords(
  const int & side_ordinal,
  const int & npoints,
  const double *side_pcoords,
  double *elem_pcoords)
{
  switch (side_ordinal) {
  case 0:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = 0.5*side_pcoords[2*i+0];
      elem_pcoords[i*3+1] = -0.5;
      elem_pcoords[i*3+2] = 0.5*side_pcoords[2*i+1];
    }
    break;
  case 1:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = 0.5;
      elem_pcoords[i*3+1] = 0.5*side_pcoords[2*i+0];
      elem_pcoords[i*3+2] = 0.5*side_pcoords[2*i+1];
    }
    break;
  case 2:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = -0.5*side_pcoords[2*i+0];
      elem_pcoords[i*3+1] = 0.5;
      elem_pcoords[i*3+2] = 0.5*side_pcoords[2*i+1];
    }
    break;
  case 3:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = -0.5;
      elem_pcoords[i*3+1] = 0.5*side_pcoords[2*i+1];
      elem_pcoords[i*3+2] = 0.5*side_pcoords[2*i+0];
    }
    break;
  case 4:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = 0.5*side_pcoords[2*i+1];
      elem_pcoords[i*3+1] = 0.5*side_pcoords[2*i+0];
      elem_pcoords[i*3+2] = -0.5;
    }
    break;
  case 5:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = 0.5*side_pcoords[2*i+0];
      elem_pcoords[i*3+1] = 0.5*side_pcoords[2*i+1];
      elem_pcoords[i*3+2] = 0.5;
    }
    break;
  default:
    throw std::runtime_error("HexSCS::sideMap invalid ordinal");
  }
}

//--------------------------------------------------------------------------
//-------- parametric_distance ---------------------------------------------
//--------------------------------------------------------------------------
double HexSCS::parametric_distance(const std::vector<double> &x)
{
  std::vector<double> y(3);
  for (int i=0; i<3; ++i) {
    y[i] = std::fabs(x[i]);
  }

  double d = 0;
  for (int i=0; i<3; ++i) {
    if (d < y[i]) {
      d = y[i];
    }
  }
  return d;
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
HexahedralP2Element::HexahedralP2Element()
  : MasterElement(),
    scsDist_(std::sqrt(3.0)/3.0),
    useGLLGLL_(false),
    nodes1D_(3),
    numQuad_(2)
{
  nDim_ = 3;
  nodesPerElement_ = nodes1D_ * nodes1D_ * nodes1D_;

  // map the standard stk node numbering to a tensor-product style node numbering (i.e. node (m,l,k) -> m+npe*l+npe^2*k)
  stkNodeMap_ = {
                   0,  8,  1, // bottom front edge
                  11, 21,  9, // bottom mid-front edge
                   3, 10,  2, // bottom back edge
                  12, 25, 13, // mid-top front edge
                  23, 20, 24, // mid-top mid-front edge
                  15, 26, 14, // mid-top back edge
                   4, 16,  5, // top front edge
                  19, 22, 17, // top mid-front edge
                   7, 18,  6  // top back edge
                };

  // a padded list of the scs locations
  scsEndLoc_ = { -1.0, -scsDist_, scsDist_, 1.0 };
}

//--------------------------------------------------------------------------
//-------- set_quadrature_rule ---------------------------------------------
//--------------------------------------------------------------------------
void
HexahedralP2Element::set_quadrature_rule()
{
  gaussAbscissaeShift_ = {-1.0,-1.0,0.0,0.0,+1.0,+1.0};

  if (useGLLGLL_ && numQuad_ < 3) {
    throw std::runtime_error("Need at least 3 points for gllgll quadrature");
  }

  switch (numQuad_) {
    case 1:
    {
      gaussAbscissae_ = { 0.0 };
      gaussWeight_ = { 1.0 };
      break;
    }
    case 2:
    {
      gaussAbscissae_ = { -std::sqrt(3.0)/3.0, std::sqrt(3.0)/3.0 };
      gaussWeight_ = { 0.5, 0.5 };
      break;
    }
    case 3:
    {
      if (!useGLLGLL_) {
        gaussAbscissae_ = { -std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0) };
        gaussWeight_ = { 5.0/18.0, 4.0/9.0,  5.0/18.0 };
      }
      else {
        gaussAbscissae1D_ = {-1.0, 0.0, +1.0};
        GLLGLL_quadrature_weights();
        gaussAbscissae_ = {
            gaussAbscissae1D_[0], gaussAbscissae1D_[1], gaussAbscissae1D_[2],
            gaussAbscissae1D_[0], gaussAbscissae1D_[1], gaussAbscissae1D_[2],
            gaussAbscissae1D_[0], gaussAbscissae1D_[1], gaussAbscissae1D_[2]
        };
      }
      break;
    }
    case 4:
    {
      if (!useGLLGLL_) {
        gaussAbscissae_ = {
            -std::sqrt(3.0/7.0+2.0/7.0*std::sqrt(6.0/5.0)),
            -std::sqrt(3.0/7.0-2.0/7.0*std::sqrt(6.0/5.0)),
            +std::sqrt(3.0/7.0-2.0/7.0*std::sqrt(6.0/5.0)),
            +std::sqrt(3.0/7.0+2.0/7.0*std::sqrt(6.0/5.0))
        };

        gaussWeight_ = {
            (18.0-std::sqrt(30.0))/72.0,
            (18.0+std::sqrt(30.0))/72.0,
            (18.0+std::sqrt(30.0))/72.0,
            (18.0-std::sqrt(30.0))/72.0
        };
      }
      else {
        gaussAbscissae1D_ = { -1.0, -std::sqrt(5.0)/5.0, std::sqrt(5.0)/5.0, +1.0};
        GLLGLL_quadrature_weights();

        gaussAbscissae_ = {
            gaussAbscissae1D_[0], gaussAbscissae1D_[1], gaussAbscissae1D_[2], gaussAbscissae1D_[3],
            gaussAbscissae1D_[0], gaussAbscissae1D_[1], gaussAbscissae1D_[2], gaussAbscissae1D_[3],
            gaussAbscissae1D_[0], gaussAbscissae1D_[1], gaussAbscissae1D_[2], gaussAbscissae1D_[3]
        };
      }
      break;
    }
    default:
      throw std::runtime_error("Quadrature rule not implemented");
  }
}

//--------------------------------------------------------------------------
//-------- GLL_quadrature_weights ----------------------------------------
//--------------------------------------------------------------------------
void
HexahedralP2Element::GLLGLL_quadrature_weights()
{
  // for fixed abscissae, use a moment-matching equation to compute weights

  const int nrows = gaussAbscissae1D_.size();
  const int nrhs = nodes1D_;

  Teuchos::SerialDenseMatrix<int, double> weightLHS(nrows, nrows);
  for (int j = 0; j < nrows; ++j) {
    for (int i = 0; i < nrows; ++i) {
      weightLHS(i, j) = std::pow(gaussAbscissae1D_[j], i);
    }
  }

  // each node has a separate RHS
  Teuchos::SerialDenseMatrix<int, double> weightRHS(nrows, nrhs);
  for (int j = 0; j < nrhs; ++j) {
    for (int i = 0; i < nrows; ++i) {
      weightRHS(i, j) = (std::pow(scsEndLoc_[j + 1], i + 1)
                       - std::pow(scsEndLoc_[j], i + 1)) / (i + 1.0);
    }
  }

  Teuchos::SerialDenseSolver<int, double> solver;
  Teuchos::SerialDenseMatrix<int, double> quadratureWeights(nrows, nrhs);
  solver.setMatrix(Teuchos::rcp(&weightLHS, false));
  solver.setVectors(
    Teuchos::rcp(&quadratureWeights, false),
    Teuchos::rcp(&weightRHS, false)
  );
  solver.solve();

  gaussWeight_.resize(nrows * nrhs);
  for (int j = 0; j < nrows; ++j) {
    for (int i = 0; i < nrhs; ++i) {
      gaussWeight_[i + j * nrhs] = quadratureWeights(i, j); //transpose
    }
  }
}

//--------------------------------------------------------------------------
//-------- tensor_product_node_map -----------------------------------------
//--------------------------------------------------------------------------
int
HexahedralP2Element::tensor_product_node_map(int i, int j, int k) const
{
   return stkNodeMap_[i+j*nodes1D_+k*nodes1D_*nodes1D_];
}

//--------------------------------------------------------------------------
//-------- gauss_point_location --------------------------------------------
//--------------------------------------------------------------------------
double
HexahedralP2Element::gauss_point_location(
  int nodeOrdinal,
  int gaussPointOrdinal) const
{
  double location1D;
  if (!useGLLGLL_) {
    location1D = isoparametric_mapping( scsEndLoc_[nodeOrdinal+1],
                                        scsEndLoc_[nodeOrdinal],
                                        gaussAbscissae_[gaussPointOrdinal] );
  }
  else {
    location1D = gaussAbscissae_[nodeOrdinal*numQuad_ + gaussPointOrdinal];
  }
   return location1D;
}

//--------------------------------------------------------------------------
//-------- gauss_point_location --------------------------------------------
//--------------------------------------------------------------------------
double
HexahedralP2Element::shifted_gauss_point_location(
  int nodeOrdinal,
  int gaussPointOrdinal) const
{
  return gaussAbscissaeShift_[nodeOrdinal*numQuad_ + gaussPointOrdinal];
}

//--------------------------------------------------------------------------
//-------- tensor_product_weight -------------------------------------------
//--------------------------------------------------------------------------

double
HexahedralP2Element::tensor_product_weight(
  int s1Node, int s2Node, int s3Node,
  int s1Ip, int s2Ip, int s3Ip) const
{
  // volume integration
  double weight;
  if (!useGLLGLL_) {
    const double Ls1 = scsEndLoc_[s1Node+1]-scsEndLoc_[s1Node];
    const double Ls2 = scsEndLoc_[s2Node+1]-scsEndLoc_[s2Node];
    const double Ls3 = scsEndLoc_[s3Node+1]-scsEndLoc_[s3Node];
    const double isoparametricArea = Ls1 * Ls2 * Ls3;

    weight = isoparametricArea
           * gaussWeight_[s1Ip]
           * gaussWeight_[s2Ip]
           * gaussWeight_[s3Ip];
   }
   else {
     weight = gaussWeight_[s1Node*numQuad_+s1Ip]
            * gaussWeight_[s2Node*numQuad_+s2Ip]
            * gaussWeight_[s3Node*numQuad_+s3Ip];

   }
   return weight;
}

//--------------------------------------------------------------------------
//-------- tensor_product_weight -------------------------------------------
//--------------------------------------------------------------------------
double
HexahedralP2Element::tensor_product_weight(
  int s1Node, int s2Node,
  int s1Ip, int s2Ip) const
{
  // surface integration
  double weight;
  if (!useGLLGLL_) {
    const double Ls1 = scsEndLoc_[s1Node+1]-scsEndLoc_[s1Node];
    const double Ls2 = scsEndLoc_[s2Node+1]-scsEndLoc_[s2Node];
    const double isoparametricArea = Ls1 * Ls2;

    weight = isoparametricArea * gaussWeight_[s1Ip] * gaussWeight_[s2Ip];
   }
   else {
     weight = gaussWeight_[s1Node*numQuad_+s1Ip]
            * gaussWeight_[s2Node*numQuad_+s2Ip];

   }
   return weight;
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
HexahedralP2Element::shape_fcn(double* shpfc)
{
  for (int ip = 0; ip < numIntPoints_ * nodesPerElement_; ++ip) {
    shpfc[ip] = shapeFunctions_[ip];
  }
}
//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
HexahedralP2Element::shifted_shape_fcn(double* shpfc)
{
  for (int ip = 0; ip < numIntPoints_ * nodesPerElement_; ++ip) {
    shpfc[ip] = shapeFunctionsShift_[ip];
  }
}


//--------------------------------------------------------------------------
//-------- eval_shape_functions_at_ips -------------------------------------
//--------------------------------------------------------------------------
void
HexahedralP2Element::eval_shape_functions_at_ips()
{
  shapeFunctions_.resize(numIntPoints_*nodesPerElement_);
  hex27_shape_fcn(numIntPoints_, intgLoc_.data(), shapeFunctions_.data());
}

//--------------------------------------------------------------------------
//-------- eval_shape_derivs_at_ips ----------------------------------------
//--------------------------------------------------------------------------
void
HexahedralP2Element::eval_shape_derivs_at_ips()
{
  shapeDerivs_.resize(numIntPoints_*nodesPerElement_*nDim_);
  hex27_shape_deriv(numIntPoints_, intgLoc_.data(), shapeDerivs_.data());
}

//--------------------------------------------------------------------------
//-------- eval_shape_functions_at_shifted_ips -----------------------------
//--------------------------------------------------------------------------
void
HexahedralP2Element::eval_shape_functions_at_shifted_ips()
{
  shapeFunctionsShift_.resize(numIntPoints_*nodesPerElement_);
  hex27_shape_fcn(numIntPoints_, intgLocShift_.data(), shapeFunctionsShift_.data());
}

//--------------------------------------------------------------------------
//-------- eval_shape_derivs_at_ips ----------------------------------------
//--------------------------------------------------------------------------
void
HexahedralP2Element::eval_shape_derivs_at_shifted_ips()
{
  shapeDerivsShift_.resize(numIntPoints_*nodesPerElement_*nDim_);
  hex27_shape_deriv(numIntPoints_, intgLocShift_.data(), shapeDerivsShift_.data());
}

//--------------------------------------------------------------------------
//-------- eval_shape_derivs_at_face_ips -----------------------------------
//--------------------------------------------------------------------------
void
HexahedralP2Element::eval_shape_derivs_at_face_ips()
{
  expFaceShapeDerivs_.resize(numIntPoints_*nodesPerElement_*nDim_);
  hex27_shape_deriv(
    numIntPoints_,
    intgExpFace_.data(),
    expFaceShapeDerivs_.data()
  );
}

//--------------------------------------------------------------------------
//-------- hex27_shape_fcn -------------------------------------------------
//--------------------------------------------------------------------------
void
HexahedralP2Element::hex27_shape_fcn(
  int numIntPoints,
  const double *intgLoc,
  double *shpfc) const
{
  const double one = 1.0;
  const double half = 1.0/2.0;
  const double one4th = 1.0/4.0;
  const double one8th = 1.0/8.0;

  for ( int ip = 0; ip < numIntPoints; ++ip ) {
    int ip_offset = nodesPerElement_*ip; // nodes per element is always 27
    int vector_offset = nDim_*ip;

    const double s = intgLoc[vector_offset+0];
    const double t = intgLoc[vector_offset+1];
    const double u = intgLoc[vector_offset+2];

    const double stu = s * t * u;
    const double  st  = s * t;
    const double  su  = s * u;
    const double  tu  = t * u;

    const double one_m_s = one - s;
    const double one_p_s = one + s;
    const double one_m_t = one - t;
    const double  one_p_t = one + t;
    const double  one_m_u = one - u;
    const double  one_p_u = one + u;

    const double  one_m_ss = one - s * s;
    const double one_m_tt = one - t * t;
    const double one_m_uu = one - u * u;

    shpfc[ip_offset+0]  = -one8th * stu * one_m_s  * one_m_t  * one_m_u;
    shpfc[ip_offset+1]  =  one8th * stu * one_p_s  * one_m_t  * one_m_u;
    shpfc[ip_offset+2]  = -one8th * stu * one_p_s  * one_p_t  * one_m_u;
    shpfc[ip_offset+3]  =  one8th * stu * one_m_s  * one_p_t  * one_m_u;
    shpfc[ip_offset+4]  =  one8th * stu * one_m_s  * one_m_t  * one_p_u;
    shpfc[ip_offset+5]  = -one8th * stu * one_p_s  * one_m_t  * one_p_u;
    shpfc[ip_offset+6]  =  one8th * stu * one_p_s  * one_p_t  * one_p_u;
    shpfc[ip_offset+7]  = -one8th * stu * one_m_s  * one_p_t  * one_p_u;
    shpfc[ip_offset+8]  =  one4th * tu  * one_m_ss * one_m_t  * one_m_u;
    shpfc[ip_offset+9]  = -one4th * su  * one_p_s  * one_m_tt * one_m_u;
    shpfc[ip_offset+10] = -one4th * tu  * one_m_ss * one_p_t  * one_m_u;
    shpfc[ip_offset+11] =  one4th * su  * one_m_s  * one_m_tt * one_m_u;
    shpfc[ip_offset+12] =  one4th * st  * one_m_s  * one_m_t  * one_m_uu;
    shpfc[ip_offset+13] = -one4th * st  * one_p_s  * one_m_t  * one_m_uu;
    shpfc[ip_offset+14] =  one4th * st  * one_p_s  * one_p_t  * one_m_uu;
    shpfc[ip_offset+15] = -one4th * st  * one_m_s  * one_p_t  * one_m_uu;
    shpfc[ip_offset+16] = -one4th * tu  * one_m_ss * one_m_t  * one_p_u;
    shpfc[ip_offset+17] =  one4th * su  * one_p_s  * one_m_tt * one_p_u;
    shpfc[ip_offset+18] =  one4th * tu  * one_m_ss * one_p_t  * one_p_u;
    shpfc[ip_offset+19] = -one4th * su  * one_m_s  * one_m_tt * one_p_u;
    shpfc[ip_offset+20] =                 one_m_ss * one_m_tt * one_m_uu;
    shpfc[ip_offset+21] =   -half * u   * one_m_ss * one_m_tt * one_m_u;
    shpfc[ip_offset+22] =    half * u   * one_m_ss * one_m_tt * one_p_u;
    shpfc[ip_offset+23] =   -half * s   * one_m_s  * one_m_tt * one_m_uu;
    shpfc[ip_offset+24] =    half * s   * one_p_s  * one_m_tt * one_m_uu;
    shpfc[ip_offset+25] =   -half * t   * one_m_ss * one_m_t  * one_m_uu;
    shpfc[ip_offset+26] =    half * t   * one_m_ss * one_p_t  * one_m_uu;
  }
}

//--------------------------------------------------------------------------
//-------- hex27_shape_deriv -----------------------------------------------
//--------------------------------------------------------------------------
void
HexahedralP2Element::hex27_shape_deriv(
  int numIntPoints,
  const double *intgLoc,
  double *shapeDerivs) const
{
  const double half = 1.0/2.0;
  const double one4th = 1.0/4.0;
  const double one8th = 1.0/8.0;
  const double two = 2.0;

  for ( int ip = 0; ip < numIntPoints; ++ip ) {
    const int vector_offset = nDim_ * ip;
    const int ip_offset  = nDim_ * nodesPerElement_ * ip;
    int node; int offset;

    const double s = intgLoc[vector_offset+0];
    const double t = intgLoc[vector_offset+1];
    const double u = intgLoc[vector_offset+2];

    const double stu = s * t * u;
    const double st  = s * t;
    const double su  = s * u;
    const double tu  = t * u;

    const double one_m_s = 1.0 - s;
    const double one_p_s = 1.0 + s;
    const double one_m_t = 1.0 - t;
    const double one_p_t = 1.0 + t;
    const double one_m_u = 1.0 - u;
    const double one_p_u = 1.0 + u;

    const double one_m_ss = 1.0 - s * s;
    const double one_m_tt = 1.0 - t * t;
    const double one_m_uu = 1.0 - u * u;

    const double one_m_2s = 1.0 - 2.0 * s;
    const double one_m_2t = 1.0 - 2.0 * t;
    const double one_m_2u = 1.0 - 2.0 * u;

    const double one_p_2s = 1.0 + 2.0 * s;
    const double one_p_2t = 1.0 + 2.0 * t;
    const double one_p_2u = 1.0 + 2.0 * u;

    node = 0;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = -one8th * tu * one_m_2s * one_m_t * one_m_u;
    shapeDerivs[offset + 1] = -one8th * su * one_m_s * one_m_2t * one_m_u;
    shapeDerivs[offset + 2] = -one8th * st * one_m_s * one_m_t * one_m_2u;

    node = 1;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = one8th * tu * one_p_2s * one_m_t * one_m_u;
    shapeDerivs[offset + 1] = one8th * su * one_p_s * one_m_2t * one_m_u;
    shapeDerivs[offset + 2] = one8th * st * one_p_s * one_m_t * one_m_2u;

    node = 2;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = -one8th * tu * one_p_2s * one_p_t * one_m_u;
    shapeDerivs[offset + 1] = -one8th * su * one_p_s * one_p_2t * one_m_u;
    shapeDerivs[offset + 2] = -one8th * st * one_p_s * one_p_t * one_m_2u;

    node = 3;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = one8th * tu * one_m_2s * one_p_t * one_m_u;
    shapeDerivs[offset + 1] = one8th * su * one_m_s * one_p_2t * one_m_u;
    shapeDerivs[offset + 2] = one8th * st * one_m_s * one_p_t * one_m_2u;

    node = 4;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = one8th * tu * one_m_2s * one_m_t * one_p_u;
    shapeDerivs[offset + 1] = one8th * su * one_m_s * one_m_2t * one_p_u;
    shapeDerivs[offset + 2] = one8th * st * one_m_s * one_m_t * one_p_2u;

    node = 5;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = -one8th * tu * one_p_2s * one_m_t * one_p_u;
    shapeDerivs[offset + 1] = -one8th * su * one_p_s * one_m_2t * one_p_u;
    shapeDerivs[offset + 2] = -one8th * st * one_p_s * one_m_t * one_p_2u;

    node = 6;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = one8th * tu * one_p_2s * one_p_t * one_p_u;
    shapeDerivs[offset + 1] = one8th * su * one_p_s * one_p_2t * one_p_u;
    shapeDerivs[offset + 2] = one8th * st * one_p_s * one_p_t * one_p_2u;

    node = 7;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = -one8th * tu * one_m_2s * one_p_t * one_p_u;
    shapeDerivs[offset + 1] = -one8th * su * one_m_s * one_p_2t * one_p_u;
    shapeDerivs[offset + 2] = -one8th * st * one_m_s * one_p_t * one_p_2u;

    node = 8;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = -half * stu * one_m_t * one_m_u;
    shapeDerivs[offset + 1] = one4th * u * one_m_ss * one_m_2t * one_m_u;
    shapeDerivs[offset + 2] = one4th * t * one_m_ss * one_m_t * one_m_2u;

    node = 9;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = -one4th * u * one_p_2s * one_m_tt * one_m_u;
    shapeDerivs[offset + 1] = half * stu * one_p_s * one_m_u;
    shapeDerivs[offset + 2] = -one4th * s * one_p_s * one_m_tt * one_m_2u;

    node = 10;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = half * stu * one_p_t * one_m_u;
    shapeDerivs[offset + 1] = -one4th * u * one_m_ss * one_p_2t * one_m_u;
    shapeDerivs[offset + 2] = -one4th * t * one_m_ss * one_p_t * one_m_2u;

    node = 11;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = one4th * u * one_m_2s * one_m_tt * one_m_u;
    shapeDerivs[offset + 1] = -half * stu * one_m_s * one_m_u;
    shapeDerivs[offset + 2] = one4th * s * one_m_s * one_m_tt * one_m_2u;

    node = 12;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = one4th * t * one_m_2s * one_m_t * one_m_uu;
    shapeDerivs[offset + 1] = one4th * s * one_m_s * one_m_2t * one_m_uu;
    shapeDerivs[offset + 2] = -half * stu * one_m_s * one_m_t;

    node = 13;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = -one4th * t * one_p_2s * one_m_t * one_m_uu;
    shapeDerivs[offset + 1] = -one4th * s * one_p_s * one_m_2t * one_m_uu;
    shapeDerivs[offset + 2] = half * stu * one_p_s * one_m_t;

    node = 14;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = one4th * t * one_p_2s * one_p_t * one_m_uu;
    shapeDerivs[offset + 1] = one4th * s * one_p_s * one_p_2t * one_m_uu;
    shapeDerivs[offset + 2] = -half * stu * one_p_s * one_p_t;

    node = 15;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = -one4th * t * one_m_2s * one_p_t * one_m_uu;
    shapeDerivs[offset + 1] = -one4th * s * one_m_s * one_p_2t * one_m_uu;
    shapeDerivs[offset + 2] = half * stu * one_m_s * one_p_t;

    node = 16;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = half * stu * one_m_t * one_p_u;
    shapeDerivs[offset + 1] = -one4th * u * one_m_ss * one_m_2t * one_p_u;
    shapeDerivs[offset + 2] = -one4th * t * one_m_ss * one_m_t * one_p_2u;

    node = 17;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = one4th * u * one_p_2s * one_m_tt * one_p_u;
    shapeDerivs[offset + 1] = -half * stu * one_p_s * one_p_u;
    shapeDerivs[offset + 2] = one4th * s * one_p_s * one_m_tt * one_p_2u;

    node = 18;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = -half * stu * one_p_t * one_p_u;
    shapeDerivs[offset + 1] = one4th * u * one_m_ss * one_p_2t * one_p_u;
    shapeDerivs[offset + 2] = one4th * t * one_m_ss * one_p_t * one_p_2u;

    node = 19;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = -one4th * u * one_m_2s * one_m_tt * one_p_u;
    shapeDerivs[offset + 1] = half * stu * one_m_s * one_p_u;
    shapeDerivs[offset + 2] = -one4th * s * one_m_s * one_m_tt * one_p_2u;

    node = 20;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = -two * s * one_m_tt * one_m_uu;
    shapeDerivs[offset + 1] = -two * t * one_m_ss * one_m_uu;
    shapeDerivs[offset + 2] = -two * u * one_m_ss * one_m_tt;

    node = 21;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = su * one_m_tt * one_m_u;
    shapeDerivs[offset + 1] = tu * one_m_ss * one_m_u;
    shapeDerivs[offset + 2] = -half * one_m_ss * one_m_tt * one_m_2u;

    node = 22;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = -su * one_m_tt * one_p_u;
    shapeDerivs[offset + 1] = -tu * one_m_ss * one_p_u;
    shapeDerivs[offset + 2] = half * one_m_ss * one_m_tt * one_p_2u;

    node = 23;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = -half * one_m_2s * one_m_tt * one_m_uu;
    shapeDerivs[offset + 1] = st * one_m_s * one_m_uu;
    shapeDerivs[offset + 2] = su * one_m_s * one_m_tt;

    node = 24;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = half * one_p_2s * one_m_tt * one_m_uu;
    shapeDerivs[offset + 1] = -st * one_p_s * one_m_uu;
    shapeDerivs[offset + 2] = -su * one_p_s * one_m_tt;

    node = 25;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = st * one_m_t * one_m_uu;
    shapeDerivs[offset + 1] = -half * one_m_ss * one_m_2t * one_m_uu;
    shapeDerivs[offset + 2] = tu * one_m_ss * one_m_t;

    node = 26;
    offset = ip_offset + nDim_ * node;
    shapeDerivs[offset + 0] = -st * one_p_t * one_m_uu;
    shapeDerivs[offset + 1] = half * one_m_ss * one_p_2t * one_m_uu;
    shapeDerivs[offset + 2] = -tu * one_m_ss * one_p_t;
  }
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Hex27SCV::Hex27SCV()
  : HexahedralP2Element()
{
  // set up the one-dimensional quadrature rule
  set_quadrature_rule();

  // set up integration rule and relevant maps for scvs
  set_interior_info();

  // compute and save shape functions and derivatives at ips
  eval_shape_functions_at_ips();
  eval_shape_derivs_at_ips();
  eval_shape_functions_at_shifted_ips();
  eval_shape_derivs_at_shifted_ips();
}

//--------------------------------------------------------------------------
//-------- set_interior_info -----------------------------------------------
//--------------------------------------------------------------------------
void
Hex27SCV::set_interior_info()
{
  //1D integration rule per sub-control volume
  numIntPoints_ = (nodes1D_ * nodes1D_  * nodes1D_) * ( numQuad_ * numQuad_ * numQuad_); // 216

  // define ip node mappings
  ipNodeMap_.resize(numIntPoints_);
  intgLoc_.resize(numIntPoints_*nDim_);
  intgLocShift_.resize(numIntPoints_*nDim_);
  ipWeight_.resize(numIntPoints_);

  // tensor product nodes (3x3x3) x tensor product quadrature (2x2x2)
  int vector_index = 0; int scalar_index = 0;
  for (int n = 0; n < nodes1D_; ++n) {
    for (int m = 0; m < nodes1D_; ++m) {
      for (int l = 0; l < nodes1D_; ++l) {

        // current node number
        const int nodeNumber = tensor_product_node_map(l,m,n);

        //tensor-product quadrature for a particular sub-cv
        for (int k = 0; k < numQuad_; ++k) {
          for (int j = 0; j < numQuad_; ++j) {
            for (int i = 0; i < numQuad_; ++i) {
              //integration point location
              intgLoc_[vector_index]     = gauss_point_location(l,i);
              intgLoc_[vector_index + 1] = gauss_point_location(m,j);
              intgLoc_[vector_index + 2] = gauss_point_location(n,k);

              intgLocShift_[vector_index]     = shifted_gauss_point_location(l,i);
              intgLocShift_[vector_index + 1] = shifted_gauss_point_location(m,j);
              intgLocShift_[vector_index + 2] = shifted_gauss_point_location(n,k);

              //weight
              ipWeight_[scalar_index] = tensor_product_weight(l,m,n,i,j,k);

              //sub-control volume association
              ipNodeMap_[scalar_index] = nodeNumber;

              // increment indices
              ++scalar_index;
              vector_index += nDim_;
            }
          }
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Hex27SCV::ipNodeMap(
  int /*ordinal*/)
{
  // define scv->node mappings
  return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void Hex27SCV::determinant(
  const int nelem,
  const double *coords,
  double *volume,
  double *error)
{
//  std::vector<double> jacobian(nDim_*nDim_);
  for (int k = 0; k < nelem; ++k) {
    const int scalar_elem_offset = numIntPoints_ * k;
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;

      //weighted jacobian determinant
      const double det_j = jacobian_determinant(&coords[coord_elem_offset],&shapeDerivs_[grad_offset]);

      //apply weight and store to volume
      volume[scalar_elem_offset + ip] = ipWeight_[ip] * det_j;

      //flag error
      if (det_j <= 0.0) {
        *error = 1.0;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- jacobian_determinant---------------------------------------------
//--------------------------------------------------------------------------
double
Hex27SCV::jacobian_determinant(
  const double *elemNodalCoords,
  const double *shapeDerivs) const
{
  double dx_ds1 = 0.0;  double dx_ds2 = 0.0; double dx_ds3 = 0.0;
  double dy_ds1 = 0.0;  double dy_ds2 = 0.0; double dy_ds3 = 0.0;
  double dz_ds1 = 0.0;  double dz_ds2 = 0.0; double dz_ds3 = 0.0;
  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = nDim_ * node;

    const double xCoord = elemNodalCoords[vector_offset+0];
    const double yCoord = elemNodalCoords[vector_offset+1];
    const double zCoord = elemNodalCoords[vector_offset+2];

    const double dn_ds1 = shapeDerivs[vector_offset+0];
    const double dn_ds2 = shapeDerivs[vector_offset+1];
    const double dn_ds3 = shapeDerivs[vector_offset+2];

    dx_ds1 += dn_ds1 * xCoord;
    dx_ds2 += dn_ds2 * xCoord;
    dx_ds3 += dn_ds3 * xCoord;

    dy_ds1 += dn_ds1 * yCoord;
    dy_ds2 += dn_ds2 * yCoord;
    dy_ds3 += dn_ds3 * yCoord;

    dz_ds1 += dn_ds1 * zCoord;
    dz_ds2 += dn_ds2 * zCoord;
    dz_ds3 += dn_ds3 * zCoord;
  }

  const double det_j = dx_ds1 * ( dy_ds2 * dz_ds3 - dz_ds2 * dy_ds3 )
                     + dy_ds1 * ( dz_ds2 * dx_ds3 - dx_ds2 * dz_ds3 )
                     + dz_ds1 * ( dx_ds2 * dy_ds3 - dy_ds2 * dx_ds3 );

  return det_j;
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Hex27SCS::Hex27SCS()
  : HexahedralP2Element()
{
  // set up the one-dimensional quadrature rule
  set_quadrature_rule();

  // set up integration rule and relevant maps on scs
  set_interior_info();

  // set up integration rule and relevant maps on faces
  set_boundary_info();

  // compute and save shape functions and derivatives at ips
  eval_shape_functions_at_ips();
  eval_shape_derivs_at_ips();
  eval_shape_derivs_at_face_ips();
  eval_shape_functions_at_shifted_ips();
  eval_shape_derivs_at_shifted_ips();
}

//--------------------------------------------------------------------------
//-------- set_interior_info -----------------------------------------------
//--------------------------------------------------------------------------
void
Hex27SCS::set_interior_info()
{
  const int surfacesPerDirection = nodes1D_ - 1; // 2
  const int ipsPerSurface = (numQuad_*numQuad_)*(nodes1D_*nodes1D_); // 36
  const int numSurfaces = surfacesPerDirection * nDim_; // 6

  numIntPoints_ = numSurfaces*ipsPerSurface; // 216
  const int numVectorPoints = numIntPoints_*nDim_; // 648

  // define L/R mappings
  lrscv_.resize(2*numIntPoints_); // size = 432

  // standard integration location
  intgLoc_.resize(numVectorPoints);

  // shifted
  intgLocShift_.resize(numVectorPoints);

  // Save quadrature weight and directionality information
  ipInfo_.resize(numIntPoints_);

  // a list of the scs locations in 1D
  const std::vector<double> scsLoc = { -scsDist_, scsDist_ };

  // correct orientation of area vector
  const std::vector<double> orientation = {-1.0, +1.0};

  // specify integration point locations in a dimension-by-dimension manner
  //u direction: bottom-top (0-1)
  int vector_index = 0; int lrscv_index = 0; int scalar_index = 0;
  for (int m = 0; m < surfacesPerDirection; ++m) {
    for (int l = 0; l < nodes1D_; ++l) {
      for (int k = 0; k < nodes1D_; ++k) {

        int leftNode; int rightNode;
        if (m == 0) {
          leftNode = tensor_product_node_map(k,l,m);
          rightNode = tensor_product_node_map(k,l,m+1);
        }
        else {
          leftNode = tensor_product_node_map(k,l,m+1);
          rightNode = tensor_product_node_map(k,l,m);
        }

        for (int j = 0; j < numQuad_; ++j) {
          for (int i = 0; i < numQuad_; ++i) {
            lrscv_[lrscv_index]     = leftNode;
            lrscv_[lrscv_index + 1] = rightNode;

            intgLoc_[vector_index]     = gauss_point_location(k,i);
            intgLoc_[vector_index + 1] = gauss_point_location(l,j);
            intgLoc_[vector_index + 2] = scsLoc[m];

            intgLocShift_[vector_index]     = shifted_gauss_point_location(k,i);
            intgLocShift_[vector_index + 1] = shifted_gauss_point_location(l,j);
            intgLocShift_[vector_index + 2] = scsLoc[m];

            //compute the quadrature weight
            ipInfo_[scalar_index].weight = orientation[m] * tensor_product_weight(k,l,i,j);

            //direction
            ipInfo_[scalar_index].direction = Jacobian::U_DIRECTION;

            ++scalar_index;
            lrscv_index += 2;
            vector_index += nDim_;
          }
        }
      }
    }
  }

  //t direction: front-back (2-3)
  for (int m = 0; m < surfacesPerDirection; ++m) {
    for (int l = 0; l < nodes1D_; ++l) {
      for (int k = 0; k < nodes1D_; ++k) {

        int leftNode; int rightNode;
        if (m == 0) {
          leftNode = tensor_product_node_map(k,m,l);
          rightNode = tensor_product_node_map(k,m+1,l);
        }
        else {
          leftNode = tensor_product_node_map(k,m+1,l);
          rightNode = tensor_product_node_map(k,m,l);
        }

        for (int j = 0; j < numQuad_; ++j) {
          for (int i = 0; i < numQuad_; ++i) {
            lrscv_[lrscv_index]     = leftNode;
            lrscv_[lrscv_index + 1] = rightNode;

            intgLoc_[vector_index]     = gauss_point_location(k,i);
            intgLoc_[vector_index + 1] = scsLoc[m];
            intgLoc_[vector_index + 2] = gauss_point_location(l,j);

            intgLocShift_[vector_index]     = shifted_gauss_point_location(k,i);
            intgLocShift_[vector_index + 1] = scsLoc[m];
            intgLocShift_[vector_index + 2] = shifted_gauss_point_location(l,j);

            //compute the quadrature weight
            ipInfo_[scalar_index].weight = orientation[m] * tensor_product_weight(k,l,i,j);

            //direction
            ipInfo_[scalar_index].direction = Jacobian::T_DIRECTION;

            ++scalar_index;
            lrscv_index += 2;
            vector_index += nDim_;
          }
        }
      }
    }
  }

  //s direction: left-right (4-5)
  for (int m = 0; m < surfacesPerDirection; ++m) {
    for (int l = 0; l < nodes1D_; ++l) {
      for (int k = 0; k < nodes1D_; ++k) {

        int leftNode; int rightNode;
        if (m == 0) {
          leftNode = tensor_product_node_map(m,k,l);
          rightNode = tensor_product_node_map(m+1,k,l);
        }
        else {
          leftNode = tensor_product_node_map(m+1,k,l);
          rightNode = tensor_product_node_map(m,k,l);
        }

        for (int j = 0; j < numQuad_; ++j) {
          for (int i = 0; i < numQuad_; ++i) {
            lrscv_[lrscv_index]     = leftNode;
            lrscv_[lrscv_index + 1] = rightNode;

            intgLoc_[vector_index]     = scsLoc[m];
            intgLoc_[vector_index + 1] = gauss_point_location(k,i);
            intgLoc_[vector_index + 2] = gauss_point_location(l,j);

            intgLocShift_[vector_index]     = scsLoc[m];
            intgLocShift_[vector_index + 1] = shifted_gauss_point_location(k,i);
            intgLocShift_[vector_index + 2] = shifted_gauss_point_location(l,j);

            //compute the quadrature weight
            ipInfo_[scalar_index].weight = -orientation[m] * tensor_product_weight(k,l,i,j);

            //direction
            ipInfo_[scalar_index].direction = Jacobian::S_DIRECTION;

            ++scalar_index;
            lrscv_index += 2;
            vector_index += nDim_;
          }
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- set_boundary_info -----------------------------------------------
//--------------------------------------------------------------------------
void
Hex27SCS::set_boundary_info()
{
  const int numFaces = 2 * nDim_; // 6
  const int nodesPerFace = nodes1D_ * nodes1D_; // 9
  ipsPerFace_ = nodesPerFace * (numQuad_ * numQuad_); // 36
  const int numFaceIps = numFaces * ipsPerFace_; // 216 = numIntPoints_ for this element

  oppFace_.resize(numFaceIps);
  ipNodeMap_.resize(numFaceIps);
  oppNode_.resize(numFaceIps);
  intgExpFace_.resize(numFaceIps*nDim_); // size = 648

  // face ordinal to tensor-product style node ordering
  const std::vector<int> stkFaceNodeMap = {
                                            0,  8,  1, 12, 25, 13,  4, 16,  5, // face 0(2): front face (cclockwise)
                                            1,  9,  2, 13, 24, 14,  5, 17,  6, // face 1(5): right face (cclockwise)
                                            3, 10,  2, 15, 26, 14,  7, 18,  6, // face 2(3): back face  (clockwise)
                                            0, 11,  3, 12, 23, 15,  4, 19,  7, // face 3(4): left face  (clockwise)
                                            0,  8,  1, 11, 21, 9,   3, 10,  2, // face 4(0): bottom face (clockwise)
                                            4, 16,  5, 19, 22,  17, 7, 18,  6  // face 5(1): top face (cclockwise)
                                          };


  // tensor-product style access to the map
  auto face_node_number = [=] (int i, int j, int faceOrdinal)
  {
    return stkFaceNodeMap[i + nodes1D_ * j + nodesPerFace * faceOrdinal];
  };

  // map face ip ordinal to nearest sub-control surface ip ordinal
  // sub-control surface renumbering
  const std::vector<int> faceToSurface = { 2, 5, 3, 4, 0, 1 };
  auto opp_face_map = [=] ( int k, int l, int i, int j, int face_index)
  {
    int face_offset = faceToSurface[face_index] * ipsPerFace_;

    int node_index = k + nodes1D_ * l;
    int node_offset = node_index * (numQuad_ * numQuad_);

    int ip_index = face_offset+node_offset+i+numQuad_*j;

    return ip_index;
  };

  // location of the faces in the correct order
  const std::vector<double> faceLoc = {-1.0, +1.0, +1.0, -1.0, -1.0, +1.0};

  // Set points face-by-face
  int vector_index = 0; int scalar_index = 0; int faceOrdinal = 0;

  // front face: t = -1.0: counter-clockwise
  faceOrdinal = 0;
  for (int l = 0; l < nodes1D_; ++l) {
    for (int k = 0; k < nodes1D_; ++k) {
      const int nearNode = face_node_number(k,l,faceOrdinal);
      int oppNode = tensor_product_node_map(k,1,l);

      //tensor-product quadrature for a particular sub-cv
      for (int j = 0; j < numQuad_; ++j) {
        for (int i = 0; i < numQuad_; ++i) {
          // set maps
          ipNodeMap_[scalar_index] = nearNode;
          oppNode_[scalar_index] = oppNode;
          oppFace_[scalar_index] = opp_face_map(k,l,i,j,faceOrdinal);

          //integration point location
          intgExpFace_[vector_index]     = intgLoc_[oppFace_[scalar_index]*nDim_+0];
          intgExpFace_[vector_index + 1] = faceLoc[faceOrdinal];
          intgExpFace_[vector_index + 2] = intgLoc_[oppFace_[scalar_index]*nDim_+2];

          // increment indices
          ++scalar_index;
          vector_index += nDim_;
        }
      }
    }
  }

  // right face: s = +1.0: counter-clockwise
  faceOrdinal = 1;
  for (int l = 0; l < nodes1D_; ++l) {
    for (int k = 0; k < nodes1D_; ++k) {
      const int nearNode = face_node_number(k,l,faceOrdinal);
      int oppNode = tensor_product_node_map(1,k,l);

      //tensor-product quadrature for a particular sub-cv
      for (int j = 0; j < numQuad_; ++j) {
        for (int i = 0; i < numQuad_; ++i) {
          // set maps
          ipNodeMap_[scalar_index] = nearNode;
          oppNode_[scalar_index] = oppNode;
          oppFace_[scalar_index] = opp_face_map(k,l,i,j,faceOrdinal);

          //integration point location
          intgExpFace_[vector_index]     = faceLoc[faceOrdinal];
          intgExpFace_[vector_index + 1] = intgLoc_[oppFace_[scalar_index]*nDim_+1];
          intgExpFace_[vector_index + 2] = intgLoc_[oppFace_[scalar_index]*nDim_+2];

          // increment indices
          ++scalar_index;
          vector_index += nDim_;
        }
      }
    }
  }

  // back face: s = +1.0: s-direction reversed
  faceOrdinal = 2;
  for (int l = 0; l < nodes1D_; ++l) {
    for (int k = nodes1D_-1; k >= 0; --k) {
      const int nearNode = face_node_number(k,l,faceOrdinal);
      int oppNode = tensor_product_node_map(k,1,l);

      //tensor-product quadrature for a particular sub-cv
      for (int j = 0; j < numQuad_; ++j) {
        for (int i = numQuad_-1; i >= 0; --i) {
          // set maps
          ipNodeMap_[scalar_index] = nearNode;
          oppNode_[scalar_index] = oppNode;
          oppFace_[scalar_index] = opp_face_map(k,l,i,j,faceOrdinal);

          //integration point location
          intgExpFace_[vector_index]     = intgLoc_[oppFace_[scalar_index]*nDim_+0];
          intgExpFace_[vector_index + 1] = faceLoc[faceOrdinal];
          intgExpFace_[vector_index + 2] = intgLoc_[oppFace_[scalar_index]*nDim_+2];

          // increment indices
          ++scalar_index;
          vector_index += nDim_;
        }
      }
    }
  }

  //left face: x = -1.0 swapped t and u
  faceOrdinal = 3;
  for (int l = 0; l < nodes1D_; ++l) {
    for (int k = 0; k < nodes1D_; ++k) {
      const int nearNode = face_node_number(l,k,faceOrdinal);
      int oppNode = tensor_product_node_map(1,l,k);

      //tensor-product quadrature for a particular sub-cv
      for (int j = 0; j < numQuad_; ++j) {
        for (int i = 0; i < numQuad_; ++i) {
          // set maps
          ipNodeMap_[scalar_index] = nearNode;
          oppNode_[scalar_index]   = oppNode;
          oppFace_[scalar_index]   = opp_face_map(l,k,j,i,faceOrdinal);

          //integration point location
          intgExpFace_[vector_index]     = faceLoc[faceOrdinal];
          intgExpFace_[vector_index + 1] = intgLoc_[oppFace_[scalar_index]*nDim_+1];
          intgExpFace_[vector_index + 2] = intgLoc_[oppFace_[scalar_index]*nDim_+2];

          // increment indices
          ++scalar_index;
          vector_index += nDim_;
        }
      }
    }
  }

  //bottom face: u = -1.0: swapped s and t
  faceOrdinal = 4;
  for (int l = 0; l < nodes1D_; ++l) {
    for (int k = 0; k < nodes1D_; ++k) {
      const int nearNode = face_node_number(l,k,faceOrdinal);
      int oppNode = tensor_product_node_map(l,k,1);

      //tensor-product quadrature for a particular sub-cv
      for (int j = 0; j < numQuad_; ++j) {
        for (int i = 0; i < numQuad_; ++i) {
          // set maps
          ipNodeMap_[scalar_index] = nearNode;
          oppNode_[scalar_index] = oppNode;
          oppFace_[scalar_index] = opp_face_map(l,k,j,i,faceOrdinal);

          //integration point location
          intgExpFace_[vector_index]     = intgLoc_[oppFace_[scalar_index]*nDim_+0];
          intgExpFace_[vector_index + 1] = intgLoc_[oppFace_[scalar_index]*nDim_+1];
          intgExpFace_[vector_index + 2] = faceLoc[faceOrdinal];

          // increment indices
          ++scalar_index;
          vector_index += nDim_;
        }
      }
    }
  }

  //top face: u = +1.0: counter-clockwise
  faceOrdinal = 5;
  for (int l = 0; l < nodes1D_; ++l) {
    for (int k = 0; k < nodes1D_; ++k) {
      const int nearNode = face_node_number(k,l,faceOrdinal);
      int oppNode = tensor_product_node_map(k,l,1);

      //tensor-product quadrature for a particular sub-cv
      for (int j = 0; j < numQuad_; ++j) {
        for (int i = 0; i < numQuad_; ++i) {
          // set maps
          ipNodeMap_[scalar_index] = nearNode;
          oppNode_[scalar_index] = oppNode;
          oppFace_[scalar_index] = opp_face_map(k,l,i,j,faceOrdinal);

          //integration point location
          intgExpFace_[vector_index]     = intgLoc_[oppFace_[scalar_index]*nDim_+0];
          intgExpFace_[vector_index + 1] = intgLoc_[oppFace_[scalar_index]*nDim_+1];
          intgExpFace_[vector_index + 2] = faceLoc[faceOrdinal];

          // increment indices
          ++scalar_index;
          vector_index += nDim_;
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- adjacentNodes ---------------------------------------------------
//--------------------------------------------------------------------------
const int *
Hex27SCS::adjacentNodes()
{
  // define L/R mappings
  return &lrscv_[0];
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Hex27SCS::ipNodeMap(
  int ordinal)
{
  // define ip->node mappings for each face (ordinal);
  return &ipNodeMap_[ordinal*ipsPerFace_];
}

//--------------------------------------------------------------------------
//-------- opposingNodes --------------------------------------------------
//--------------------------------------------------------------------------
int
Hex27SCS::opposingNodes(
  const int ordinal,
  const int node)
{
  return oppNode_[ordinal*ipsPerFace_+node];
}

//--------------------------------------------------------------------------
//-------- opposingFace --------------------------------------------------
//--------------------------------------------------------------------------
int
Hex27SCS::opposingFace(
  const int ordinal,
  const int node)
{
  return oppFace_[ordinal*ipsPerFace_+node];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void
Hex27SCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  //returns the normal vector x_t x x_u for constant s curves
  //returns the normal vector x_u x x_s for constant t curves
  //returns the normal vector x_s x x_t for constant u curves
  std::array<double,3> areaVector;

  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
        const int vector_elem_offset = nDim_ * numIntPoints_ * k;

    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;
      const int offset = nDim_ * ip + vector_elem_offset;

      //compute area vector for this ip
      area_vector( ipInfo_[ip].direction,
                   &coords[coord_elem_offset],
                   &shapeDerivs_[grad_offset],
                   areaVector.data() );

      // apply quadrature weight and orientation (combined as weight)
      for (int j = 0; j < nDim_; ++j) {
        areav[offset+j]  = ipInfo_[ip].weight * areaVector[j];
      }
    }
  }

  *error = 0;
}

//--------------------------------------------------------------------------
//-------- area_vector -----------------------------------------------------
//--------------------------------------------------------------------------
void
Hex27SCS::area_vector(
  const Jacobian::Direction direction,
  const double *elemNodalCoords,
  double *shapeDeriv,
  double *areaVector) const
{

  int s1Component; int s2Component;
  switch (direction) {
    case Jacobian::S_DIRECTION:
      s1Component = static_cast<int>(Jacobian::T_DIRECTION);
      s2Component = static_cast<int>(Jacobian::U_DIRECTION);
      break;
    case Jacobian::T_DIRECTION:
      s1Component = static_cast<int>(Jacobian::S_DIRECTION);
      s2Component = static_cast<int>(Jacobian::U_DIRECTION);
      break;
    case Jacobian::U_DIRECTION:
      s1Component = static_cast<int>(Jacobian::T_DIRECTION);
      s2Component = static_cast<int>(Jacobian::S_DIRECTION);
      break;
    default:
      throw std::runtime_error("Not a valid direction for this element!");
  }

  // return the normal area vector given shape derivatives dnds OR dndt
  double dx_ds1 = 0.0; double dy_ds1 = 0.0; double dz_ds1 = 0.0;
  double dx_ds2 = 0.0; double dy_ds2 = 0.0; double dz_ds2 = 0.0;

  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = nDim_ * node;
    const double xCoord = elemNodalCoords[vector_offset+0];
    const double yCoord = elemNodalCoords[vector_offset+1];
    const double zCoord = elemNodalCoords[vector_offset+2];

    const double dn_ds1 = shapeDeriv[vector_offset+s1Component];
    const double dn_ds2 = shapeDeriv[vector_offset+s2Component];

    dx_ds1 += dn_ds1 * xCoord;
    dx_ds2 += dn_ds2 * xCoord;

    dy_ds1 += dn_ds1 * yCoord;
    dy_ds2 += dn_ds2 * yCoord;

    dz_ds1 += dn_ds1 * zCoord;
    dz_ds2 += dn_ds2 * zCoord;
  }

  //cross product
  areaVector[0] = dy_ds1*dz_ds2 - dz_ds1*dy_ds2;
  areaVector[1] = dz_ds1*dx_ds2 - dx_ds1*dz_ds2;
  areaVector[2] = dx_ds1*dy_ds2 - dy_ds1*dx_ds2;
}

//--------------------------------------------------------------------------
//-------- grad_op ---------------------------------------------------------
//--------------------------------------------------------------------------
void Hex27SCS::grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  *error = 0.0;

  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    const int scalar_elem_offset = numIntPoints_ * k;
    const int grad_elem_offset = numIntPoints_ * nDim_ * nodesPerElement_ * k;

    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;
      const int offset = grad_offset + grad_elem_offset;

      for (int j = 0; j < nodesPerElement_ * nDim_; ++j) {
        deriv[offset + j] = shapeDerivs_[grad_offset +j];
      }

      gradient( &coords[coord_elem_offset],
                &shapeDerivs_[grad_offset],
                &gradop[offset],
                &det_j[scalar_elem_offset+ip] );

      if (det_j[ip] <= 0.0) {
        *error = 1.0;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- shifted_grad_op -------------------------------------------------
//--------------------------------------------------------------------------
void Hex27SCS::shifted_grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  *error = 0.0;

  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    const int scalar_elem_offset = numIntPoints_ * k;
    const int grad_elem_offset = numIntPoints_ * nDim_ * nodesPerElement_ * k;

    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;
      const int offset = grad_offset + grad_elem_offset;

      for (int j = 0; j < nodesPerElement_ * nDim_; ++j) {
        deriv[offset + j] = shapeDerivsShift_[grad_offset +j];
      }

      gradient( &coords[coord_elem_offset],
                &shapeDerivsShift_[grad_offset],
                &gradop[offset],
                &det_j[scalar_elem_offset+ip] );

      if (det_j[ip] <= 0.0) {
        *error = 1.0;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- face_grad_op ----------------------------------------------------
//--------------------------------------------------------------------------
void Hex27SCS::face_grad_op(
  const int nelem,
  const int face_ordinal,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  *error = 0.0;

  const int face_offset =  nDim_ * ipsPerFace_ * nodesPerElement_ * face_ordinal;
  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    const int scalar_elem_offset = ipsPerFace_ * k;
    const int grad_elem_offset = ipsPerFace_ * nDim_ * nodesPerElement_ * k;

    for (int ip = 0; ip < ipsPerFace_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;
      const int offset = grad_offset + grad_elem_offset;

      gradient( &coords[coord_elem_offset],
                &expFaceShapeDerivs_[face_offset+grad_offset],
                &gradop[offset],
                &det_j[scalar_elem_offset+ip] );

      if (det_j[ip] <= 0.0) {
        *error = 1.0;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- gradient --------------------------------------------------------
//--------------------------------------------------------------------------
void
Hex27SCS::gradient(
  const double* elemNodalCoords,
  const double* shapeDeriv,
  double* grad,
  double* det_j) const
{
  double dx_ds1 = 0.0;  double dx_ds2 = 0.0; double dx_ds3 = 0.0;
  double dy_ds1 = 0.0;  double dy_ds2 = 0.0; double dy_ds3 = 0.0;
  double dz_ds1 = 0.0;  double dz_ds2 = 0.0; double dz_ds3 = 0.0;

  //compute Jacobian
  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = nDim_ * node;

    const double xCoord = elemNodalCoords[vector_offset + 0];
    const double yCoord = elemNodalCoords[vector_offset + 1];
    const double zCoord = elemNodalCoords[vector_offset + 2];

    const double dn_ds1 = shapeDeriv[vector_offset + 0];
    const double dn_ds2 = shapeDeriv[vector_offset + 1];
    const double dn_ds3 = shapeDeriv[vector_offset + 2];

    dx_ds1 += dn_ds1 * xCoord;
    dx_ds2 += dn_ds2 * xCoord;
    dx_ds3 += dn_ds3 * xCoord;

    dy_ds1 += dn_ds1 * yCoord;
    dy_ds2 += dn_ds2 * yCoord;
    dy_ds3 += dn_ds3 * yCoord;

    dz_ds1 += dn_ds1 * zCoord;
    dz_ds2 += dn_ds2 * zCoord;
    dz_ds3 += dn_ds3 * zCoord;
  }

  *det_j = dx_ds1 * ( dy_ds2 * dz_ds3 - dz_ds2 * dy_ds3 )
         + dy_ds1 * ( dz_ds2 * dx_ds3 - dx_ds2 * dz_ds3 )
         + dz_ds1 * ( dx_ds2 * dy_ds3 - dy_ds2 * dx_ds3 );

  const double inv_det_j = (*det_j > 0.0) ? 1.0 / (*det_j) : 0.0;

  const double ds1_dx = inv_det_j*(dy_ds2 * dz_ds3 - dz_ds2 * dy_ds3);
  const double ds2_dx = inv_det_j*(dz_ds1 * dy_ds3 - dy_ds1 * dz_ds3);
  const double ds3_dx = inv_det_j*(dy_ds1 * dz_ds2 - dz_ds1 * dy_ds2);

  const double ds1_dy = inv_det_j*(dz_ds2 * dx_ds3 - dx_ds2 * dz_ds3);
  const double ds2_dy = inv_det_j*(dx_ds1 * dz_ds3 - dz_ds1 * dx_ds3);
  const double ds3_dy = inv_det_j*(dz_ds1 * dx_ds2 - dx_ds1 * dz_ds2);

  const double ds1_dz = inv_det_j*(dx_ds2 * dy_ds3 - dy_ds2 * dx_ds3);
  const double ds2_dz = inv_det_j*(dy_ds1 * dx_ds3 - dx_ds1 * dy_ds3);
  const double ds3_dz = inv_det_j*(dx_ds1 * dy_ds2 - dy_ds1 * dx_ds2);

  // metrics
  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = nDim_ * node;

    const double dn_ds1 = shapeDeriv[vector_offset + 0];
    const double dn_ds2 = shapeDeriv[vector_offset + 1];
    const double dn_ds3 = shapeDeriv[vector_offset + 2];

    grad[vector_offset + 0] = dn_ds1 * ds1_dx + dn_ds2 * ds2_dx + dn_ds3 * ds3_dx;
    grad[vector_offset + 1] = dn_ds1 * ds1_dy + dn_ds2 * ds2_dy + dn_ds3 * ds3_dy;
    grad[vector_offset + 2] = dn_ds1 * ds1_dz + dn_ds2 * ds2_dz + dn_ds3 * ds3_dz;
  }
}

//--------------------------------------------------------------------------
//-------- gij -------------------------------------------------------------
//--------------------------------------------------------------------------
void Hex27SCS::gij(
  const double *coords,
  double *gupperij,
  double *glowerij,
  double *deriv)
{
  SIERRA_FORTRAN(threed_gij)
    ( &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gupperij, glowerij);
}

//--------------------------------------------------------------------------
//-------- general_face_grad_op --------------------------------------------
//--------------------------------------------------------------------------
void 
Hex27SCS::general_face_grad_op(
  const int face_ordinal,
  const double *isoParCoord,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  const int ipsPerFace = 1;
  double faceShapeFuncDerivs[81];

  hex27_shape_deriv(
    ipsPerFace,
    isoParCoord,
    faceShapeFuncDerivs);

  gradient( &coords[0],
            faceShapeFuncDerivs,
            &gradop[0],
            &det_j[0] );
  
  if (det_j[0] <= 0.0) {
    *error = 1.0;
  }
}

//--------------------------------------------------------------------------
//-------- sidePcoords_to_elemPcoords --------------------------------------
//--------------------------------------------------------------------------
void 
Hex27SCS::sidePcoords_to_elemPcoords(
  const int & side_ordinal,
  const int & npoints,
  const double *side_pcoords,
  double *elem_pcoords)
{
  // each ME are -1:1, e.g., hex27:quad93d
  switch (side_ordinal) {
  case 0:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = side_pcoords[2*i+0];
      elem_pcoords[i*3+1] = -1.0;
      elem_pcoords[i*3+2] = side_pcoords[2*i+1];
    }
    break;
  case 1:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = 1.0;
      elem_pcoords[i*3+1] = side_pcoords[2*i+0];
      elem_pcoords[i*3+2] = side_pcoords[2*i+1];
    }
    break;
  case 2:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = -side_pcoords[2*i+0];
      elem_pcoords[i*3+1] = 1.0;
      elem_pcoords[i*3+2] = side_pcoords[2*i+1];
    }
    break;
  case 3:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = -1.0;
      elem_pcoords[i*3+1] = side_pcoords[2*i+1];
      elem_pcoords[i*3+2] = side_pcoords[2*i+0];
    }
    break;
  case 4:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = side_pcoords[2*i+1];
      elem_pcoords[i*3+1] = side_pcoords[2*i+0];
      elem_pcoords[i*3+2] = -1.0;
    }
    break;
  case 5:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = side_pcoords[2*i+0];
      elem_pcoords[i*3+1] = side_pcoords[2*i+1];
      elem_pcoords[i*3+2] = 1.0;
    }
    break;
  default:
    throw std::runtime_error("HexSCS::sideMap invalid ordinal");
  }
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TetSCV::TetSCV()
  : MasterElement()
{
  nDim_ = 3;
  nodesPerElement_ = 4;
  numIntPoints_ = 4;

  // define ip node mappings
  ipNodeMap_.resize(4);
  ipNodeMap_[0] = 0; ipNodeMap_[1] = 1; ipNodeMap_[2] = 2; ipNodeMap_[3] = 3;

  // standard integration location
  intgLoc_.resize(12);
  const double seventeen96ths = 17.0/96.0;
  const double fourfive96ths  = 45.0/96.0;
  intgLoc_[0] = seventeen96ths; intgLoc_[1]  = seventeen96ths; intgLoc_[2]  = seventeen96ths; // vol 1
  intgLoc_[3] = fourfive96ths;  intgLoc_[4]  = seventeen96ths; intgLoc_[5]  = seventeen96ths; // vol 2
  intgLoc_[6] = seventeen96ths; intgLoc_[7]  = fourfive96ths;  intgLoc_[8]  = seventeen96ths; // vol 3
  intgLoc_[9] = seventeen96ths; intgLoc_[10] = seventeen96ths; intgLoc_[11] = fourfive96ths;  // vol 4

  // shifted
  intgLocShift_.resize(12);
  intgLocShift_[0] = 0.0; intgLocShift_[1]  = 0.0;  intgLocShift_[2] = 0.0;
  intgLocShift_[3] = 1.0; intgLocShift_[4]  = 0.0;  intgLocShift_[5] = 0.0;
  intgLocShift_[6] = 0.0; intgLocShift_[7]  = 1.0;  intgLocShift_[8] = 0.0;
  intgLocShift_[9] = 0.0; intgLocShift_[10] = 0.0; intgLocShift_[11] = 1.0;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
TetSCV::~TetSCV()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
TetSCV::ipNodeMap(
  int /*ordinal*/)
{
  // define scv->node mappings
  return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void TetSCV::determinant(
  const int nelem,
  const double *coords,
  double *volume,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(tet_scv_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords,
      volume, error, &lerr );
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
TetSCV::shape_fcn(double *shpfc)
{
  tet_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
TetSCV::shifted_shape_fcn(double *shpfc)
{
  tet_shape_fcn(numIntPoints_, &intgLocShift_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- tet_shape_fcn ---------------------------------------------------
//--------------------------------------------------------------------------
void
TetSCV::tet_shape_fcn(
  const int  &npts,
  const double *par_coord, 
  double *shape_fcn)
{
  for (int j = 0; j < npts; ++j ) {
    const int fourj = 4*j;
    const int k = 3*j;
    const double xi = par_coord[k];
    const double eta = par_coord[k+1];
    const double zeta = par_coord[k+2];
    shape_fcn[fourj] = 1.0 - xi - eta - zeta;
    shape_fcn[1 + fourj] = xi;
    shape_fcn[2 + fourj] = eta;
    shape_fcn[3 + fourj] = zeta;
  }
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TetSCS::TetSCS()
  : MasterElement()
{
  nDim_ = 3;
  nodesPerElement_ = 4;
  numIntPoints_ = 6;

  // define L/R mappings
  lrscv_.resize(12);
  lrscv_[0]  = 0; lrscv_[1]  = 1;
  lrscv_[2]  = 1; lrscv_[3]  = 2;
  lrscv_[4]  = 0; lrscv_[5]  = 2;
  lrscv_[6]  = 0; lrscv_[7]  = 3;
  lrscv_[8]  = 1; lrscv_[9]  = 3;
  lrscv_[10] = 2; lrscv_[11] = 3;

  // define opposing node
  oppNode_.resize(12);
  // face 0
  oppNode_[0] = 2; oppNode_[1] = 2;  oppNode_[2] = 2;
  // face 1
  oppNode_[3] = 0; oppNode_[4] = 0;  oppNode_[5] = 0;
  // face 2
  oppNode_[6] = 1; oppNode_[7] = 1;  oppNode_[8] = 1;
  // face 3
  oppNode_[9] = 3; oppNode_[10] = 3; oppNode_[11] = 3;

  // define opposing face
  oppFace_.resize(12);
  // face 0
  oppFace_[0]  = 2; oppFace_[1] = 1;  oppFace_[2] = 5;
  // face 1
  oppFace_[3]  = 0; oppFace_[4] = 2;  oppFace_[5] = 3;
  // face 2
  oppFace_[6]  = 0; oppFace_[7] = 4;  oppFace_[8] = 1;
  // face 3
  oppFace_[9]  = 3; oppFace_[10] = 5; oppFace_[11] = 4;

  // standard integration location
  intgLoc_.resize(18);
  const double seventeen48ths = 17.0/48.0;
  const double seven48ths = 7.0/48.0;
  intgLoc_[0]  =  seventeen48ths; intgLoc_[1]  = seven48ths;     intgLoc_[2]  = seven48ths; // surf 1    1->2
  intgLoc_[3]  =  seventeen48ths; intgLoc_[4]  = seventeen48ths; intgLoc_[5]  = seven48ths; // surf 2    2->3
  intgLoc_[6]  =  seven48ths;     intgLoc_[7]  = seventeen48ths; intgLoc_[8]  = seven48ths; // surf 3    1->3
  intgLoc_[9]  =  seven48ths ;    intgLoc_[10] = seven48ths;     intgLoc_[11] = seventeen48ths; // surf 4    1->4
  intgLoc_[12] =  seventeen48ths; intgLoc_[13] = seven48ths;     intgLoc_[14] = seventeen48ths; // surf 5    2->4
  intgLoc_[15] =  seven48ths;     intgLoc_[16] = seventeen48ths; intgLoc_[17] = seventeen48ths; // surf 6    3->4

  // shifted
  intgLocShift_.resize(36);
  intgLocShift_[0]  =  0.50; intgLocShift_[1]  =  0.00; intgLocShift_[2]  =  0.00; // surf 1    1->2
  intgLocShift_[3]  =  0.50; intgLocShift_[4]  =  0.50; intgLocShift_[5]  =  0.00; // surf 2    2->3
  intgLocShift_[6]  =  0.00; intgLocShift_[7]  =  0.50; intgLocShift_[8]  =  0.00; // surf 3    1->3
  intgLocShift_[9]  =  0.00; intgLocShift_[10] =  0.00; intgLocShift_[11] =  0.50; // surf 4    1->4
  intgLocShift_[12] =  0.50; intgLocShift_[13] =  0.00; intgLocShift_[14] =  0.50; // surf 5    2->4
  intgLocShift_[15] =  0.00; intgLocShift_[16] =  0.50; intgLocShift_[17] =  0.50; // surf 6    3->4

  // exposed face
  intgExpFace_.resize(36);
  const double five24ths = 5.0/24.0;
  const double seven12ths = 7.0/12.0;
  // face 0; nodes 0,1,3: scs 0, 1, 2
  intgExpFace_[0]  = five24ths;  intgExpFace_[1]  =  0.00; intgExpFace_[2]  = five24ths;
  intgExpFace_[3]  = seven12ths; intgExpFace_[4]  =  0.00; intgExpFace_[5]  = five24ths;
  intgExpFace_[6]  = five24ths;  intgExpFace_[7]  =  0.00; intgExpFace_[8]  = seven12ths;
  // face 1; nodes 1,2,3; scs 0, 1, 2
  intgExpFace_[9]  = seven12ths; intgExpFace_[10] = five24ths;  intgExpFace_[11] = five24ths;
  intgExpFace_[12] = five24ths;  intgExpFace_[13] = seven12ths; intgExpFace_[14] = five24ths;
  intgExpFace_[15] = five24ths;  intgExpFace_[16] = five24ths;  intgExpFace_[17] = seven12ths;
  // face 2; nodes 0,3,2; scs 0, 1, 2
  intgExpFace_[18] =  0.00;      intgExpFace_[19] = five24ths;  intgExpFace_[20] = five24ths;
  intgExpFace_[21] =  0.00;      intgExpFace_[22] = five24ths;  intgExpFace_[23] = seven12ths;
  intgExpFace_[24] =  0.00;      intgExpFace_[25] = seven12ths; intgExpFace_[26] = five24ths;
  //face 3; nodes 0, 2, 1; scs 0, 1, 2
  intgExpFace_[27] = five24ths;  intgExpFace_[28] = five24ths;  intgExpFace_[29] =  0.00;
  intgExpFace_[30] = seven12ths; intgExpFace_[31] = five24ths;  intgExpFace_[32] =  0.00;
  intgExpFace_[33] = five24ths;  intgExpFace_[34] = seven12ths; intgExpFace_[35] =  0.00;

  // boundary integration point ip node mapping (ip on an ordinal to local node number)
  ipNodeMap_.resize(12); // 3 ips * 4 faces
  // face 0;
  ipNodeMap_[0] = 0;  ipNodeMap_[1] = 1;  ipNodeMap_[2] = 3; 
  // face 1; 
  ipNodeMap_[3] = 1;  ipNodeMap_[4] = 2;  ipNodeMap_[5] = 3; 
  // face 2;
  ipNodeMap_[6] = 0;  ipNodeMap_[7] = 3;  ipNodeMap_[8] = 2;  
  // face 3;
  ipNodeMap_[9] = 0; ipNodeMap_[10] = 2; ipNodeMap_[11] = 1;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
TetSCS::~TetSCS()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
TetSCS::ipNodeMap(
  int ordinal)
{
  // define ip->node mappings for each face (ordinal); 
  return &ipNodeMap_[ordinal*3];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void TetSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  SIERRA_FORTRAN(tet_scs_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords, areav );

  // all is always well; no error checking
  *error = 0;
}

//--------------------------------------------------------------------------
//-------- grad_op ---------------------------------------------------------
//--------------------------------------------------------------------------
void TetSCS::grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(tet_derivative)
    ( &numIntPoints_, deriv );
  
  SIERRA_FORTRAN(tet_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );

  if ( lerr )
    std::cout << "sorry, negative TetSCS volume.." << std::endl;  
}

//--------------------------------------------------------------------------
//-------- shifted_grad_op -------------------------------------------------
//--------------------------------------------------------------------------
void TetSCS::shifted_grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(tet_derivative)
    ( &numIntPoints_, deriv );

  SIERRA_FORTRAN(tet_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );

  if ( lerr )
    std::cout << "sorry, negative TetSCS volume.." << std::endl;
}

//--------------------------------------------------------------------------
//-------- face_grad_op ----------------------------------------------------
//--------------------------------------------------------------------------
void TetSCS::face_grad_op(
  const int nelem,
  const int /*face_ordinal*/,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  int lerr = 0;
  int npf = 3;

  const int nface = 1;
  double dpsi[12];
  double grad[12];

  for ( int n=0; n<nelem; n++ ) {

    for ( int k=0; k<npf; k++ ) {

      // derivatives are constant
      SIERRA_FORTRAN(tet_derivative)
        ( &nface, dpsi );

      SIERRA_FORTRAN(tet_gradient_operator)
        ( &nface,
          &nodesPerElement_,
          &nface,
          dpsi,
          &coords[12*n], grad, &det_j[npf*n+k], error, &lerr );

      if ( lerr )
        std::cout << "sorry, issue with face_grad_op.." << std::endl;

      for ( int j=0; j<12; j++) {
        gradop[k*nelem*12+n*12+j] = grad[j];
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- guij ------------------------------------------------------------
//--------------------------------------------------------------------------
void TetSCS::gij(
  const double *coords,
  double *gupperij,
  double *glowerij,
  double *deriv)
{
  SIERRA_FORTRAN(threed_gij)
    ( &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gupperij, glowerij);
}

//--------------------------------------------------------------------------
//-------- adjacentNodes ---------------------------------------------------
//--------------------------------------------------------------------------
const int *
TetSCS::adjacentNodes()
{
  // define L/R mappings
  return &lrscv_[0];
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
TetSCS::shape_fcn(double *shpfc)
{
  tet_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
TetSCS::shifted_shape_fcn(double *shpfc)
{
  tet_shape_fcn(numIntPoints_, &intgLocShift_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- tet_shape_fcn ---------------------------------------------------
//--------------------------------------------------------------------------
void
TetSCS::tet_shape_fcn(
  const int  &npts,
  const double *par_coord, 
  double *shape_fcn)
{
  for (int j = 0; j < npts; ++j ) {
    const int fourj = 4*j;
    const int k = 3*j;
    const double xi = par_coord[k];
    const double eta = par_coord[k+1];
    const double zeta = par_coord[k+2];
    shape_fcn[fourj] = 1.0 - xi - eta - zeta;
    shape_fcn[1 + fourj] = xi;
    shape_fcn[2 + fourj] = eta;
    shape_fcn[3 + fourj] = zeta;
  }
}

//--------------------------------------------------------------------------
//-------- opposingNodes --------------------------------------------------
//--------------------------------------------------------------------------
int
TetSCS::opposingNodes(
  const int ordinal,
  const int node)
{
  return oppNode_[ordinal*3+node];
}

//--------------------------------------------------------------------------
//-------- opposingFace --------------------------------------------------
//--------------------------------------------------------------------------
int
TetSCS::opposingFace(
  const int ordinal,
  const int node)
{
  return oppFace_[ordinal*3+node];
}

//--------------------------------------------------------------------------
//-------- isInElement -----------------------------------------------------
//--------------------------------------------------------------------------
double
TetSCS::isInElement(
    const double * elem_nodal_coor,
    const double * point_coor,
	  double * par_coor ) 
{
  // load up the element coordinates
  const double x1 = elem_nodal_coor[0];
  const double x2 = elem_nodal_coor[1];
  const double x3 = elem_nodal_coor[2];
  const double x4 = elem_nodal_coor[3];

  const double y1 = elem_nodal_coor[4];
  const double y2 = elem_nodal_coor[5];
  const double y3 = elem_nodal_coor[6];
  const double y4 = elem_nodal_coor[7];

  const double z1 = elem_nodal_coor[8];
  const double z2 = elem_nodal_coor[9];
  const double z3 = elem_nodal_coor[10];
  const double z4 = elem_nodal_coor[11];

  // determinant of matrix M in eqn x-x1 = M*xi
  const double det =
    (x2 - x1)*( (y3 - y1)*(z4 - z1) - (y4 - y1)*(z3 - z1) ) -
    (x3 - x1)*( (y2 - y1)*(z4 - z1) - (y4 - y1)*(z2 - z1) ) +
    (x4 - x1)*( (y2 - y1)*(z3 - z1) - (y3 - y1)*(z2 - z1) );

  const double invDet = 1.0/det;

  // matrix entries in inverse of M

  const double m11 = y3*z4 - y1*z4 - y4*z3+y1*z3+y4*z1 - y3*z1;
  const double m12 =  - (x3*z4 - x1*z4 - x4*z3+x1*z3+x4*z1 - x3*z1);
  const double m13 = x3*y4 - x1*y4 - x4*y3+x1*y3+x4*y1 - x3*y1;
  const double m21 =  - (y2*z4 - y1*z4 - y4*z2+y1*z2+y4*z1 - y2*z1);
  const double m22 = x2*z4 - x1*z4 - x4*z2+x1*z2+x4*z1 - x2*z1;
  const double m23 =  - (x2*y4 - x1*y4 - x4*y2+x1*y2+x4*y1 - x2*y1);
  const double m31 = y2*z3 - y1*z3 - y3*z2+y1*z2+y3*z1 - y2*z1;
  const double m32 =  - (x2*z3 - x1*z3 - x3*z2+x1*z2+x3*z1 - x2*z1);
  const double m33 = x2*y3 - x1*y3 - x3*y2+x1*y2+x3*y1 - x2*y1;

  const double xx1 = point_coor[0] - x1;
  const double yy1 = point_coor[1] - y1;
  const double zz1 = point_coor[2] - z1;

  // solve for parametric coordinates

  const double xi   = invDet*(m11*xx1 + m12*yy1 + m13*zz1);
  const double eta  = invDet*(m21*xx1 + m22*yy1 + m23*zz1);
  const double zeta = invDet*(m31*xx1 + m32*yy1 + m33*zz1);

  // if volume coordinates are negative, point is outside the tet

  par_coor[0] = xi;
  par_coor[1] = eta;
  par_coor[2] = zeta;

  std::vector<double> x(3);
  x[0]=par_coor[0];
  x[1]=par_coor[1];
  x[2]=par_coor[2];

  const double dist = parametric_distance(x);

  return dist;
}

//--------------------------------------------------------------------------
//-------- interpolatePoint ------------------------------------------------
//--------------------------------------------------------------------------
void
TetSCS::interpolatePoint(
    const int  & ncomp_field,
    const double * par_coord,           // (3)
    const double * field,               // (4,ncomp_field)
	  double * result ) // (ncomp_field)
{
  const double xi   = par_coord[0];
  const double eta  = par_coord[1];
  const double zeta = par_coord[2];

  const double psi1 = 1.0 - xi - eta - zeta;

  for(int i = 0; i < ncomp_field; ++i) {
    const int fourI = 4*i;

    const double f1 = field[fourI];
    const double f2 = field[1 + fourI];
    const double f3 = field[2 + fourI];
    const double f4 = field[3 + fourI];

    result[i] = f1*psi1 + f2*xi + f3*eta + f4*zeta;
  }
}

//--------------------------------------------------------------------------
//-------- general_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
TetSCS::general_shape_fcn(
  const int numIp,
  const double *isoParCoord,
  double *shpfc)
{
  tet_shape_fcn(numIp, &isoParCoord[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- general_face_grad_op --------------------------------------------
//--------------------------------------------------------------------------
void 
TetSCS::general_face_grad_op(
  const int /*face_ordinal*/,
  const double */*isoParCoord*/,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  int lerr = 0;

  const int nface = 1;
  double dpsi[12];
  double grad[12];

  // derivatives are constant
  SIERRA_FORTRAN(tet_derivative)
    ( &nface, dpsi );

  SIERRA_FORTRAN(tet_gradient_operator)
    ( &nface,
      &nodesPerElement_,
      &nface,
      dpsi,
      &coords[0], grad, &det_j[0], error, &lerr );
  
  if ( lerr )
    throw std::runtime_error("TetSCS::general_face_grad_op issue");
 
  for ( int j=0; j<12; j++) {
    gradop[j] = grad[j];
  }
}

//--------------------------------------------------------------------------
//-------- sidePcoords_to_elemPcoords --------------------------------------
//--------------------------------------------------------------------------
void 
TetSCS::sidePcoords_to_elemPcoords(
  const int & side_ordinal,
  const int & npoints,
  const double *side_pcoords,
  double *elem_pcoords)
{
  switch (side_ordinal) {
  case 0:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = side_pcoords[2*i+0];
      elem_pcoords[i*3+1] = 0.0;
      elem_pcoords[i*3+2] = side_pcoords[2*i+1];
    }
    break;
  case 1:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = 1.0 - side_pcoords[2*i+0] - side_pcoords[2*i+1];
      elem_pcoords[i*3+1] = side_pcoords[2*i+0];
      elem_pcoords[i*3+2] = side_pcoords[2*i+1];
    }
    break;
  case 2:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = 0.0;
      elem_pcoords[i*3+1] = side_pcoords[2*i+1];
      elem_pcoords[i*3+2] = side_pcoords[2*i+0];
    }
    break;
  case 3:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*3+0] = side_pcoords[2*i+1];
      elem_pcoords[i*3+1] = side_pcoords[2*i+0];
      elem_pcoords[i*3+2] = 0.0;
    }
    break;
  default:
    throw std::runtime_error("TetSCS::sideMap invalid ordinal");
  }
  return;
}

//--------------------------------------------------------------------------
//-------- parametric_distance ---------------------------------------------
//--------------------------------------------------------------------------
double
TetSCS::parametric_distance(const std::vector<double> &x)
{
  const double X=x[0] - 1./4.;
  const double Y=x[1] - 1./4.;
  const double Z=x[2] - 1./4.;
  const double dist0 = -4*X;
  const double dist1 = -4*Y;
  const double dist2 = -4*Z;
  const double dist3 =  4*(X+Y+Z);
  const double dist  = std::max(std::max(dist0,dist1),std::max(dist2,dist3));
  return dist;
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
PyrSCV::PyrSCV()
  : MasterElement()
{
  nDim_ = 3;
  nodesPerElement_ = 5;
  numIntPoints_ = 5; 

  // define ip node mappings
  ipNodeMap_.resize(5);
  ipNodeMap_[0] = 0; ipNodeMap_[1] = 1; ipNodeMap_[2] = 2; ipNodeMap_[3] = 3;
  ipNodeMap_[4] = 4;

  // standard integration location
  intgLoc_.resize(15);
  intgLoc_[0]  = -19.0/48.0; intgLoc_[1]  = -19.0/48.0; intgLoc_[2]  = 41.0/240.0;  // vol 0
  intgLoc_[3]  =  19.0/48.0; intgLoc_[4]  = -19.0/48.0; intgLoc_[5]  = 41.0/240.0;  // vol 1
  intgLoc_[6]  =  19.0/48.0; intgLoc_[7]  =  19.0/48.0; intgLoc_[8]  = 41.0/240.0;  // vol 2
  intgLoc_[9]  = -19.0/48.0; intgLoc_[10] =  19.0/48.0; intgLoc_[11] = 41.0/240.0;  // vol 3
  intgLoc_[12] =   0.0;      intgLoc_[13] =   0.0;      intgLoc_[14] = 0.6 ;        // vol 4

  // shifted
  intgLocShift_.resize(15);
  intgLocShift_[0]  = -1.0; intgLocShift_[1]  = -1.0; intgLocShift_[2]  = 0.0;  // vol 0
  intgLocShift_[3]  =  1.0; intgLocShift_[4]  = -1.0; intgLocShift_[5]  = 0.0;  // vol 1
  intgLocShift_[6]  =  1.0; intgLocShift_[7]  =  1.0; intgLocShift_[8]  = 0.0;  // vol 2
  intgLocShift_[9]  = -1.0; intgLocShift_[10] =  1.0; intgLocShift_[11] = 0.0;  // vol 3
  intgLocShift_[12] =  0.0; intgLocShift_[13] =  0.0; intgLocShift_[14] = 1.0;  // vol 4
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
PyrSCV::~PyrSCV()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
PyrSCV::ipNodeMap(
  int /*ordinal*/)
{
  // define scv->node mappings
  return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void PyrSCV::determinant(
  const int nelem,
  const double *coords,
  double *volume,
  double *error)
{

  int lerr = 0;

  SIERRA_FORTRAN(pyr_scv_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords,
      volume, error, &lerr );
}


//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
PyrSCV::shape_fcn(double *shpfc)
{
  pyr_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
PyrSCV::shifted_shape_fcn(double *shpfc)
{
  pyr_shape_fcn(numIntPoints_, &intgLocShift_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- pyr_shape_fcn ---------------------------------------------------
//--------------------------------------------------------------------------
void
PyrSCV::pyr_shape_fcn(
  const int  &npts,
  const double *par_coord, 
  double *shape_fcn)
{
  const double one  = 1.0;
  for ( int j = 0; j < npts; ++j ) {
    const int fivej = 5*j;
    const int k     = 3*j;
    const double r    = par_coord[k+0];
    const double s    = par_coord[k+1];
    const double t    = par_coord[k+2];

    shape_fcn[0 + fivej] = 0.25*(1.0-r)*(1.0-s)*(one-t);
    shape_fcn[1 + fivej] = 0.25*(1.0+r)*(1.0-s)*(one-t);
    shape_fcn[2 + fivej] = 0.25*(1.0+r)*(1.0+s)*(one-t);
    shape_fcn[3 + fivej] = 0.25*(1.0-r)*(1.0+s)*(one-t);
    shape_fcn[4 + fivej] = t;
  }
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
PyrSCS::PyrSCS()
  : MasterElement()
{
  nDim_ = 3;
  nodesPerElement_ = 5;
  numIntPoints_ = 8;

  // define L/R mappings
  lrscv_.resize(16);
  lrscv_[0]  = 0; lrscv_[1]  = 1;
  lrscv_[2]  = 1; lrscv_[3]  = 2;
  lrscv_[4]  = 2; lrscv_[5]  = 3;
  lrscv_[6]  = 0; lrscv_[7]  = 3;
  lrscv_[8]  = 0; lrscv_[9]  = 4;
  lrscv_[10] = 1; lrscv_[11] = 4;
  lrscv_[12] = 2; lrscv_[13] = 4;
  lrscv_[14] = 3; lrscv_[15] = 4;
  
  // define opposing node
  oppNode_.resize(20);
  // face 0; nodes 0,1,4
  oppNode_[0] = 3; oppNode_[1] = 2; oppNode_[2] = 4; oppNode_[3] = -1;
  // face 1; nodes 1,2,4
  oppNode_[4] = 0; oppNode_[5] = 3; oppNode_[6] = 4; oppNode_[7] = -1;
  // face 2; nodes 2,3,4
  oppNode_[8] = 1; oppNode_[9] = 4; oppNode_[10] = 2; oppNode_[11] = -1;
  // face 3; nodes 0,4,3
  oppNode_[12] = 1; oppNode_[13] = 5; oppNode_[14] = 6; oppNode_[15] = 2;
  // face 4; nodes 1,3,2,1
  oppNode_[16] = 4; oppNode_[17] = 4; oppNode_[18] = 4; oppNode_[19] = 4;

  // define opposing face

  // standard integration location
  intgLoc_.resize(24);
  const double fortyFiveHundredFourths = 45.0/104.0;
  const double fortyOneHundredTwentyths = 41.0/120.0;
  const double sevenFiftyTwoths = 7.0/52.0;
  const double sevenTwentyFourths = 7.0/24.0;
  intgLoc_[0]  =  0.00;                     intgLoc_[1]  = -fortyFiveHundredFourths; intgLoc_[2]  = sevenFiftyTwoths; // surf 1    1->2
  intgLoc_[3]  =  fortyFiveHundredFourths;  intgLoc_[4]  = 0.00;                     intgLoc_[5]  = sevenFiftyTwoths; // surf 2    2->3
  intgLoc_[6]  =  0.00;                     intgLoc_[7]  = fortyFiveHundredFourths;  intgLoc_[8]  = sevenFiftyTwoths; // surf 3    3->4
  intgLoc_[9]  =  -fortyFiveHundredFourths; intgLoc_[10] = 0.0;                      intgLoc_[11] = sevenFiftyTwoths; // surf 4    1->4
  intgLoc_[12] =  -sevenTwentyFourths;      intgLoc_[13] = -sevenTwentyFourths;      intgLoc_[14] = fortyOneHundredTwentyths; // surf 5    1->5
  intgLoc_[15] =  sevenTwentyFourths;       intgLoc_[16] = -sevenTwentyFourths;      intgLoc_[17] = fortyOneHundredTwentyths; // surf 6    2->5
  intgLoc_[18] =  sevenTwentyFourths;       intgLoc_[19] = sevenTwentyFourths;       intgLoc_[20] = fortyOneHundredTwentyths; // surf 7    3->5
  intgLoc_[21] =  -sevenTwentyFourths;      intgLoc_[22] = sevenTwentyFourths;       intgLoc_[23] = fortyOneHundredTwentyths; // surf 8    4->5

  // shifted
  intgLocShift_.resize(24);
  intgLocShift_[0]  =  0.00; intgLocShift_[1]  = -1.00; intgLocShift_[2]  =  0.00; // surf 1    1->2
  intgLocShift_[3]  =  1.00; intgLocShift_[4]  =  0.00; intgLocShift_[5]  =  0.00; // surf 2    2->3
  intgLocShift_[6]  =  0.00; intgLocShift_[7]  =  1.00; intgLocShift_[8]  =  0.00; // surf 3    3->4
  intgLocShift_[9]  = -1.00; intgLocShift_[10] =  0.00; intgLocShift_[11] =  0.00; // surf 4    1->4
  intgLocShift_[12] = -0.50; intgLocShift_[13] = -0.50; intgLocShift_[14] =  0.50; // surf 5    1->5
  intgLocShift_[15] =  0.50; intgLocShift_[16] = -0.50; intgLocShift_[17] =  0.50; // surf 6    2->5
  intgLocShift_[18] =  0.50; intgLocShift_[19] =  0.50; intgLocShift_[20] =  0.50; // surf 7    3->5
  intgLocShift_[21] = -0.50; intgLocShift_[22] =  0.50; intgLocShift_[23] =  0.50; // surf 8    4->5

  // exposed face; n/a
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
PyrSCS::~PyrSCS()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void PyrSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  SIERRA_FORTRAN(pyr_scs_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords, areav );

  // all is always well; no error checking
  *error = 0;
}

//--------------------------------------------------------------------------
//-------- grad_op ---------------------------------------------------------
//--------------------------------------------------------------------------
void PyrSCS::grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  pyr_derivative(numIntPoints_, &intgLoc_[0], deriv);
  
  SIERRA_FORTRAN(pyr_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );

  if ( lerr )
    std::cout << "sorry, negative PyrSCS volume.." << std::endl;
}

//--------------------------------------------------------------------------
//-------- shifted_grad_op -------------------------------------------------
//--------------------------------------------------------------------------
void PyrSCS::shifted_grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  pyr_derivative(numIntPoints_, &intgLocShift_[0], deriv);

  SIERRA_FORTRAN(pyr_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );

  if ( lerr )
    std::cout << "sorry, negative PyrSCS volume.." << std::endl;
}

//--------------------------------------------------------------------------
//-------- pyr_derivative --------------------------------------------------
//--------------------------------------------------------------------------
void PyrSCS::pyr_derivative(
  const int npts,
  const double *intgLoc,
  double *deriv)
{
  // d3d(c,s,j) = deriv[c + 3*(s + 5*j)] = deriv[c+3s+15j]

  for ( int j = 0; j < npts; ++j) {
    const int k = j*3;
    const int p = 15*j;

    double r = intgLoc[k+0];
    double s = intgLoc[k+1];
    double t = intgLoc[k+2];

    deriv[0+3*0+p] =-0.25*(1.0-s)*(1.0-t);  // d(N_1)/ d(r) = deriv[0]
    deriv[1+3*0+p] =-0.25*(1.0-r)*(1.0-t);  // d(N_1)/ d(s) = deriv[1]
    deriv[2+3*0+p] =-0.25*(1.0-r)*(1.0-s);  // d(N_1)/ d(t) = deriv[2]

    deriv[0+3*1+p] = 0.25*(1.0-s)*(1.0-t);  // d(N_2)/ d(r) = deriv[0+3]
    deriv[1+3*1+p] =-0.25*(1.0+r)*(1.0-t);  // d(N_2)/ d(s) = deriv[1+3]
    deriv[2+3*1+p] =-0.25*(1.0+r)*(1.0-s);  // d(N_2)/ d(t) = deriv[2+3]

    deriv[0+3*2+p] = 0.25*(1.0+s)*(1.0-t);  // d(N_3)/ d(r) = deriv[0+6]
    deriv[1+3*2+p] = 0.25*(1.0+r)*(1.0-t);  // d(N_3)/ d(s) = deriv[1+6]
    deriv[2+3*2+p] =-0.25*(1.0+r)*(1.0+s);  // d(N_3)/ d(t) = deriv[2+6]

    deriv[0+3*3+p] =-0.25*(1.0+s)*(1.0-t);  // d(N_4)/ d(r) = deriv[0+9]
    deriv[1+3*3+p] = 0.25*(1.0-r)*(1.0-t);  // d(N_4)/ d(s) = deriv[1+9]
    deriv[2+3*3+p] =-0.25*(1.0-r)*(1.0+s);  // d(N_4)/ d(t) = deriv[2+9]

    deriv[0+3*4+p] = 0.0;                   // d(N_5)/ d(r) = deriv[0+12]
    deriv[1+3*4+p] = 0.0;                   // d(N_5)/ d(s) = deriv[1+12]
    deriv[2+3*4+p] = 1.0;                   // d(N_5)/ d(t) = deriv[2+12]
  }
}

//--------------------------------------------------------------------------
//-------- guij ------------------------------------------------------------
//--------------------------------------------------------------------------
void PyrSCS::gij(
  const double *coords,
  double *gupperij,
  double *glowerij,
  double *deriv)
{
  SIERRA_FORTRAN(threed_gij)
    ( &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gupperij, glowerij);
}

//--------------------------------------------------------------------------
//-------- adjacentNodes ---------------------------------------------------
//--------------------------------------------------------------------------
const int *
PyrSCS::adjacentNodes()
{
  // define L/R mappings
  return &lrscv_[0];
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
PyrSCS::shape_fcn(double *shpfc)
{
  pyr_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
PyrSCS::shifted_shape_fcn(double *shpfc)
{
  pyr_shape_fcn(numIntPoints_, &intgLocShift_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- pyr_shape_fcn ---------------------------------------------------
//--------------------------------------------------------------------------
void
PyrSCS::pyr_shape_fcn(
  const int  &npts,
  const double *par_coord, 
  double *shape_fcn)
{
  const double one  = 1.0;
  for ( int j = 0; j < npts; ++j ) {
    const int fivej = 5*j;
    const int k     = 3*j;
    const double r    = par_coord[k+0];
    const double s    = par_coord[k+1];
    const double t    = par_coord[k+2];

    shape_fcn[0 + fivej] = 0.25*(1.0-r)*(1.0-s)*(one-t);
    shape_fcn[1 + fivej] = 0.25*(1.0+r)*(1.0-s)*(one-t);
    shape_fcn[2 + fivej] = 0.25*(1.0+r)*(1.0+s)*(one-t);
    shape_fcn[3 + fivej] = 0.25*(1.0-r)*(1.0+s)*(one-t);
    shape_fcn[4 + fivej] = t;
  }
}

//--------------------------------------------------------------------------
//-------- opposingNodes --------------------------------------------------
//--------------------------------------------------------------------------
int
PyrSCS::opposingNodes(
  const int ordinal,
  const int node)
{
  return oppNode_[ordinal*4+node];
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
WedSCV::WedSCV()
  : MasterElement()
{
  nDim_ = 3;
  nodesPerElement_ = 6;
  numIntPoints_ = 6;

  // define ip node mappings
  ipNodeMap_.resize(6);
  ipNodeMap_[0] = 0; ipNodeMap_[1] = 1; ipNodeMap_[2] = 2;
  ipNodeMap_[3] = 3; ipNodeMap_[4] = 4; ipNodeMap_[5] = 5;

  // standard integration location
  intgLoc_.resize(18);
  const double seven12th = 7.0/12.0;
  const double five24th = 5.0/24.0;
  intgLoc_[0]  = five24th;  intgLoc_[1]  = five24th;  intgLoc_[2]  = -0.5; // vol 0
  intgLoc_[3]  = seven12th; intgLoc_[4]  = five24th;  intgLoc_[5]  = -0.5; // vol 1
  intgLoc_[6]  = five24th;  intgLoc_[7]  = seven12th; intgLoc_[8]  = -0.5; // vol 2
  intgLoc_[9]  = five24th;  intgLoc_[10] = five24th;  intgLoc_[11] = 0.5;  // vol 3
  intgLoc_[12] = seven12th; intgLoc_[13] = five24th;  intgLoc_[14] = 0.5;  // vol 4
  intgLoc_[15] = five24th;  intgLoc_[16] = seven12th; intgLoc_[17] = 0.5;  // vol 5

  // shifted
  intgLocShift_.resize(18);
  intgLocShift_[0]  = 0.0;  intgLocShift_[1]  = 0.0; intgLocShift_[2]  = -1.0; // vol 0
  intgLocShift_[3]  = 1.0;  intgLocShift_[4]  = 0.0; intgLocShift_[5]  = -1.0; // vol 1
  intgLocShift_[6]  = 0.0;  intgLocShift_[7]  = 1.0; intgLocShift_[8]  = -1.0; // vol 2
  intgLocShift_[9]  = 0.0;  intgLocShift_[10] = 0.0; intgLocShift_[11] =  1.0; // vol 3
  intgLocShift_[12] = 1.0;  intgLocShift_[13] = 0.0; intgLocShift_[14] = 1.0;  // vol 4
  intgLocShift_[15] = 0.0;  intgLocShift_[16] = 1.0; intgLocShift_[17] = 1.0;  // vol 5
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
WedSCV::~WedSCV()
{
  // does nothing
}
//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
WedSCV::ipNodeMap(
  int /*ordinal*/)
{
  // define scv->node mappings
  return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void WedSCV::determinant(
  const int nelem,
  const double *coords,
  double *volume,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(wed_scv_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords,
      volume, error, &lerr );
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
WedSCV::shape_fcn(double *shpfc)
{
  wedge_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
WedSCV::shifted_shape_fcn(double *shpfc)
{
  wedge_shape_fcn(numIntPoints_, &intgLocShift_[0], shpfc);
}


//--------------------------------------------------------------------------
//-------- wedge_shape_fcn -------------------------------------------------
//--------------------------------------------------------------------------
void
WedSCV::wedge_shape_fcn(
  const int  &npts,
  const double *isoParCoord, 
  double *shape_fcn)
{
  for (int j = 0; j < npts; ++j ) {
    int sixj = 6 * j;
    int k    = 3 * j;
    double r    = isoParCoord[k];
    double s    = isoParCoord[k + 1];
    double t    = 1.0 - r - s;
    double xi   = isoParCoord[k + 2];
    shape_fcn[    sixj] = 0.5 * t * (1.0 - xi);
    shape_fcn[1 + sixj] = 0.5 * r * (1.0 - xi);
    shape_fcn[2 + sixj] = 0.5 * s * (1.0 - xi);
    shape_fcn[3 + sixj] = 0.5 * t * (1.0 + xi);
    shape_fcn[4 + sixj] = 0.5 * r * (1.0 + xi);
    shape_fcn[5 + sixj] = 0.5 * s * (1.0 + xi);
  }
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
WedSCS::WedSCS()
  : MasterElement()
{
  nDim_ = 3;
  nodesPerElement_ = 6;
  numIntPoints_ = 9;

  // define L/R mappings
  lrscv_.resize(18);
  lrscv_[0]  = 0; lrscv_[1]  = 1;
  lrscv_[2]  = 1; lrscv_[3]  = 2;
  lrscv_[4]  = 0; lrscv_[5]  = 2;
  lrscv_[6]  = 3; lrscv_[7]  = 4;
  lrscv_[8]  = 4; lrscv_[9]  = 5;
  lrscv_[10] = 3; lrscv_[11] = 5;
  lrscv_[12] = 0; lrscv_[13] = 3;
  lrscv_[14] = 1; lrscv_[15] = 4;
  lrscv_[16] = 2; lrscv_[17] = 5;

  // define opposing node
  oppNode_.resize(20);
  // face 0; nodes 0,1,4,3
  oppNode_[0] = 2; oppNode_[1] = 2; oppNode_[2] = 5; oppNode_[3] = 5;
  // face 1; nodes 1,2,5,4
  oppNode_[4] = 0; oppNode_[5] = 0; oppNode_[6] = 3; oppNode_[7] = 3;
  // face 2; nodes 0,3,5,2
  oppNode_[8] = 1; oppNode_[9] = 4; oppNode_[10] = 4; oppNode_[11] = 1;
  // face 3; nodes 0,2,1
  oppNode_[12] = 3; oppNode_[13] = 5; oppNode_[14] = 4; oppNode_[15] = -1;
  // face 4; nodes 3,4,5
  oppNode_[16] = 0; oppNode_[17] = 1; oppNode_[18] = 2; oppNode_[19] = -1;

  // define opposing face
  oppFace_.resize(20);
  // face 0; nodes 0, 1, 4, 3
  oppFace_[0]  = 2; oppFace_[1]  = 1; oppFace_[2]  = 4;  oppFace_[3]  = 5;
  // face 1; nodes 1,2,5,4
  oppFace_[4]  = 0; oppFace_[5]  = 2; oppFace_[6]  = 5;  oppFace_[7]  = 3;
  // face 2, nodes 0,3,5,2
  oppFace_[8]  = 0; oppFace_[9]  = 3; oppFace_[10] = 4;  oppFace_[11] = 1;
  // face 3, nodes 0,2,1
  oppFace_[12] = 6; oppFace_[13] = 8; oppFace_[14] = 7;  oppFace_[15] = -1;
  //face 4, nodes 3,4,5
  oppFace_[16] = 6; oppFace_[17] = 7; oppFace_[18] = 8;  oppFace_[19] = -1;
 
  // standard integration location
  const double oneSixth = 1.0/6.0;
  const double five12th = 5.0/12.0;
  const double seven12th = 7.0/12.0;
  const double five24th = 5.0/24.0;
  intgLoc_.resize(27);    
  intgLoc_[0]  =  five12th;  intgLoc_[1]  = oneSixth;  intgLoc_[2]  = -0.50; // surf 1    1->2
  intgLoc_[3]  =  five12th;  intgLoc_[4]  = five12th;  intgLoc_[5]  = -0.50; // surf 2    2->3
  intgLoc_[6]  =  oneSixth;  intgLoc_[7]  = five12th;  intgLoc_[8]  = -0.50; // surf 3    1->3
  intgLoc_[9]  =  five12th;  intgLoc_[10] = oneSixth;  intgLoc_[11] =  0.50; // surf 4    4->5
  intgLoc_[12] =  five12th;  intgLoc_[13] = five12th;  intgLoc_[14] =  0.50; // surf 5    5->6
  intgLoc_[15] =  oneSixth;  intgLoc_[16] = five12th;  intgLoc_[17] =  0.50; // surf 6    4->6
  intgLoc_[18] =  five24th;  intgLoc_[19] = five24th;  intgLoc_[20] =  0.00; // surf 7    1->4
  intgLoc_[21] =  seven12th; intgLoc_[22] = five24th;  intgLoc_[23] =  0.00; // surf 8    2->5
  intgLoc_[24] =  five24th;  intgLoc_[25] = seven12th; intgLoc_[26] =  0.00; // surf 9    3->6

  // shifted
  intgLocShift_.resize(27);
  intgLocShift_[0]  =  0.50; intgLocShift_[1]  =  0.00; intgLocShift_[2]  = -1.00; // surf 1    1->2
  intgLocShift_[3]  =  0.50; intgLocShift_[4]  =  0.50; intgLocShift_[5]  = -1.00; // surf 2    2->3
  intgLocShift_[6]  =  0.00; intgLocShift_[7]  =  0.50; intgLocShift_[8]  = -1.00; // surf 3    1->3
  intgLocShift_[9]  =  0.50; intgLocShift_[10] =  0.00; intgLocShift_[11] =  1.00; // surf 4    4->5
  intgLocShift_[12] =  0.50; intgLocShift_[13] =  0.50; intgLocShift_[14] =  1.00; // surf 5    5->6
  intgLocShift_[15] =  0.00; intgLocShift_[16] =  0.50; intgLocShift_[17] =  1.00; // surf 6    4->6
  intgLocShift_[18] =  0.00; intgLocShift_[19] =  0.00; intgLocShift_[20] =  0.00; // surf 7    1->4
  intgLocShift_[21] =  1.00; intgLocShift_[22] =  0.00; intgLocShift_[23] =  0.00; // surf 8    2->5
  intgLocShift_[24] =  0.00; intgLocShift_[25] =  1.00; intgLocShift_[26] =  0.00; // surf 9    3->6

  // exposed face
  intgExpFace_.resize(60);
  intgExpFace_[0] = 0.25;       intgExpFace_[1]  = 0.0;       intgExpFace_[2] = -0.5;  // surf 0; nodes 0,1,4,3
  intgExpFace_[3] = 0.75;       intgExpFace_[4]  = 0.0;       intgExpFace_[5] = -0.5;  // face 0, surf 1
  intgExpFace_[6] = 0.75;       intgExpFace_[7]  = 0.0;       intgExpFace_[8] =  0.5;  // face 0, surf 2
  intgExpFace_[9] = 0.25;       intgExpFace_[10] = 0.0;       intgExpFace_[11] = 0.5;  // face 0, surf 3
  intgExpFace_[12] = 0.75;      intgExpFace_[13] = 0.25;      intgExpFace_[14] = -0.5; // surf 1; nodes 1,2,5,4
  intgExpFace_[15] = 0.25;      intgExpFace_[16] = 0.75;      intgExpFace_[17] = -0.5; // face 1, surf 1
  intgExpFace_[18] = 0.25;      intgExpFace_[19] = 0.75;      intgExpFace_[20] =  0.5; // face 1, surf 2
  intgExpFace_[21] = 0.75;      intgExpFace_[22] = 0.25;      intgExpFace_[23] =  0.5; // face 1, surf 3
  intgExpFace_[24] = 0.0;       intgExpFace_[25] = 0.25;      intgExpFace_[26] = -0.5; // surf 2; nodes 0,3,5,2
  intgExpFace_[27] = 0.0;       intgExpFace_[28] = 0.25;      intgExpFace_[29] =  0.5; // face 2, surf 1
  intgExpFace_[30] = 0.0;       intgExpFace_[31] = 0.75;      intgExpFace_[32] =  0.5; // face 2, surf 2
  intgExpFace_[33] = 0.0;       intgExpFace_[34] = 0.75;      intgExpFace_[35] = -0.5; // face 2, surf 3
  intgExpFace_[36] = five24th;  intgExpFace_[37] = five24th;  intgExpFace_[38] = -1.0; // surf 3; nodes 0,2,1
  intgExpFace_[39] = five24th;  intgExpFace_[40] = seven12th; intgExpFace_[41] = -1.0; // face 3, surf 1
  intgExpFace_[42] = seven12th; intgExpFace_[43] = five24th;  intgExpFace_[44] = -1.0; // face 3, surf 2
  intgExpFace_[45] = 0.0;       intgExpFace_[46] = 0.0;       intgExpFace_[47] =  0.0; // (blank)
  intgExpFace_[48] = five24th;  intgExpFace_[49] = five24th;  intgExpFace_[50] = 1.0;  // surf 4; nodes 3,4,5
  intgExpFace_[51] = seven12th; intgExpFace_[52] = five24th;  intgExpFace_[53] = 1.0;  // face 4, surf 1
  intgExpFace_[54] = five24th;  intgExpFace_[55] = seven12th; intgExpFace_[56] = 1.0;  // face 4, surf 2
  intgExpFace_[57] = 0.0;       intgExpFace_[58] = 0.0;       intgExpFace_[59] = 0.0;  // (blank)

  // boundary integration point ip node mapping (ip on an ordinal to local node number)
  ipNodeMap_.resize(20); // 4 ips (pick quad) * 5 faces
  // face 0;
  ipNodeMap_[0] = 0;  ipNodeMap_[1] = 1;  ipNodeMap_[2] = 4;  ipNodeMap_[3] = 3; 
  // face 1; 
  ipNodeMap_[4] = 1;  ipNodeMap_[5] = 2;  ipNodeMap_[6] = 5;  ipNodeMap_[7] = 4; 
  // face 2;
  ipNodeMap_[8] = 0;  ipNodeMap_[9] = 3;  ipNodeMap_[10] = 5; ipNodeMap_[11] = 2; 
  // face 3;
  ipNodeMap_[12] = 0; ipNodeMap_[13] = 1; ipNodeMap_[14] = 1; ipNodeMap_[15] = 0; //empty 
  // face 4;
  ipNodeMap_[16] = 3; ipNodeMap_[17] = 4; ipNodeMap_[18] = 5; ipNodeMap_[19] = 0; // empty 
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
WedSCS::~WedSCS()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
WedSCS::ipNodeMap(
  int ordinal)
{
  // define ip->node mappings for each face (ordinal); 
  return &ipNodeMap_[ordinal*4];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void WedSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  SIERRA_FORTRAN(wed_scs_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords, areav );

  // all is always well; no error checking
  *error = 0;
}

//--------------------------------------------------------------------------
//-------- grad_op ---------------------------------------------------------
//--------------------------------------------------------------------------
void WedSCS::grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  wedge_derivative(numIntPoints_, &intgLoc_[0], deriv);

  SIERRA_FORTRAN(wed_gradient_operator) (
      &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );

  if ( lerr )
    std::cout << "sorry, negative WedSCS volume.." << std::endl;
}

//--------------------------------------------------------------------------
//-------- shifted_grad_op -------------------------------------------------
//--------------------------------------------------------------------------
void WedSCS::shifted_grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  wedge_derivative(numIntPoints_, &intgLocShift_[0], deriv);

  SIERRA_FORTRAN(wed_gradient_operator) (
      &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );

  if ( lerr )
    std::cout << "sorry, negative WedSCS volume.." << std::endl;
}

//--------------------------------------------------------------------------
//-------- wedge_derivative --------------------------------------------------
//--------------------------------------------------------------------------
void WedSCS::wedge_derivative(
  const int npts,
  const double *intgLoc,
  double *deriv)
{
  // d3d(c,s,j) = deriv[c + 3*(s + 6*j)] = deriv[c+3s+18j]

  for (int  j = 0; j < npts; ++j) {

    int k  = j*3;
    const int p = 18*j;

    const double r  = intgLoc[k];
    const double s  = intgLoc[k+1];
    const double t  = 1.0 - r - s;
    const double xi = intgLoc[k + 2];

    deriv[0+3*0+p] = -0.5 * (1.0 - xi);  // d(N_1)/ d(r)  = deriv[0]
    deriv[1+3*0+p] = -0.5 * (1.0 - xi);  // d(N_1)/ d(s)  = deriv[1]
    deriv[2+3*0+p] = -0.5 * t;           // d(N_1)/ d(xi) = deriv[2]

    deriv[0+3*1+p] =  0.5 * (1.0 - xi);  // d(N_2)/ d(r)  = deriv[0 + 3]
    deriv[1+3*1+p] =  0.0;               // d(N_2)/ d(s)  = deriv[1 + 3]
    deriv[2+3*1+p] = -0.5 * r;           // d(N_2)/ d(xi) = deriv[2 + 3]

    deriv[0+3*2+p] =  0.0;               // d(N_3)/ d(r)  = deriv[0 + 6]
    deriv[1+3*2+p] =  0.5 * (1.0 - xi);  // d(N_3)/ d(s)  = deriv[1 + 6]
    deriv[2+3*2+p] = -0.5 * s;           // d(N_3)/ d(xi) = deriv[2 + 6]

    deriv[0+3*3+p] = -0.5 * (1.0 + xi);  // d(N_4)/ d(r)  = deriv[0 + 9]
    deriv[1+3*3+p] = -0.5 * (1.0 + xi);  // d(N_4)/ d(s)  = deriv[1 + 9]
    deriv[2+3*3+p] =  0.5 * t;           // d(N_4)/ d(xi) = deriv[2 + 9]

    deriv[0+3*4+p] =  0.5 * (1.0 + xi);  // d(N_5)/ d(r)  = deriv[0 + 12]
    deriv[1+3*4+p] =  0.0;               // d(N_5)/ d(s)  = deriv[1 + 12]
    deriv[2+3*4+p] =  0.5 * r;           // d(N_5)/ d(xi) = deriv[2 + 12]

    deriv[0+3*5+p] =  0.0;               // d(N_6)/ d(r)  = deriv[0 + 15]
    deriv[1+3*5+p] =  0.5 * (1.0 + xi);  // d(N_6)/ d(s)  = deriv[1 + 15]
    deriv[2+3*5+p] =  0.5 * s;           // d(N_6)/ d(xi) = deriv[2 + 15]
  }
}

//--------------------------------------------------------------------------
//-------- face_grad_op ----------------------------------------------------
//--------------------------------------------------------------------------
void 
WedSCS::face_grad_op(
  const int nelem,
  const int face_ordinal,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  int lerr = 0;
  const int ndim = 3;
  const int nface = 1;
  double dpsi[18], grad[18];

  // nodes per face... ordinal 0, 1, 2 are quad faces, 3 and 4 are tri
  const int npf = (face_ordinal < 3 ) ? 4 : 3;

  for ( int n = 0; n < nelem; ++n ) {

    for ( int k=0; k<npf; ++k ) {
      
      const int row = 12*face_ordinal +k*ndim;
      wedge_derivative(nface, &intgExpFace_[row], dpsi);
      
      SIERRA_FORTRAN(wed_gradient_operator) (
          &nface,
          &nodesPerElement_,
          &nface,
          dpsi,
          &coords[18*n], grad, &det_j[npf*n+k], error, &lerr );
      
      if ( lerr )
        std::cout << "problem with EwedSCS::face_grad" << std::endl;
      
      for ( int j=0; j<18; j++) {
        gradop[k*nelem*18+n*18+j] = grad[j];
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- guij ------------------------------------------------------------
//--------------------------------------------------------------------------
void WedSCS::gij(
  const double *coords,
  double *gupperij,
  double *glowerij,
  double *deriv)
{
  SIERRA_FORTRAN(threed_gij)
    ( &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gupperij, glowerij);
}

//--------------------------------------------------------------------------
//-------- adjacentNodes ---------------------------------------------------
//--------------------------------------------------------------------------
const int *
WedSCS::adjacentNodes()
{
  // define L/R mappings
  return &lrscv_[0];
}

//--------------------------------------------------------------------------
//-------- opposingNodes --------------------------------------------------
//--------------------------------------------------------------------------
int
WedSCS::opposingNodes(
  const int ordinal,
  const int node)
{
  return oppNode_[ordinal*4+node];
}


//--------------------------------------------------------------------------
//-------- opposingFace --------------------------------------------------
//--------------------------------------------------------------------------
int
WedSCS::opposingFace(
  const int ordinal,
  const int node)
{
  return oppFace_[ordinal*4+node];
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
WedSCS::shape_fcn(double *shpfc)
{
  wedge_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
WedSCS::shifted_shape_fcn(double *shpfc)
{
  wedge_shape_fcn(numIntPoints_, &intgLocShift_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- isInElement -----------------------------------------------------
//--------------------------------------------------------------------------
double
WedSCS::isInElement(
  const double *elemNodalCoord,
  const double *pointCoord,
  double *isoParCoord )
{
  const double isInElemConverged = 1.0e-16;

  // ------------------------------------------------------------------
  // Pentahedron master element space is (r,s,xi):
  // r=([0,1]), s=([0,1]), xi=([-1,+1])
  // Use natural coordinates to determine if point is in pentahedron.
  // ------------------------------------------------------------------

  // Translate element so that (x,y,z) coordinates of first node are (0,0,0)

  double x[] = {0.0,
                elemNodalCoord[ 1] - elemNodalCoord[ 0],
                elemNodalCoord[ 2] - elemNodalCoord[ 0],
                elemNodalCoord[ 3] - elemNodalCoord[ 0],
                elemNodalCoord[ 4] - elemNodalCoord[ 0],
                elemNodalCoord[ 5] - elemNodalCoord[ 0] };
  double y[] = {0.0,
                elemNodalCoord[ 7] - elemNodalCoord[ 6],
                elemNodalCoord[ 8] - elemNodalCoord[ 6],
                elemNodalCoord[ 9] - elemNodalCoord[ 6],
                elemNodalCoord[10] - elemNodalCoord[ 6],
                elemNodalCoord[11] - elemNodalCoord[ 6] };
  double z[] = {0.0,
                elemNodalCoord[13] - elemNodalCoord[12],
                elemNodalCoord[14] - elemNodalCoord[12],
                elemNodalCoord[15] - elemNodalCoord[12],
                elemNodalCoord[16] - elemNodalCoord[12],
                elemNodalCoord[17] - elemNodalCoord[12] };

  // (xp,yp,zp) is the point to be mapped into (r,s,xi) coordinate system.
  // This point must also be translated as above.

  double xp = pointCoord[0] - elemNodalCoord[ 0];
  double yp = pointCoord[1] - elemNodalCoord[ 6];
  double zp = pointCoord[2] - elemNodalCoord[12];

  // Newton-Raphson iteration for (r,s,xi)
  double j[3][3];
  double jinv[3][3];
  double f[3];
  double shapefct[6];
  double rnew   = 1.0 / 3.0; // initial guess (centroid)
  double snew   = 1.0 / 3.0;
  double xinew  = 0.0;
  double rcur   = rnew;
  double scur   = snew;
  double xicur  = xinew;
  double xidiff[] = { 1.0, 1.0, 1.0 };

  double shp_func_deriv[18];
  double current_pc[3];

  const int MAX_NR_ITER = 20;
  int i = 0;
  do
  {
    current_pc[0] = rcur  = rnew;
    current_pc[1] = scur  = snew;
    current_pc[2] = xicur = xinew;

    // Build Jacobian and Invert

    //aj(1,1)=( dN/dr  ) * x[]
    //aj(1,2)=( dN/ds  ) * x[]
    //aj(1,3)=( dN/dxi ) * x[]
    //aj(2,1)=( dN/dr  ) * y[]
    //aj(2,2)=( dN/ds  ) * y[]
    //aj(2,3)=( dN/dxi ) * y[]
    //aj(3,1)=( dN/dr  ) * z[]
    //aj(3,2)=( dN/ds  ) * z[]
    //aj(3,3)=( dN/dxi ) * z[]

    wedge_derivative(1, current_pc, shp_func_deriv);

    for (int row = 0; row != 3; ++row)
      for (int col = 0; col != 3; ++col)
	j[row][col] = 0.0;

    for (int k = 1; k != 6; ++k)
    {
      j[0][0] -= shp_func_deriv[k * 3 + 0] * x[k];
      j[0][1] -= shp_func_deriv[k * 3 + 1] * x[k];
      j[0][2] -= shp_func_deriv[k * 3 + 2] * x[k];

      j[1][0] -= shp_func_deriv[k * 3 + 0] * y[k];
      j[1][1] -= shp_func_deriv[k * 3 + 1] * y[k];
      j[1][2] -= shp_func_deriv[k * 3 + 2] * y[k];

      j[2][0] -= shp_func_deriv[k * 3 + 0] * z[k];
      j[2][1] -= shp_func_deriv[k * 3 + 1] * z[k];
      j[2][2] -= shp_func_deriv[k * 3 + 2] * z[k];
    }

    const double jdet =   j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1])
		      - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
		      + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);

    jinv[0][0] =  (j[1][1] * j[2][2] - j[1][2] * j[2][1]) / jdet;
    jinv[0][1] = -(j[0][1] * j[2][2] - j[2][1] * j[0][2]) / jdet;
    jinv[0][2] =  (j[1][2] * j[0][1] - j[0][2] * j[1][1]) / jdet;
    jinv[1][0] = -(j[1][0] * j[2][2] - j[2][0] * j[1][2]) / jdet;
    jinv[1][1] =  (j[0][0] * j[2][2] - j[0][2] * j[2][0]) / jdet;
    jinv[1][2] = -(j[0][0] * j[1][2] - j[1][0] * j[0][2]) / jdet;
    jinv[2][0] =  (j[1][0] * j[2][1] - j[2][0] * j[1][1]) / jdet;
    jinv[2][1] = -(j[0][0] * j[2][1] - j[2][0] * j[0][1]) / jdet;
    jinv[2][2] =  (j[0][0] * j[1][1] - j[0][1] * j[1][0]) / jdet;

    wedge_shape_fcn(1, current_pc, shapefct);

    // x[0] = y[0] = z[0] = 0 by construction
    f[0] = xp - (shapefct[1] * x[1] +
		 shapefct[2] * x[2] +
		 shapefct[3] * x[3] +
		 shapefct[4] * x[4] +
		 shapefct[5] * x[5]);
    f[1] = yp - (shapefct[1] * y[1] +
		 shapefct[2] * y[2] +
		 shapefct[3] * y[3] +
		 shapefct[4] * y[4] +
		 shapefct[5] * y[5]);
    f[2] = zp - (shapefct[1] * z[1] +
		 shapefct[2] * z[2] +
		 shapefct[3] * z[3] +
		 shapefct[4] * z[4] +
		 shapefct[5] * z[5]);

    rnew  = rcur  - (f[0] * jinv[0][0] + f[1] * jinv[0][1] + f[2] * jinv[0][2]);
    snew  = scur  - (f[0] * jinv[1][0] + f[1] * jinv[1][1] + f[2] * jinv[1][2]);
    xinew = xicur - (f[0] * jinv[2][0] + f[1] * jinv[2][1] + f[2] * jinv[2][2]);

    xidiff[0] = rnew  - rcur;
    xidiff[1] = snew  - scur;
    xidiff[2] = xinew - xicur;
  }
  while (!within_tolerance(vector_norm_sq(xidiff,3), isInElemConverged) && ++i != MAX_NR_ITER);

  isoParCoord[0] = isoParCoord[1] = isoParCoord[2] = std::numeric_limits<double>::max();
  double dist = std::numeric_limits<double>::max();

  if (i < MAX_NR_ITER)
  {
    isoParCoord[0] = rnew;
    isoParCoord[1] = snew;
    isoParCoord[2] = xinew;
    std::vector<double> xx = { isoParCoord[0], isoParCoord[1], isoParCoord[2] };

    dist = parametric_distance(xx);
  }
  return dist;
}

//--------------------------------------------------------------------------
//-------- interpolatePoint ------------------------------------------------
//--------------------------------------------------------------------------
void
WedSCS::interpolatePoint(
  const int &nComp,
  const double *isoParCoord,
  const double *field,
  double *result )
{
  double shapefct[6];

  wedge_shape_fcn(1, isoParCoord, shapefct);

  for ( int i = 0; i < nComp; i++)
  {
    // Base 'field array' index for i_th component
    int b = 6 * i;

    result[i] = 0.0;

    for (int j = 0; j != 6; ++j)
      result[i] += shapefct[j] * field[b + j];
  }
}

//--------------------------------------------------------------------------
//-------- wedge_shape_fcn -------------------------------------------------
//--------------------------------------------------------------------------
void
WedSCS::wedge_shape_fcn(
  const int  &npts,
  const double *isoParCoord, 
  double *shape_fcn)
{
  for (int j = 0; j < npts; ++j ) {
    int sixj = 6 * j;
    int k    = 3 * j;
    double r    = isoParCoord[k];
    double s    = isoParCoord[k + 1];
    double t    = 1.0 - r - s;
    double xi   = isoParCoord[k + 2];
    shape_fcn[    sixj] = 0.5 * t * (1.0 - xi);
    shape_fcn[1 + sixj] = 0.5 * r * (1.0 - xi);
    shape_fcn[2 + sixj] = 0.5 * s * (1.0 - xi);
    shape_fcn[3 + sixj] = 0.5 * t * (1.0 + xi);
    shape_fcn[4 + sixj] = 0.5 * r * (1.0 + xi);
    shape_fcn[5 + sixj] = 0.5 * s * (1.0 + xi);
  }
}

//--------------------------------------------------------------------------
//-------- parametric_distance ---------------------------------------------
//--------------------------------------------------------------------------
double
WedSCS::parametric_distance(const double X, const double Y)
{
  const double dist0 = -3*X;
  const double dist1 = -3*Y;
  const double dist2 =  3*(X+Y);
  const double dist = std::max(std::max(dist0,dist1),dist2);
  return dist;
}

//--------------------------------------------------------------------------
//-------- parametric_distance ---------------------------------------------
//--------------------------------------------------------------------------
double
WedSCS::parametric_distance(const std::vector<double> &x)
{
  const double X = x[0] - 1./3.;
  const double Y = x[1] - 1./3.;
  const double Z = x[2] ;
  const double dist_t = parametric_distance(X,Y);
  const double dist_z = std::fabs(Z);
  const double dist = std::max(dist_z, dist_t);
  return dist;
}

//--------------------------------------------------------------------------
//-------- general_face_grad_op --------------------------------------------
//--------------------------------------------------------------------------
void 
WedSCS::general_face_grad_op(
  const int face_ordinal,
  const double *isoParCoord,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  int lerr = 0;
  const int nface = 1;
  double dpsi[18], grad[18];
      
  wedge_derivative(nface, &isoParCoord[0], dpsi);
      
  SIERRA_FORTRAN(wed_gradient_operator) 
    ( &nface,
      &nodesPerElement_,
      &nface,
      dpsi,
      &coords[0], grad, &det_j[0], error, &lerr );
      
  if ( lerr )
    std::cout << "problem with EwedSCS::general_face_grad" << std::endl;
  
  for ( int j=0; j<18; ++j) {
    gradop[j] = grad[j];
  }  
}

//--------------------------------------------------------------------------
//-------- sidePcoords_to_elemPcoords --------------------------------------
//--------------------------------------------------------------------------
void 
WedSCS::sidePcoords_to_elemPcoords(
  const int & side_ordinal,
  const int & npoints,
  const double *side_pcoords,
  double *elem_pcoords)
{
  switch (side_ordinal) {
  case 0:
    for (int i=0; i<npoints; i++) {//face0:quad: (x,y) -> (0.5*(1 + x),0,y)
      elem_pcoords[i*3+0] = 0.5*(1.0+side_pcoords[2*i+0]);
      elem_pcoords[i*3+1] = 0.0;
      elem_pcoords[i*3+2] = side_pcoords[2*i+1];
    }
    break;
  case 1:
    for (int i=0; i<npoints; i++) {//face1:quad: (x,y) -> (0.5*(1-y),0.5*(1 + y),x)
      elem_pcoords[i*3+0] = 0.5*(1.0-side_pcoords[2*i+1]);
      elem_pcoords[i*3+1] = 0.5*(1.0+side_pcoords[2*i+1]);
      elem_pcoords[i*3+2] = side_pcoords[2*i+0];
    }
    break;
  case 2:
    for (int i=0; i<npoints; i++) {//face2:quad: (x,y) -> (0,0.5*(1 + x),y)
      elem_pcoords[i*3+0] = 0.0;
      elem_pcoords[i*3+1] = 0.5*(1.0+side_pcoords[2*i+0]);
      elem_pcoords[i*3+2] = side_pcoords[2*i+1];
    }
    break;
  case 3:
    for (int i=0; i<npoints; i++) {//face3:tri: (x,y) -> (x,y,-1)
      elem_pcoords[i*3+0] = side_pcoords[2*i+0];
      elem_pcoords[i*3+1] = side_pcoords[2*i+1];
      elem_pcoords[i*3+2] = -1.0;
    }
    break;
  case 4:
    for (int i=0; i<npoints; i++) {//face4:tri: (x,y) -> (x,y,+1 )
      elem_pcoords[i*3+0] = side_pcoords[2*i+0];
      elem_pcoords[i*3+1] = side_pcoords[2*i+1];
      elem_pcoords[i*3+2] = 1.0;
    }
    break;
  default:
    throw std::runtime_error("WedSCS::sidePcoords_to_elemPcoords invalid ordinal");
  }
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Quad2DSCV::Quad2DSCV()
  : MasterElement()
{
  nDim_ = 2;
  nodesPerElement_ = 4;
  numIntPoints_ = 4;

  // define ip node mappings
  ipNodeMap_.resize(4);
  ipNodeMap_[0] = 0; ipNodeMap_[1] = 1; ipNodeMap_[2] = 2; ipNodeMap_[3] = 3;

  // standard integration location
  intgLoc_.resize(8);    
  intgLoc_[0]  = -0.25; intgLoc_[1]  = -0.25; 
  intgLoc_[2]  = +0.25; intgLoc_[3]  = -0.25; 
  intgLoc_[4]  = +0.25; intgLoc_[5]  = +0.25; 
  intgLoc_[6]  = -0.25; intgLoc_[7]  = +0.25; 

  // shifted integration location
  intgLocShift_.resize(8);    
  intgLocShift_[0]  = -0.50; intgLocShift_[1]  = -0.50; 
  intgLocShift_[2]  = +0.50; intgLocShift_[3]  = -0.50; 
  intgLocShift_[4]  = +0.50; intgLocShift_[5]  = +0.50; 
  intgLocShift_[6]  = -0.50; intgLocShift_[7]  = +0.50; 
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Quad2DSCV::~Quad2DSCV()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad2DSCV::ipNodeMap(
  int /*ordinal*/)
{
  // define scv->node mappings
  return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void Quad2DSCV::determinant(
  const int nelem,
  const double *coords,
  double *volume,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(quad_scv_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords,
      volume, error, &lerr );
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Quad2DSCV::shape_fcn(double *shpfc)
{
  quad_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Quad2DSCV::shifted_shape_fcn(double *shpfc)
{
  quad_shape_fcn(numIntPoints_, &intgLocShift_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- quad_shape_fcn ---------------------------------------------------
//--------------------------------------------------------------------------
void
Quad2DSCV::quad_shape_fcn(
  const int  &npts,
  const double *isoParCoord, 
  double *shape_fcn)
{
  for (int j = 0; j < npts; ++j ) {
    const int fourj = 4*j;
    const int k = 2*j;
    const double s1 = isoParCoord[k];
    const double s2 = isoParCoord[k+1];
    shape_fcn[    fourj] = 1.0/4.0 + 0.5*(-s1 - s2 ) + s1*s2;
    shape_fcn[1 + fourj] = 1.0/4.0 + 0.5*( s1 - s2 ) - s1*s2;
    shape_fcn[2 + fourj] = 1.0/4.0 + 0.5*( s1 + s2 ) + s1*s2;
    shape_fcn[3 + fourj] = 1.0/4.0 + 0.5*(-s1 + s2 ) - s1*s2;
  }
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Quad2DSCS::Quad2DSCS()
  : MasterElement()
{
  nDim_ = 2;
  nodesPerElement_ = 4;
  numIntPoints_ = 4;
  scaleToStandardIsoFac_ = 2.0;

  // define L/R mappings
  lrscv_.resize(8);
  lrscv_[0]  = 0; lrscv_[1]  = 1;
  lrscv_[2]  = 1; lrscv_[3]  = 2;
  lrscv_[4]  = 2; lrscv_[5]  = 3;
  lrscv_[6]  = 0; lrscv_[7]  = 3;
  
  // define opposing node
  oppNode_.resize(8);
  // face 0; nodes 0,1
  oppNode_[0] = 3; oppNode_[1] = 2;
  // face 1; nodes 1,2
  oppNode_[2] = 0; oppNode_[3] = 3;
  // face 2; nodes 2,3
  oppNode_[4] = 1; oppNode_[5] = 0;
  // face 3; nodes 3,0
  oppNode_[6] = 2; oppNode_[7] = 1;

  // define opposing face
  oppFace_.resize(8);
  // face 0
  oppFace_[0]  = 3; oppFace_[1] = 1;
  // face 1
  oppFace_[2]  = 0; oppFace_[3] = 2;
  // face 2
  oppFace_[4]  = 1; oppFace_[5] = 3;  
  // face 3
  oppFace_[6]  = 2; oppFace_[7] = 0;

  // standard integration location
  intgLoc_.resize(8);    
  intgLoc_[0] =  0.00; intgLoc_[1] = -0.25; // surf 1; 1->2
  intgLoc_[2] =  0.25; intgLoc_[3] =  0.00; // surf 2; 2->3
  intgLoc_[4] =  0.00; intgLoc_[5] =  0.25; // surf 3; 3->4
  intgLoc_[6] = -0.25; intgLoc_[7] =  0.00; // surf 3; 1->5

  // shifted
  intgLocShift_.resize(8);
  intgLocShift_[0] =  0.00; intgLocShift_[1] = -0.50;
  intgLocShift_[2] =  0.50; intgLocShift_[3] =  0.00;
  intgLocShift_[4] =  0.00; intgLocShift_[5] =  0.50;
  intgLocShift_[6] = -0.50; intgLocShift_[7] =  0.00;

  // exposed face
  intgExpFace_.resize(16);
  // face 0; scs 0, 1; nodes 0,1
  intgExpFace_[0]  = -0.25; intgExpFace_[1]  = -0.50; 
  intgExpFace_[2]  =  0.25; intgExpFace_[3]  = -0.50;
  // face 1; scs 0, 1; nodes 1,2
  intgExpFace_[4]  =  0.50; intgExpFace_[5]  = -0.25;
  intgExpFace_[6]  =  0.50; intgExpFace_[7]  =  0.25;
  // face 2, surf 0, 1; nodes 2,3
  intgExpFace_[8]  =  0.25; intgExpFace_[9]  =  0.50;
  intgExpFace_[10] = -0.25; intgExpFace_[11] =  0.50;
  // face 3, surf 0, 1; nodes 3,0
  intgExpFace_[12] = -0.50; intgExpFace_[13] =  0.25;
  intgExpFace_[14] = -0.50; intgExpFace_[15] = -0.25; 

  // boundary integration point ip node mapping (ip on an ordinal to local node number)
  ipNodeMap_.resize(8); // 2 ips * 4 faces
  // face 0;
  ipNodeMap_[0] = 0;  ipNodeMap_[1] = 1;  
  // face 1;
  ipNodeMap_[2] = 1;  ipNodeMap_[3] = 2; 
  // face 2; 
  ipNodeMap_[4] = 2;  ipNodeMap_[5] = 3;  
  // face 3;
  ipNodeMap_[6] = 3;  ipNodeMap_[7] = 0; 
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Quad2DSCS::~Quad2DSCS()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad2DSCS::ipNodeMap(
  int ordinal)
{
  // define ip->node mappings for each face (ordinal); 
  return &ipNodeMap_[ordinal*2];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void Quad2DSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  SIERRA_FORTRAN(quad_scs_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords, areav );

  // all is always well; no error checking
  *error = 0;
}

//--------------------------------------------------------------------------
//-------- grad_op ---------------------------------------------------------
//--------------------------------------------------------------------------
void Quad2DSCS::grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(quad_derivative)
    ( &numIntPoints_, &intgLoc_[0], deriv );
  
  SIERRA_FORTRAN(quad_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );
  
  if ( lerr )
    std::cout << "sorry, negative Quad2DSCS volume.." << std::endl;  
}

//--------------------------------------------------------------------------
//-------- shifted_grad_op -------------------------------------------------
//--------------------------------------------------------------------------
void Quad2DSCS::shifted_grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(quad_derivative)
    ( &numIntPoints_, &intgLocShift_[0], deriv );

  SIERRA_FORTRAN(quad_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );

  if ( lerr )
    std::cout << "sorry, negative Quad2DSCS volume.." << std::endl;
}

//--------------------------------------------------------------------------
//-------- face_grad_op ----------------------------------------------------
//--------------------------------------------------------------------------
void Quad2DSCS::face_grad_op(
  const int nelem,
  const int face_ordinal,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  int lerr = 0;
  int npf = 2;
  int ndim = 2;

  const int nface = 1;
  double dpsi[8];
  double grad[8];

  for ( int n=0; n<nelem; n++ ) {
    
    for ( int k=0; k<npf; k++ ) {
      
      const int row = 4*face_ordinal + k*ndim;
      SIERRA_FORTRAN(quad_derivative)
        ( &nface, &intgExpFace_[row], dpsi );
      
      SIERRA_FORTRAN(quad_gradient_operator)
        ( &nface,
          &nodesPerElement_,
          &nface,
          dpsi,
          &coords[8*n], grad, &det_j[npf*n+k], error, &lerr );
      
      if ( lerr )
        std::cout << "sorry, issue with face_grad_op.." << std::endl;
      
      for ( int j=0; j<8; j++) {
        gradop[k*nelem*8+n*8+j] = grad[j];
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- gij -------------------------------------------------------------
//--------------------------------------------------------------------------
void Quad2DSCS::gij(
  const double *coords,
  double *gupperij,
  double *glowerij,
  double *deriv)
{
  SIERRA_FORTRAN(twod_gij)
    ( &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gupperij, glowerij);
}

//--------------------------------------------------------------------------
//-------- adjacentNodes ---------------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad2DSCS::adjacentNodes()
{
  // define L/R mappings
  return &lrscv_[0];
}

//--------------------------------------------------------------------------
//-------- opposingNodes --------------------------------------------------
//--------------------------------------------------------------------------
int
Quad2DSCS::opposingNodes(
  const int ordinal,
  const int node)
{
  return oppNode_[ordinal*2+node];
}

//--------------------------------------------------------------------------
//-------- opposingFace --------------------------------------------------
//--------------------------------------------------------------------------
int
Quad2DSCS::opposingFace(
  const int ordinal,
  const int node)
{
  return oppFace_[ordinal*2+node];
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Quad2DSCS::shape_fcn(double *shpfc)
{
  quad_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Quad2DSCS::shifted_shape_fcn(double *shpfc)
{
  quad_shape_fcn(numIntPoints_, &intgLocShift_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- quad_shape_fcn ---------------------------------------------------
//--------------------------------------------------------------------------
void
Quad2DSCS::quad_shape_fcn(
  const int  &npts,
  const double *isoParCoord, 
  double *shape_fcn)
{
  for (int j = 0; j < npts; ++j ) {
    const int fourj = 4*j;
    const int k = 2*j;
    const double s1 = isoParCoord[k];
    const double s2 = isoParCoord[k+1];
    shape_fcn[    fourj] = 1.0/4.0 + 0.5*(-s1 - s2 ) + s1*s2;
    shape_fcn[1 + fourj] = 1.0/4.0 + 0.5*( s1 - s2 ) - s1*s2;
    shape_fcn[2 + fourj] = 1.0/4.0 + 0.5*( s1 + s2 ) + s1*s2;
    shape_fcn[3 + fourj] = 1.0/4.0 + 0.5*(-s1 + s2 ) - s1*s2;
  }
}

//--------------------------------------------------------------------------
//-------- isInElement -----------------------------------------------------
//--------------------------------------------------------------------------
double
Quad2DSCS::isInElement(
  const double *elemNodalCoord,
  const double *pointCoord,
  double *isoParCoord )
{
  // square of the desired norm, 1.0e-8
  const double isInElemConverged = 1.0e-16;
  const int maxNonlinearIter = 10;

  // -1:1 isoparametric range
  
  // Translate element so that (x,y) coordinates of the first node are (0,0)
  double x[4] = {0.,
                 elemNodalCoord[1] - elemNodalCoord[0],
                 elemNodalCoord[2] - elemNodalCoord[0],
                 elemNodalCoord[3] - elemNodalCoord[0] };
  double y[4] = {0.,
                 elemNodalCoord[5] - elemNodalCoord[4],
                 elemNodalCoord[6] - elemNodalCoord[4],
                 elemNodalCoord[7] - elemNodalCoord[4] };
  
  // (xp,yp) is the point at which we're searching for (xi,eta)
  // (must translate this also)
  
  double xp = pointCoord[0] - elemNodalCoord[0];
  double yp = pointCoord[1] - elemNodalCoord[4];
  
  // Newton-Raphson iteration for (xi,eta)
  double j[4];
  double f[2];
  double shapefct[4];
  
  double xinew = 0.5;     // initial guess
  double etanew = 0.5;
  
  double xicur = 0.5;
  double etacur = 0.5;
  
  double xidiff[2] = { 1.0, 1.0};
  int i = 0;
  
  bool converged = false;
  
  do {
    xicur = xinew;
    etacur = etanew;

    j[0]=  0.25*(1.00-etacur)*x[1]
	  +0.25*(1.00+etacur)*x[2]
	  -0.25*(1.00+etacur)*x[3];

    j[1]= -0.25*(1.00+xicur)*x[1]
	  +0.25*(1.00+xicur)*x[2]
	  +0.25*(1.00-xicur)*x[3];

    j[2]=  0.25*(1.00-etacur)*y[1]
	  +0.25*(1.00+etacur)*y[2]
	  -0.25*(1.00+etacur)*y[3];

    j[3]= -0.25*(1.00+xicur)*y[1]
	  +0.25*(1.00+xicur)*y[2]
	  +0.25*(1.00-xicur)*y[3];

    double jdet = j[0]*j[3] - j[1]*j[2];

    shapefct[0]=0.25*(1.00-etacur)*(1.00-xicur);
    shapefct[1]=0.25*(1.00-etacur)*(1.00+xicur);
    shapefct[2]=0.25*(1.00+etacur)*(1.00+xicur);
    shapefct[3]=0.25*(1.00+etacur)*(1.00-xicur);


    f[0] = (shapefct[1]*x[1]+shapefct[2]*x[2]+shapefct[3]*x[3]) - xp;
    f[1] = (shapefct[1]*y[1]+shapefct[2]*y[2]+shapefct[3]*y[3]) - yp;

    xinew  = xicur  - ( f[0]*j[3] - f[1]*j[1])/jdet;
    etanew = etacur - (-f[0]*j[2] + f[1]*j[0])/jdet;

    xidiff[0] = xinew  - xicur;
    xidiff[1] = etanew - etacur;
    
    double vectorNorm = xidiff[0]*xidiff[0] + xidiff[1]*xidiff[1];
    converged = (vectorNorm < isInElemConverged);
  }  while ( !converged && (++i < maxNonlinearIter) );

  // set a bad value
  isoParCoord[0] = isoParCoord[1] = 1.0e6;
  double dist = 1.0e6;
  if ( i < maxNonlinearIter ) {
    isoParCoord[0] = xinew;
    isoParCoord[1] = etanew;
    dist = (std::abs(xinew) > std::abs(etanew))
      ? std::abs(xinew) : std::abs(etanew);
  }
  return dist;
}

//--------------------------------------------------------------------------
//-------- interpolatePoint ------------------------------------------------
//--------------------------------------------------------------------------
void
Quad2DSCS::interpolatePoint(
  const int &nComp,
  const double *isoParCoord,
  const double *field,
  double *result )
{
  // -1:1 isoparametric range
  const double xi   = isoParCoord[0];
  const double eta  = isoParCoord[1];

  for ( int i = 0; i < nComp; i++ )
  {
    // Base 'field array' index for ith component
    int b = 4*i;

    result[i] = 0.250000000000000 * (
      (1.00000000000000-eta) * (1.00000000000000-xi ) * field[b+0] +
      (1.00000000000000-eta) * (1.00000000000000+xi ) * field[b+1] +
      (1.00000000000000+eta) * (1.00000000000000+xi ) * field[b+2] +
      (1.00000000000000+eta) * (1.00000000000000-xi ) * field[b+3] ) ;
  }  
}

//--------------------------------------------------------------------------
//-------- general_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Quad2DSCS::general_shape_fcn(
  const int numIp,
  const double *isoParCoord,
  double *shpfc)
{
  // -1:1 isoparametric range
  const double npe = nodesPerElement_;
  for ( int ip = 0; ip < numIp; ++ip ) {
    
    const int rowIpc = 2*ip;
    const int rowSfc = npe*ip;
    
    const double s1 = isoParCoord[rowIpc];
    const double s2 = isoParCoord[rowIpc+1];
    shpfc[rowSfc  ] = 0.25*(1.0-s1)*(1.0-s2);
    shpfc[rowSfc+1] = 0.25*(1.0+s1)*(1.0-s2);
    shpfc[rowSfc+2] = 0.25*(1.0+s1)*(1.0+s2);
    shpfc[rowSfc+3] = 0.25*(1.0-s1)*(1.0+s2);
    
  }
}

//--------------------------------------------------------------------------
//-------- general_face_grad_op --------------------------------------------
//--------------------------------------------------------------------------
void 
Quad2DSCS::general_face_grad_op(
  const int face_ordinal,
  const double *isoParCoord,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  int lerr = 0;
  const int nface = 1;

  double dpsi[8];
  double grad[8];

  SIERRA_FORTRAN(quad_derivative)
    ( &nface, &isoParCoord[0], dpsi );
      
  SIERRA_FORTRAN(quad_gradient_operator)
    ( &nface,
      &nodesPerElement_,
      &nface,
      dpsi,
      &coords[0], grad, &det_j[0], error, &lerr );
  
  if ( lerr )
    std::cout << "Quad2DSCS::general_face_grad_op: issue.." << std::endl;
  
  for ( int j=0; j<8; j++) {
    gradop[j] = grad[j];
  }
}

//--------------------------------------------------------------------------
//-------- sidePcoords_to_elemPcoords --------------------------------------
//--------------------------------------------------------------------------
void 
Quad2DSCS::sidePcoords_to_elemPcoords(
  const int & side_ordinal,
  const int & npoints,
  const double *side_pcoords,
  double *elem_pcoords)
{
  switch (side_ordinal) {
  case 0:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*2+0] = 0.5*side_pcoords[i];
      elem_pcoords[i*2+1] = -0.5;
    }
    break;
  case 1:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*2+0] = 0.5;
      elem_pcoords[i*2+1] = 0.5*side_pcoords[i];
    }
    break;
  case 2:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*2+0] = -0.5*side_pcoords[i];
      elem_pcoords[i*2+1] = 0.5;
    }
    break;
  case 3:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*2+0] = -0.5;
      elem_pcoords[i*2+1] = -0.5*side_pcoords[i];
    }
    break;
  default:
    throw std::runtime_error("Quad2DSCS::sideMap invalid ordinal");
  }
}

//--------------------------------------------------------------------------
//-------- constructor------------------------------------------------------
//--------------------------------------------------------------------------
QuadrilateralP2Element::QuadrilateralP2Element()
  : MasterElement(),
    scsDist_(std::sqrt(3.0)/3.0),
    useGLLGLL_(false),
    nodes1D_(3),
    numQuad_(2)
{
  //TODO: Avoid redundant gradient computations when useGLLGLL_ is enabled

  nDim_ = 2;
  nodesPerElement_ = nodes1D_ * nodes1D_;

  // map the standard stk (refinement consistent) node numbering
  // to a tensor-product style node numbering (i.e. node (m,l,k) -> m+npe*l+npe^2*k)
  stkNodeMap_ = {
                  0, 4, 1, // bottom row of nodes
                  7, 8, 5, // middle row of nodes
                  3, 6, 2  // top row of nodes
                };

  // a padded list of scs locations
  scsEndLoc_ = { -1.0, -scsDist_, scsDist_, +1.0 };
}

//--------------------------------------------------------------------------
//-------- set_quadrature_rule ---------------------------------------------
//--------------------------------------------------------------------------
void
QuadrilateralP2Element::set_quadrature_rule()
{
  gaussAbscissaeShift_ = {-1.0,-1.0,0.0,0.0,+1.0,+1.0};

  if (useGLLGLL_ && numQuad_ < 3) {
    throw std::runtime_error("Need at least 3 points for gllgll quadrature");
  }

  switch (numQuad_) {
    case 1:
    {
      gaussAbscissae_ = { 0.0 };
      gaussWeight_ = { 1.0 };
      break;
    }
    case 2:
    {
      gaussAbscissae_ = { -std::sqrt(3.0)/3.0, std::sqrt(3.0)/3.0 };
      gaussWeight_ = { 0.5, 0.5 };
      break;
    }
    case 3:
    {
      if (!useGLLGLL_) {
        gaussAbscissae_ = { -std::sqrt(3.0/5.0), 0.0, std::sqrt(3.0/5.0) };
        gaussWeight_ = { 5.0/18.0, 4.0/9.0,  5.0/18.0 };
      }
      else {
        gaussAbscissae1D_ = {-1.0, 0.0, +1.0};
        GLLGLL_quadrature_weights();
        gaussAbscissae_ = {
            gaussAbscissae1D_[0], gaussAbscissae1D_[1], gaussAbscissae1D_[2],
            gaussAbscissae1D_[0], gaussAbscissae1D_[1], gaussAbscissae1D_[2],
            gaussAbscissae1D_[0], gaussAbscissae1D_[1], gaussAbscissae1D_[2]
        };
      }
      break;
    }
    case 4:
    {
      if (!useGLLGLL_) {
        gaussAbscissae_ = {
            -std::sqrt(3.0/7.0+2.0/7.0*std::sqrt(6.0/5.0)),
            -std::sqrt(3.0/7.0-2.0/7.0*std::sqrt(6.0/5.0)),
            +std::sqrt(3.0/7.0-2.0/7.0*std::sqrt(6.0/5.0)),
            +std::sqrt(3.0/7.0+2.0/7.0*std::sqrt(6.0/5.0))
        };

        gaussWeight_ = {
            (18.0-std::sqrt(30.0))/72.0,
            (18.0+std::sqrt(30.0))/72.0,
            (18.0+std::sqrt(30.0))/72.0,
            (18.0-std::sqrt(30.0))/72.0
        };
      }
      else {
        gaussAbscissae1D_ = { -1.0, -std::sqrt(5.0)/5.0, std::sqrt(5.0)/5.0, +1.0};
        GLLGLL_quadrature_weights();

        gaussAbscissae_ = {
            gaussAbscissae1D_[0], gaussAbscissae1D_[1], gaussAbscissae1D_[2], gaussAbscissae1D_[3],
            gaussAbscissae1D_[0], gaussAbscissae1D_[1], gaussAbscissae1D_[2], gaussAbscissae1D_[3],
            gaussAbscissae1D_[0], gaussAbscissae1D_[1], gaussAbscissae1D_[2], gaussAbscissae1D_[3]
        };
      }
      break;
    }
    default:
      throw std::runtime_error("Quadrature rule not implemented");
  }
}

//--------------------------------------------------------------------------
//-------- GLLGLL_quadrature_weights ---------------------------------------
//--------------------------------------------------------------------------
void
QuadrilateralP2Element::GLLGLL_quadrature_weights()
{
  // for fixed abscissae, use a moment-matching equation to compute weights
  // given a heaviside weight function
  const int nrows = gaussAbscissae1D_.size();
  const int nrhs = nodes1D_;

  Teuchos::SerialDenseMatrix<int, double> weightLHS(nrows, nrows);
  for (int j = 0; j < nrows; ++j) {
    for (int i = 0; i < nrows; ++i) {
      weightLHS(i, j) = std::pow(gaussAbscissae1D_[j], i);
    }
  }

  // each node has a separate RHS
  Teuchos::SerialDenseMatrix<int, double> weightRHS(nrows, nrhs);
  for (int j = 0; j < nrhs; ++j) {
    for (int i = 0; i < nrows; ++i) {
      weightRHS(i, j) = (std::pow(scsEndLoc_[j + 1], i + 1)
                       - std::pow(scsEndLoc_[j], i + 1)) / (i + 1.0);
    }
  }

  Teuchos::SerialDenseSolver<int, double> solver;
  Teuchos::SerialDenseMatrix<int, double> quadratureWeights(nrows, nrhs);
  solver.setMatrix(Teuchos::rcp(&weightLHS, false));
  solver.setVectors(
    Teuchos::rcp(&quadratureWeights, false),
    Teuchos::rcp(&weightRHS, false)
  );
  solver.solve();

  gaussWeight_.resize(nrows * nrhs);
  for (int j = 0; j < nrows; ++j) {
    for (int i = 0; i < nrhs; ++i) {
      gaussWeight_[i + j * nrhs] = quadratureWeights(i, j); //transpose
    }
  }
}

//--------------------------------------------------------------------------
//-------- tensor_product_node_map -----------------------------------------
//--------------------------------------------------------------------------
int
QuadrilateralP2Element::tensor_product_node_map(int i, int j) const
{
   return stkNodeMap_[i+nodes1D_*j];
}

//--------------------------------------------------------------------------
//-------- gauss_point_location --------------------------------------------
//--------------------------------------------------------------------------
double
QuadrilateralP2Element::gauss_point_location(
  int nodeOrdinal,
  int gaussPointOrdinal) const
{
  double location1D;
  if (!useGLLGLL_) {
    location1D = isoparametric_mapping( scsEndLoc_[nodeOrdinal+1],
                                        scsEndLoc_[nodeOrdinal],
                                        gaussAbscissae_[gaussPointOrdinal] );
  }
  else {
    location1D = gaussAbscissae_[nodeOrdinal*numQuad_ + gaussPointOrdinal];
  }
   return location1D;
}
//--------------------------------------------------------------------------
//-------- shifted_gauss_point_location ------------------------------------
//--------------------------------------------------------------------------
double
QuadrilateralP2Element::shifted_gauss_point_location(
  int nodeOrdinal,
  int gaussPointOrdinal) const
{
  return gaussAbscissaeShift_[nodeOrdinal*numQuad_ + gaussPointOrdinal];
}
//--------------------------------------------------------------------------
//-------- tensor_product_weight -------------------------------------------
//--------------------------------------------------------------------------
double
QuadrilateralP2Element::tensor_product_weight(
  int s1Node, int s2Node,
  int s1Ip, int s2Ip) const
{
  //surface integration
  double weight;
  if (!useGLLGLL_) {
    const double Ls1 = scsEndLoc_[s1Node+1]-scsEndLoc_[s1Node];
    const double Ls2 = scsEndLoc_[s2Node+1]-scsEndLoc_[s2Node];
    const double isoparametricArea = Ls1 * Ls2;

    weight = isoparametricArea * gaussWeight_[s1Ip] * gaussWeight_[s2Ip];
   }
   else {
     // weights are node-specific and take into account the isoparametric area
     weight = gaussWeight_[s1Node+nodes1D_*s1Ip] * gaussWeight_[s2Node+nodes1D_*s2Ip];
   }
   return weight;
}

//--------------------------------------------------------------------------
//-------- tensor_product_weight -------------------------------------------
//--------------------------------------------------------------------------
double
QuadrilateralP2Element::tensor_product_weight(int s1Node, int s1Ip) const
{
  //line integration
  double weight;
  if (!useGLLGLL_) {
    const double isoparametricLength = scsEndLoc_[s1Node+1]-scsEndLoc_[s1Node];
    weight = isoparametricLength * gaussWeight_[s1Ip];
   }
   else {
     // weights are node-specific and take into account the isoparametric area
     weight = gaussWeight_[s1Node+nodes1D_*s1Ip];
   }
   return weight;
}


//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
QuadrilateralP2Element::shape_fcn(double* shpfc)
{
  for (int ni = 0; ni < numIntPoints_ * nodesPerElement_; ++ni) {
    shpfc[ni] = shapeFunctions_[ni];
  }
}

void
QuadrilateralP2Element::shifted_shape_fcn(double* shpfc)
{
  for (int ip = 0; ip < numIntPoints_ * nodesPerElement_; ++ip) {
    shpfc[ip] = shapeFunctionsShift_[ip];
  }
}
//--------------------------------------------------------------------------
//-------- eval_shape_functions_at_ips -------------------------------------
//--------------------------------------------------------------------------
void
QuadrilateralP2Element::eval_shape_functions_at_ips()
{
  shapeFunctions_.resize(numIntPoints_*nodesPerElement_);
  quad9_shape_fcn(numIntPoints_, intgLoc_.data(), shapeFunctions_.data());
}

//--------------------------------------------------------------------------
//-------- eval_shape_derivs_at_ips ----------------------------------------
//--------------------------------------------------------------------------
void
QuadrilateralP2Element::eval_shape_derivs_at_ips()
{
  shapeDerivs_.resize(numIntPoints_*nodesPerElement_*nDim_);
  quad9_shape_deriv(numIntPoints_, intgLoc_.data(), shapeDerivs_.data());
}

//--------------------------------------------------------------------------
//-------- eval_shape_functions_at_shifted_ips -----------------------------
//--------------------------------------------------------------------------
void
QuadrilateralP2Element::eval_shape_functions_at_shifted_ips()
{
  shapeFunctionsShift_.resize(numIntPoints_*nodesPerElement_);
  quad9_shape_fcn(numIntPoints_, intgLocShift_.data(), shapeFunctionsShift_.data());
}

//--------------------------------------------------------------------------
//-------- eval_shape_derivs_at_ips ----------------------------------------
//--------------------------------------------------------------------------
void
QuadrilateralP2Element::eval_shape_derivs_at_shifted_ips()
{
  shapeDerivsShift_.resize(numIntPoints_*nodesPerElement_*nDim_);
  quad9_shape_deriv(numIntPoints_, intgLocShift_.data(), shapeDerivsShift_.data());
}
//--------------------------------------------------------------------------
//-------- eval_shape_derivs_at_face_ips ----------------------------------------
//--------------------------------------------------------------------------
void
QuadrilateralP2Element::eval_shape_derivs_at_face_ips()
{
  expFaceShapeDerivs_.resize(numIntPoints_*nodesPerElement_*nDim_);
  quad9_shape_deriv(numIntPoints_, intgExpFace_.data(), expFaceShapeDerivs_.data());
}

//--------------------------------------------------------------------------
//-------- quad9_shape_fcn -------------------------------------------------
//--------------------------------------------------------------------------
void
QuadrilateralP2Element::quad9_shape_fcn(
  int  numIntPoints,
  const double *intgLoc,
  double *shpfc) const
{
  for ( int ip = 0; ip < numIntPoints; ++ip ) {
    int nineIp = nodesPerElement_ * ip; // nodes per element is always 9
    int vector_offset = nDim_ * ip;
    const double s = intgLoc[vector_offset+0];
    const double t = intgLoc[vector_offset+1];

    const double one_m_s = 1.0 - s;
    const double one_p_s = 1.0 + s;
    const double one_m_t = 1.0 - t;
    const double one_p_t = 1.0 + t;

    const double one_m_ss = 1.0 - s * s;
    const double one_m_tt = 1.0 - t * t;

    shpfc[nineIp  ] =  0.25 * s * t *  one_m_s *  one_m_t;
    shpfc[nineIp+1] = -0.25 * s * t *  one_p_s *  one_m_t;
    shpfc[nineIp+2] =  0.25 * s * t *  one_p_s *  one_p_t;
    shpfc[nineIp+3] = -0.25 * s * t *  one_m_s *  one_p_t;
    shpfc[nineIp+4] = -0.50 *     t *  one_p_s *  one_m_s * one_m_t;
    shpfc[nineIp+5] =  0.50 * s     *  one_p_t *  one_m_t * one_p_s;
    shpfc[nineIp+6] =  0.50 *     t *  one_p_s *  one_m_s * one_p_t;
    shpfc[nineIp+7] = -0.50 * s     *  one_p_t *  one_m_t * one_m_s;
    shpfc[nineIp+8] =  one_m_ss * one_m_tt;
  }
}

//--------------------------------------------------------------------------
//-------- quad9_shape_deriv -----------------------------------------------
//--------------------------------------------------------------------------
void
QuadrilateralP2Element::quad9_shape_deriv(
  int numIntPoints,
  const double *intgLoc,
  double *deriv) const
{
  for ( int ip = 0; ip < numIntPoints; ++ip ) {
    const int grad_offset = nDim_ * nodesPerElement_ * ip; // nodes per element is always 9
    const int vector_offset = nDim_ * ip;
    int node; int offset;

    const double s = intgLoc[vector_offset+0];
    const double t = intgLoc[vector_offset+1];

    const double s2 = s*s;
    const double t2 = t*t;

    node = 0;
    offset = grad_offset + nDim_ * node;
    deriv[offset+0] = 0.25 * (2.0 * s * t2 - 2.0 * s * t - t2 + t);
    deriv[offset+1] = 0.25 * (2.0 * s2 * t - 2.0 * s * t - s2 + s);

    node = 1;
    offset = grad_offset + nDim_ * node;
    deriv[offset+0] = 0.25 * (2.0 * s * t2 - 2.0 * s * t + t2 - t);
    deriv[offset+1] = 0.25 * (2.0 * s2 * t + 2.0 * s * t - s2 - s);

    node = 2;
    offset = grad_offset + nDim_ * node;
    deriv[offset+0] = 0.25 * (2.0 * s * t2 + 2.0 * s * t + t2 + t);
    deriv[offset+1] = 0.25 * (2.0 * s2 * t + 2.0 * s * t + s2 + s);

    node = 3;
    offset = grad_offset + nDim_ * node;
    deriv[offset+0] = 0.25 * (2.0 * s * t2 + 2.0 * s * t - t2 - t);
    deriv[offset+1] = 0.25 * (2.0 * s2 * t - 2.0 * s * t + s2 - s);

    node = 4;
    offset = grad_offset + nDim_ * node;
    deriv[offset+0] = -0.5 * (2.0 * s * t2 - 2.0 * s * t);
    deriv[offset+1] = -0.5 * (2.0 * s2 * t - s2 - 2.0 * t + 1.0);

    node = 5;
    offset = grad_offset + nDim_ * node;
    deriv[offset+0] = -0.5 * (2.0 * s * t2 + t2 - 2.0 * s - 1.0);
    deriv[offset+1] = -0.5 * (2.0 * s2 * t + 2.0 * s * t);

    node = 6;
    offset = grad_offset + nDim_ * node;
    deriv[offset+0] = -0.5 * (2.0 * s * t2 + 2.0 * s * t);
    deriv[offset+1] = -0.5 * (2.0 * s2 * t + s2 - 2.0 * t - 1.0);

    node = 7;
    offset = grad_offset + nDim_ * node;
    deriv[offset+0] = -0.5 * (2.0 * s * t2 - t2 - 2.0 * s + 1.0);
    deriv[offset+1] = -0.5 * (2.0 * s2 * t - 2.0 * s * t);

    node = 8;
    offset = grad_offset + nDim_ * node;
    deriv[offset+0] = 2.0 * s * t2 - 2.0 * s;
    deriv[offset+1] = 2.0 * s2 * t - 2.0 * t;
  }
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Quad92DSCV::Quad92DSCV()
: QuadrilateralP2Element()
{
  // set up the one-dimensional quadrature rule
  set_quadrature_rule();

  // set up integration rule and relevant maps for scvs
  set_interior_info();

  // compute and save shape functions and derivatives at ips
  eval_shape_functions_at_ips();
  eval_shape_derivs_at_ips();

  eval_shape_functions_at_shifted_ips();
  eval_shape_derivs_at_shifted_ips();
}

//--------------------------------------------------------------------------
//-------- set_interior_info -----------------------------------------------
//--------------------------------------------------------------------------
void
Quad92DSCV::set_interior_info()
{
  //1D integration rule per sub-control volume
  numIntPoints_ = (nodes1D_ * nodes1D_) * ( numQuad_ * numQuad_ ); // 36

  // define ip node mappings
  ipNodeMap_.resize(numIntPoints_);
  intgLoc_.resize(numIntPoints_*nDim_); // size = 72
  intgLocShift_.resize(numIntPoints_*nDim_); // size = 72
  ipWeight_.resize(numIntPoints_);

  // tensor product nodes (3x3x3) x tensor product quadrature (2x2x2)
  int vector_index = 0; int scalar_index = 0;
  for (int l = 0; l < nodes1D_; ++l) {
    for (int k = 0; k < nodes1D_; ++k) {
      const int nodeNumber = tensor_product_node_map(k,l);
      //tensor-product quadrature for a particular sub-cv
      for (int j = 0; j < numQuad_; ++j) {
        for (int i = 0; i < numQuad_; ++i) {
          //integration point location
          intgLoc_[vector_index]     = gauss_point_location(k,i);
          intgLoc_[vector_index + 1] = gauss_point_location(l,j);

          intgLocShift_[vector_index]     = shifted_gauss_point_location(k,i);
          intgLocShift_[vector_index + 1] = shifted_gauss_point_location(l,j);

          //weight
          ipWeight_[scalar_index] = tensor_product_weight(k,l,i,j);

          //sub-control volume association
          ipNodeMap_[scalar_index] = nodeNumber;

          // increment indices
          ++scalar_index;
          vector_index += nDim_;
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad92DSCV::ipNodeMap(
  int /*ordinal*/)
{
 // define scv->node mappings
 return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void Quad92DSCV::determinant(
  const int nelem,
  const double *coords,
  double *volume,
  double *error)
{
  for (int k = 0; k < nelem; ++k) {
    const int scalar_elem_offset = numIntPoints_ * k;
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;

      //weighted jacobian determinant
      const double det_j = jacobian_determinant(&coords[coord_elem_offset], &shapeDerivs_[grad_offset]);

      //apply weight and store to volume
      volume[scalar_elem_offset + ip] = ipWeight_[ip] * det_j;

      //flag error
      if (det_j <= 0.0) {
        *error = 1.0;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- jacobian_determinant --------------------------------------------
//--------------------------------------------------------------------------
double
Quad92DSCV::jacobian_determinant(
  const double *elemNodalCoords,
  const double *shapeDerivs) const
{
  double dx_ds1 = 0.0;  double dx_ds2 = 0.0;
  double dy_ds1 = 0.0;  double dy_ds2 = 0.0;

  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = node * nDim_;

    const double xCoord = elemNodalCoords[vector_offset + 0];
    const double yCoord = elemNodalCoords[vector_offset + 1];

    const double dn_ds1  = shapeDerivs[vector_offset + 0];
    const double dn_ds2  = shapeDerivs[vector_offset + 1];

    dx_ds1 += dn_ds1 * xCoord;
    dx_ds2 += dn_ds2 * xCoord;

    dy_ds1 += dn_ds1 * yCoord;
    dy_ds2 += dn_ds2 * yCoord;
  }

  const double det_j = dx_ds1 * dy_ds2 - dy_ds1 * dx_ds2;

  return det_j;
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Quad92DSCS::Quad92DSCS()
  : QuadrilateralP2Element()
{
  // set up the one-dimensional quadrature rule
  set_quadrature_rule();

  // set up integration rule and relevant maps for scs
  set_interior_info();

  // set up integration rule and relevant maps for faces
  set_boundary_info();

  // compute and save shape functions and derivatives at ips
  eval_shape_functions_at_ips();
  eval_shape_derivs_at_ips();
  eval_shape_derivs_at_face_ips();

  eval_shape_functions_at_shifted_ips();
  eval_shape_derivs_at_shifted_ips();
}

//--------------------------------------------------------------------------
//-------- set_interior_info -----------------------------------------------
//--------------------------------------------------------------------------
void
Quad92DSCS::set_interior_info()
{
  const int linesPerDirection = nodes1D_ - 1; // 2
  const int ipsPerLine = numQuad_ * nodes1D_;
  const int numLines = linesPerDirection * nDim_;

  numIntPoints_ = numLines * ipsPerLine; // 24

  // define L/R mappings
  lrscv_.resize(2*numIntPoints_); // size = 48

  // standard integration location
  intgLoc_.resize(numIntPoints_*nDim_); // size = 48

  // shifted
  intgLocShift_.resize(numIntPoints_*nDim_);

  ipInfo_.resize(numIntPoints_);

  // a list of the scs locations in 1D
  const std::vector<double> scsLoc =  { -scsDist_, scsDist_ };

  // correct orientation for area vector
  const std::vector<double> orientation = { -1.0, +1.0 };

  // specify integration point locations in a dimension-by-dimension manner

  //u-direction
  int vector_index = 0;
  int lrscv_index = 0;
  int scalar_index = 0;
  for (int m = 0; m < linesPerDirection; ++m) {
    for (int l = 0; l < nodes1D_; ++l) {

      int leftNode; int rightNode;
      if (m == 0) {
        leftNode  = tensor_product_node_map(l,m);
        rightNode = tensor_product_node_map(l,m + 1);
      }
      else {
        leftNode  = tensor_product_node_map(l,m + 1);
        rightNode = tensor_product_node_map(l,m);
      }

      for (int j = 0; j < numQuad_; ++j) {

        lrscv_[lrscv_index] = leftNode;
        lrscv_[lrscv_index + 1] = rightNode;

        intgLoc_[vector_index] = gauss_point_location(l,j);
        intgLoc_[vector_index + 1] = scsLoc[m];

        intgLocShift_[vector_index] = shifted_gauss_point_location(l,j);
        intgLocShift_[vector_index + 1] = scsLoc[m];

        //compute the quadrature weight
        ipInfo_[scalar_index].weight = orientation[m]*tensor_product_weight(l,j);

        //direction
        ipInfo_[scalar_index].direction = Jacobian::T_DIRECTION;

        ++scalar_index;
        lrscv_index += 2;
        vector_index += nDim_;
      }
    }
  }

  //t-direction
  for (int m = 0; m < linesPerDirection; ++m) {
    for (int l = 0; l < nodes1D_; ++l) {

      int leftNode; int rightNode;
      if (m == 0) {
        leftNode  = tensor_product_node_map(m,l);
        rightNode = tensor_product_node_map(m+1,l);
      }
      else {
        leftNode  = tensor_product_node_map(m+1,l);
        rightNode = tensor_product_node_map(m,l);
      }

      for (int j = 0; j < numQuad_; ++j) {

        lrscv_[lrscv_index]   = leftNode;
        lrscv_[lrscv_index+1] = rightNode;

        intgLoc_[vector_index] = scsLoc[m];
        intgLoc_[vector_index+1] = gauss_point_location(l,j);

        intgLocShift_[vector_index] = scsLoc[m];
        intgLocShift_[vector_index+1] = shifted_gauss_point_location(l,j);

        //compute the quadrature weight
        ipInfo_[scalar_index].weight = -orientation[m]*tensor_product_weight(l,j);

        //direction
        ipInfo_[scalar_index].direction = Jacobian::S_DIRECTION;

        ++scalar_index;
        lrscv_index += 2;
        vector_index += nDim_;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- set_boundary_info -----------------------------------------------
//--------------------------------------------------------------------------
void
Quad92DSCS::set_boundary_info()
{
  const int numFaces = 2*nDim_;
  const int nodesPerFace = nodes1D_;
  ipsPerFace_ = nodesPerFace*numQuad_;

  const int numFaceIps = numFaces*ipsPerFace_; // 24 -- different from numIntPoints_ for p > 2 ?

  oppFace_.resize(numFaceIps);
  ipNodeMap_.resize(numFaceIps);
  oppNode_.resize(numFaceIps);
  intgExpFace_.resize(numFaceIps*nDim_);

  const std::vector<int> stkFaceNodeMap = {
                                            0, 4, 1, //face 0, bottom face
                                            1, 5, 2, //face 1, right face
                                            2, 6, 3, //face 2, top face  -- reversed order
                                            3, 7, 0  //face 3, left face -- reversed order
                                          };

  auto face_node_number = [=] (int number,int faceOrdinal)
  {
    return stkFaceNodeMap[number+nodes1D_*faceOrdinal];
  };

  const std::vector<int> faceToLine = { 0, 3, 1, 2 };
  const std::vector<double> faceLoc = {-1.0, +1.0, +1.0, -1.0};

  int scalar_index = 0; int vector_index = 0;
  int faceOrdinal = 0; //bottom face
  int oppFaceIndex = 0;
  for (int k = 0; k < nodes1D_; ++k) {
    const int nearNode = face_node_number(k,faceOrdinal);
    int oppNode = tensor_product_node_map(k,1);

    for (int j = 0; j < numQuad_; ++j) {
      ipNodeMap_[scalar_index] = nearNode;
      oppNode_[scalar_index] = oppNode;
      oppFace_[scalar_index] = oppFaceIndex + faceToLine[faceOrdinal]*ipsPerFace_;

      intgExpFace_[vector_index]   = intgLoc_[oppFace_[scalar_index]*nDim_+0];
      intgExpFace_[vector_index+1] = faceLoc[faceOrdinal];

      ++scalar_index;
      vector_index += nDim_;
      ++oppFaceIndex;
    }
  }

  faceOrdinal = 1; //right face
  oppFaceIndex = 0;
  for (int k = 0; k < nodes1D_; ++k) {
    const int nearNode = face_node_number(k,faceOrdinal);
    int oppNode = tensor_product_node_map(1,k);

    for (int j = 0; j < numQuad_; ++j) {
      ipNodeMap_[scalar_index] = nearNode;
      oppNode_[scalar_index] = oppNode;
      oppFace_[scalar_index] = oppFaceIndex + faceToLine[faceOrdinal]*ipsPerFace_;

      intgExpFace_[vector_index]   = faceLoc[faceOrdinal];
      intgExpFace_[vector_index+1] = intgLoc_[oppFace_[scalar_index]*nDim_+1];

      ++scalar_index;
      vector_index += nDim_;
      ++oppFaceIndex;
    }
  }


  faceOrdinal = 2; //top face
  oppFaceIndex = 0;
  //NOTE: this face is reversed
  for (int k = nodes1D_-1; k >= 0; --k) {
    const int nearNode = face_node_number(nodes1D_-k-1,faceOrdinal);
    int oppNode = tensor_product_node_map(k,1);
    for (int j = 0; j < numQuad_; ++j) {
      ipNodeMap_[scalar_index] = nearNode;
      oppNode_[scalar_index] = oppNode;
      oppFace_[scalar_index] = (ipsPerFace_-1) - oppFaceIndex + faceToLine[faceOrdinal]*ipsPerFace_;

      intgExpFace_[vector_index] = intgLoc_[oppFace_[scalar_index]*nDim_+0];
      intgExpFace_[vector_index+1] = faceLoc[faceOrdinal];

      ++scalar_index;
      vector_index += nDim_;
      ++oppFaceIndex;
    }
  }

  faceOrdinal = 3; //left face
  oppFaceIndex = 0;
  //NOTE: this faces is reversed
  for (int k = nodes1D_-1; k >= 0; --k) {
    const int nearNode = face_node_number(nodes1D_-k-1,faceOrdinal);
    int oppNode = tensor_product_node_map(1,k);
    for (int j = 0; j < numQuad_; ++j) {
      ipNodeMap_[scalar_index] = nearNode;
      oppNode_[scalar_index] = oppNode;
      oppFace_[scalar_index] = (ipsPerFace_-1) - oppFaceIndex + faceToLine[faceOrdinal]*ipsPerFace_;

      intgExpFace_[vector_index]   = faceLoc[faceOrdinal];
      intgExpFace_[vector_index+1] = intgLoc_[oppFace_[scalar_index]*nDim_+1];

      ++scalar_index;
      vector_index += nDim_;
      ++oppFaceIndex;
    }
  }
}


//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad92DSCS::ipNodeMap(
  int ordinal)
{
  // define ip->node mappings for each face (ordinal); 
  return &ipNodeMap_[ordinal*ipsPerFace_];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void
Quad92DSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  //returns the normal vector (dyds,-dxds) for constant t curves
  //returns the normal vector (dydt,-dxdt) for constant s curves

  std::array<double,2> areaVector;

  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    const int vector_elem_offset = nDim_*numIntPoints_*k;

    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;
      const int offset = nDim_ * ip + vector_elem_offset;

      //compute area vector for this ip
      area_vector( ipInfo_[ip].direction,
                   &coords[coord_elem_offset],
                   &shapeDerivs_[grad_offset],
                   areaVector.data() );

      // apply quadrature weight and orientation (combined as weight)
      for (int j = 0; j < nDim_; ++j) {
        areav[offset+j]  = ipInfo_[ip].weight * areaVector[j];
      }
    }
  }

  *error = 0;
}

//--------------------------------------------------------------------------
//-------- grad_op ---------------------------------------------------------
//--------------------------------------------------------------------------
void Quad92DSCS::grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  int numShapeDerivs = numIntPoints_*nodesPerElement_*nDim_;
  for (int j = 0; j < numShapeDerivs; ++j) {
    deriv[j] = shapeDerivs_[j];
  }

  SIERRA_FORTRAN(quad_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );

  if ( lerr )
    std::cout << "sorry, negative area.." << std::endl;

}

//--------------------------------------------------------------------------
//-------- shifted_grad_op -------------------------------------------------
//--------------------------------------------------------------------------
void Quad92DSCS::shifted_grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  int numShapeDerivs = numIntPoints_*nodesPerElement_*nDim_;
  for (int j = 0; j < numShapeDerivs; ++j) {
    deriv[j] = shapeDerivsShift_[j];
  }

  SIERRA_FORTRAN(quad_gradient_operator)
  ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );

  if ( lerr )
    std::cout << "sorry, negative area.." << std::endl;
}

//--------------------------------------------------------------------------
//-------- face_grad_op ----------------------------------------------------
//--------------------------------------------------------------------------
void Quad92DSCS::face_grad_op(
  const int nelem,
  const int face_ordinal,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  int lerr = 0;

  const int nface = 1;
  const int face_offset =  nDim_ * ipsPerFace_ * nodesPerElement_ * face_ordinal;
  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    const int scalar_elem_offset = ipsPerFace_ * k;

    for (int ip = 0; ip < ipsPerFace_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;

      SIERRA_FORTRAN(quad_gradient_operator)
      ( &nface,
          &nodesPerElement_,
          &nface,
          &expFaceShapeDerivs_[face_offset+grad_offset],
          &coords[coord_elem_offset],
          &gradop[grad_offset],
          &det_j[scalar_elem_offset+ip],
          error,
          &lerr
      );

      if (det_j[ip] <= 0.0) {
        *error = 1.0;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- gij -------------------------------------------------------------
//--------------------------------------------------------------------------
void Quad92DSCS::gij(
  const double *coords,
  double *gupperij,
  double *glowerij,
  double *deriv)
{
  SIERRA_FORTRAN(twod_gij)
    ( &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gupperij, glowerij);
}

//--------------------------------------------------------------------------
//-------- adjacentNodes ---------------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad92DSCS::adjacentNodes()
{
  // define L/R mappings
  return &lrscv_[0];
}

//--------------------------------------------------------------------------
//-------- opposingNodes ---------------------------------------------------
//--------------------------------------------------------------------------
int
Quad92DSCS::opposingNodes(
  const int ordinal,
  const int node)
{
  return oppNode_[ordinal*ipsPerFace_+node];
}

//--------------------------------------------------------------------------
//-------- opposingFace ----------------------------------------------------
//--------------------------------------------------------------------------
int
Quad92DSCS::opposingFace(
  const int ordinal,
  const int node)
{
  return oppFace_[ordinal*ipsPerFace_+node];
}

//--------------------------------------------------------------------------
//-------- area_vector -----------------------------------------------------
//--------------------------------------------------------------------------
void
Quad92DSCS::area_vector(
  const Jacobian::Direction direction,
  const double *elemNodalCoords,
  double *shapeDeriv,
  double *normalVec ) const
{
  int s1Component;
  switch (direction) {
    case Jacobian::S_DIRECTION:
      s1Component = static_cast<int>(Jacobian::T_DIRECTION);
      break;
    case Jacobian::T_DIRECTION:
      s1Component = static_cast<int>(Jacobian::S_DIRECTION);
      break;
    default:
      throw std::runtime_error("Not a valid direction for this element!");
  }

  double dxdr = 0.0;  double dydr = 0.0;
  for (int node = 0; node < nodesPerElement_; ++node) {
    const int vector_offset = nDim_ * node;
    const double xCoord = elemNodalCoords[vector_offset+0];
    const double yCoord = elemNodalCoords[vector_offset+1];

    dxdr += shapeDeriv[vector_offset+s1Component] * xCoord;
    dydr += shapeDeriv[vector_offset+s1Component] * yCoord;
  }

  normalVec[0] =  dydr;
  normalVec[1] = -dxdr;
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Tri2DSCV::Tri2DSCV()
  : MasterElement()
{
  nDim_ = 2;
  nodesPerElement_ = 3;
  numIntPoints_ = 3;

  // define ip node mappings
  ipNodeMap_.resize(3);
  ipNodeMap_[0] = 0; ipNodeMap_[1] = 1; ipNodeMap_[2] = 2;

  intgLoc_ = {
      5.0/24.0, 5.0/24.0,
      7.0/12.0, 5.0/24.0,
      5.0/24.0, 7.0/12.0
  };

  intgLocShift_ = {
      0.0,  0.0,
      1.0,  0.0,
      0.0,  1.0
  };
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Tri2DSCV::~Tri2DSCV()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Tri2DSCV::ipNodeMap(
  int /*ordinal*/)
{
  // define scv->node mappings
  return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Tri2DSCV::shape_fcn(double *shpfc)
{
  tri_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Tri2DSCV::shifted_shape_fcn(double *shpfc)
{
  tri_shape_fcn(numIntPoints_, &intgLocShift_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- tri_shape_fcn ---------------------------------------------------
//--------------------------------------------------------------------------
void
Tri2DSCV::tri_shape_fcn(
  const int  &npts,
  const double *isoParCoord,
  double *shape_fcn)
{
  for (int j = 0; j < npts; ++j ) {
    const int threej = 3*j;
    const int k = 2*j;
    const double xi = isoParCoord[k];
    const double eta = isoParCoord[k+1];
    shape_fcn[threej] = 1.0 - xi - eta;
    shape_fcn[1 + threej] = xi;
    shape_fcn[2 + threej] = eta;
  }
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void Tri2DSCV::determinant(
  const int nelem,
  const double *coords,
  double *volume,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(tri_scv_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords,
      volume, error, &lerr );
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Tri2DSCS::Tri2DSCS()
  : MasterElement()
{
  nDim_ = 2;
  nodesPerElement_ = 3;
  numIntPoints_ = 3;

  // define L/R mappings
  lrscv_.resize(6);
  lrscv_[0]  = 0; lrscv_[1]  = 1;
  lrscv_[2]  = 1; lrscv_[3]  = 2;
  lrscv_[4]  = 0; lrscv_[5]  = 2;

  // define opposing node
  oppNode_.resize(6);
  // face 0; nodes 0,1
  oppNode_[0] = 2; oppNode_[1] = 2;
  // face 1; nodes 1,2
  oppNode_[2] = 0; oppNode_[3] = 0;
  // face 2; nodes 2,0
  oppNode_[4] = 1; oppNode_[5] = 1;

  // define opposing face
  oppFace_.resize(6);
  // face 0
  oppFace_[0]  = 2; oppFace_[1] = 1;
  // face 1
  oppFace_[2]  = 0; oppFace_[3] = 2;
  // face 2
  oppFace_[4]  = 1; oppFace_[5] = 0;  

  // standard integration location
  const double five12ths = 5.0/12.0;
  const double one6th = 1.0/6.0;
  intgLoc_.resize(6);    
  intgLoc_[0] = five12ths; intgLoc_[1] = one6th;    // surf 1; 0->1
  intgLoc_[2] = five12ths; intgLoc_[3] = five12ths; // surf 2; 1->3
  intgLoc_[4] = one6th;    intgLoc_[5] = five12ths; // surf 3; 0->2

  // shifted
  intgLocShift_.resize(6);
  intgLocShift_[0] = 0.50; intgLocShift_[1] = 0.00;  // surf 1; 0->1
  intgLocShift_[2] = 0.50; intgLocShift_[3] = 0.50;  // surf 1; 1->3
  intgLocShift_[4] = 0.00; intgLocShift_[5] = 0.50;  // surf 1; 0->2

  // exposed face
  intgExpFace_.resize(12);
  // face 0; scs 0, 1; nodes 0,1
  intgExpFace_[0]  = 0.25; intgExpFace_[1]  = 0.00; 
  intgExpFace_[2]  = 0.75; intgExpFace_[3]  = 0.00;
  // face 1; scs 0, 1; nodes 1,2
  intgExpFace_[4]  = 0.75; intgExpFace_[5]  = 0.25;
  intgExpFace_[6]  = 0.25; intgExpFace_[7]  = 0.75;
  // face 2, surf 0, 1; nodes 2,0
  intgExpFace_[8]  = 0.00; intgExpFace_[9]  = 0.75;
  intgExpFace_[10] = 0.00; intgExpFace_[11] = 0.25;
  
  // boundary integration point ip node mapping (ip on an ordinal to local node number)
  ipNodeMap_.resize(6); // 2 ips * 3 faces
  // face 0;
  ipNodeMap_[0] = 0;  ipNodeMap_[1] = 1; 
  // face 1; 
  ipNodeMap_[2] = 1;  ipNodeMap_[3] = 2; 
  // face 2;
  ipNodeMap_[4] = 2;  ipNodeMap_[5] = 0;  
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Tri2DSCS::~Tri2DSCS()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Tri2DSCS::ipNodeMap(
  int ordinal)
{
  // define ip->node mappings for each face (ordinal); 
  return &ipNodeMap_[ordinal*2];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void Tri2DSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  SIERRA_FORTRAN(tri_scs_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_, coords, areav );

  // all is always well; no error checking
  *error = 0;
}

//--------------------------------------------------------------------------
//-------- grad_op ---------------------------------------------------------
//--------------------------------------------------------------------------
void Tri2DSCS::grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(tri_derivative)
    ( &numIntPoints_, deriv );
  
  SIERRA_FORTRAN(tri_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );
  
  if ( lerr )
    std::cout << "sorry, negative Tri2DSCS volume.." << std::endl;
}

//--------------------------------------------------------------------------
//-------- shifted_grad_op -------------------------------------------------
//--------------------------------------------------------------------------
void Tri2DSCS::shifted_grad_op(
  const int nelem,
  const double *coords,
  double *gradop,
  double *deriv,
  double *det_j,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(tri_derivative)
    ( &numIntPoints_, deriv );

  SIERRA_FORTRAN(tri_gradient_operator)
    ( &nelem,
      &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gradop, det_j, error, &lerr );

  if ( lerr )
    std::cout << "sorry, negative Tri2DSCS volume.." << std::endl;
}

//--------------------------------------------------------------------------
//-------- face_grad_op ----------------------------------------------------
//--------------------------------------------------------------------------
void Tri2DSCS::face_grad_op(
  const int nelem,
  const int /*face_ordinal*/,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  int lerr = 0;
  int npf = 2;

  const int nface = 1;
  double dpsi[6];
  double grad[6];

  for ( int n=0; n<nelem; n++ ) {
    
    for ( int k=0; k<npf; k++ ) {
      
      // derivatives are constant
      SIERRA_FORTRAN(tri_derivative)
        ( &nface, dpsi );
      
      SIERRA_FORTRAN(tri_gradient_operator)
        ( &nface,
          &nodesPerElement_,
          &nface,
          dpsi,
          &coords[12*n], grad, &det_j[npf*n+k], error, &lerr );
      
      if ( lerr )
        std::cout << "sorry, issue with face_grad_op.." << std::endl;
      
      for ( int j=0; j<6; j++) {
        gradop[k*nelem*6+n*6+j] = grad[j];
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- gij -------------------------------------------------------------
//--------------------------------------------------------------------------
void Tri2DSCS::gij(
  const double *coords,
  double *gupperij,
  double *glowerij,
  double *deriv)
{
  SIERRA_FORTRAN(twod_gij)
    ( &nodesPerElement_,
      &numIntPoints_,
      deriv,
      coords, gupperij, glowerij);
}

//--------------------------------------------------------------------------
//-------- adjacentNodes ---------------------------------------------------
//--------------------------------------------------------------------------
const int *
Tri2DSCS::adjacentNodes()
{
  // define L/R mappings
  return &lrscv_[0];
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Tri2DSCS::shape_fcn(double *shpfc)
{
  tri_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Tri2DSCS::shifted_shape_fcn(double *shpfc)
{
  tri_shape_fcn(numIntPoints_, &intgLocShift_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- tri_shape_fcn ---------------------------------------------------
//--------------------------------------------------------------------------
void
Tri2DSCS::tri_shape_fcn(
  const int  &npts,
  const double *isoParCoord, 
  double *shape_fcn)
{
  for (int j = 0; j < npts; ++j ) {
    const int threej = 3*j;
    const int k = 2*j;
    const double xi = isoParCoord[k];
    const double eta = isoParCoord[k+1];
    shape_fcn[threej] = 1.0 - xi - eta;
    shape_fcn[1 + threej] = xi;
    shape_fcn[2 + threej] = eta;
  }
}

//--------------------------------------------------------------------------
//-------- opposingNodes --------------------------------------------------
//--------------------------------------------------------------------------
int
Tri2DSCS::opposingNodes(
  const int ordinal,
  const int node)
{
  return oppNode_[ordinal*2+node];
}

//--------------------------------------------------------------------------
//-------- opposingFace --------------------------------------------------
//--------------------------------------------------------------------------
int
Tri2DSCS::opposingFace(
  const int ordinal,
  const int node)
{
  return oppFace_[ordinal*2+node];
}

//--------------------------------------------------------------------------
//-------- isInElement -----------------------------------------------------
//--------------------------------------------------------------------------
double
Tri2DSCS::isInElement(
  const double *elemNodalCoord,
  const double *pointCoord,
  double *isoParCoord )
{
  // Translate element so that (x,y) coordinates of the
  // first node

  double x[2] = { elemNodalCoord[1] - elemNodalCoord[0],
		elemNodalCoord[2] - elemNodalCoord[0] };
  double y[2] = { elemNodalCoord[4] - elemNodalCoord[3],
		elemNodalCoord[5] - elemNodalCoord[3] };

  // Translate position vector of point in same manner

  double xp = pointCoord[0] - elemNodalCoord[0];
  double yp = pointCoord[1] - elemNodalCoord[3];

  // Set new nodal coordinates with Node 1 at origin and with new
  // x and y axes lying in the plane of the element
  double len12 = std::sqrt( x[0]*x[0] + y[0]*y[0] );
  double len13 = std::sqrt( x[1]*x[1] + y[1]*y[1] );

  double xnew[2];
  double ynew[2];

  // Use cross-product of find enclosed angle 

  const double cross = x[0]*y[1] - x[1]*y[0];

  double Area2 = std::sqrt( cross*cross );

  // find sin of angle
  double sin_theta = Area2 / ( len12 * len13 ) ;

  // find cosine of angle
  double cos_theta = (x[0]*x[1] + y[0]*y[1])/(len12 * len13);

  // nodal coordinates of nodes 2 and 3 in new system
  // (coordinates of node 1 are identically 0.0)
  double x_nod_new[2] = { len12, len13*cos_theta};
  double y_nod_new[2] = {  0.0, len13*sin_theta};

  // find direction cosines transform position vector of
  // point to be checked into new coordinate system
  // direction cosines of new x axis along side 12

  xnew[0] = x[0]/len12;
  xnew[1] = y[0]/len12;

  // direction cosines of new y-axis
  ynew[0] =  -xnew[1];
  ynew[1] =   xnew[0];

  // compute transformed coordinates of point
  // (coordinates in xnew,ynew)
  double xpnew = xnew[0]*xp + xnew[1]*yp;
  double ypnew = ynew[0]*xp + ynew[1]*yp;

  // Find parametric coordinates of point and check that
  // it lies in the element
  isoParCoord[0] = 1. - xpnew / x_nod_new[0] +
    ypnew*( x_nod_new[1] - x_nod_new[0] ) / Area2;
  isoParCoord[1] = ( xpnew*y_nod_new[1] - ypnew*x_nod_new[1] ) / Area2;

  std::vector<double> w(2);
  w[0]=isoParCoord[0];
  w[1]=isoParCoord[1];

  isoParCoord[0] = w[1];
  isoParCoord[1] = 1.0-w[0]-w[1];

  const double dist = tri_parametric_distance(w);

  return dist;

}

//--------------------------------------------------------------------------
//-------- tri_parametric_distance -----------------------------------------
//--------------------------------------------------------------------------
double
Tri2DSCS::tri_parametric_distance(
  const std::vector<double> &x)
{
  const double X=x[0] - 1./3.;
  const double Y=x[1] - 1./3.;
  const double dist0 = -3*X;
  const double dist1 = -3*Y;
  const double dist2 =  3*(X+Y);
  const double dist = std::max(std::max(dist0,dist1),dist2);
  return dist;
}

//--------------------------------------------------------------------------
//-------- interpolatePoint ------------------------------------------------
//--------------------------------------------------------------------------
void
Tri2DSCS::interpolatePoint(
  const int &nComp,
  const double *isoParCoord,
  const double *field,
  double *result )
{
  const double s = isoParCoord[0];
  const double t = isoParCoord[1];
  const double oneMinusST = 1.0 - s - t;
  for ( int i = 0; i < nComp; i++ )
  {
    const int b = 3*i;
    result[i] =  oneMinusST*field[b+0] + s*field[b+1] + t*field[b+2];
  }
}

//--------------------------------------------------------------------------
//-------- general_face_grad_op --------------------------------------------
//--------------------------------------------------------------------------
void 
Tri2DSCS::general_face_grad_op(
  const int /*face_ordinal*/,
  const double */*isoParCoord*/,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  int lerr = 0;
  
  const int nface = 1;
  double dpsi[6];
  double grad[6];
  
  // derivatives are constant
  SIERRA_FORTRAN(tri_derivative)
    ( &nface, dpsi );
      
  SIERRA_FORTRAN(tri_gradient_operator)
    ( &nface,
      &nodesPerElement_,
      &nface,
      dpsi,
      &coords[0], grad, &det_j[0], error, &lerr );
      
  if ( lerr )
    std::cout << "sorry, issue with face_grad_op.." << std::endl;
  
  for ( int j=0; j<6; j++) {
    gradop[j] = grad[j];
  }
}


//--------------------------------------------------------------------------
//-------- sidePcoords_to_elemPcoords --------------------------------------
//--------------------------------------------------------------------------
void 
Tri2DSCS::sidePcoords_to_elemPcoords(
  const int & side_ordinal,
  const int & npoints,
  const double *side_pcoords,
  double *elem_pcoords)
{
  switch (side_ordinal) {
  case 0:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*2+0] = 0.5*(1.0 + side_pcoords[i]);
      elem_pcoords[i*2+1] = 0.0;
    }
    break;
  case 1:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*2+0] = 1. - 0.5 * (side_pcoords[i] + 1.);
      elem_pcoords[i*2+1] = 0.5 * (side_pcoords[i] + 1.);
    }
    break;
  case 2:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*2+0] = 0.0;
      elem_pcoords[i*2+1] = 1. - 0.5 * (side_pcoords[i] + 1.);
    }
    break;
  default:
    throw std::runtime_error("Tri2DSCS::sideMap invalid ordinal");
  }
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Quad3DSCS::Quad3DSCS()  
  : MasterElement(),
    elemThickness_(0.1)
{
  nDim_ = 3;
  nodesPerElement_ = 4;
  numIntPoints_ = 4;
  scaleToStandardIsoFac_ = 2.0;

  // define ip node mappings; ordinal size = 1
  ipNodeMap_.resize(4);
  ipNodeMap_[0] = 0;
  ipNodeMap_[1] = 1;
  ipNodeMap_[2] = 2;
  ipNodeMap_[3] = 3;

  // standard integration location
  intgLoc_.resize(8);    
  intgLoc_[0]  = -0.25; intgLoc_[1] = -0.25; // surf 1
  intgLoc_[2]  =  0.25; intgLoc_[3] = -0.25; // surf 2
  intgLoc_[4]  =  0.25; intgLoc_[5] =  0.25; // surf 3
  intgLoc_[6]  = -0.25; intgLoc_[7] =  0.25; // surf 4

  // shifted
  intgLocShift_.resize(8);    
  intgLocShift_[0]  = -0.50; intgLocShift_[1] = -0.50; // surf 1
  intgLocShift_[2]  =  0.50; intgLocShift_[3] = -0.50; // surf 2
  intgLocShift_[4]  =  0.50; intgLocShift_[5] =  0.50; // surf 3
  intgLocShift_[6]  = -0.50; intgLocShift_[7] =  0.50; // surf 4  
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Quad3DSCS::~Quad3DSCS()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad3DSCS::ipNodeMap(
  int /*ordinal*/)
{
  // define ip->node mappings for each face (single ordinal); 
  return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void Quad3DSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(quad3d_scs_det)
    ( &nelem, coords, areav );

  // fake check
  *error = (lerr == 0) ? 0.0 : 1.0;
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Quad3DSCS::shape_fcn(double *shpfc)
{
  SIERRA_FORTRAN(quad3d_shape_fcn)
    (&numIntPoints_,&intgLoc_[0],shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Quad3DSCS::shifted_shape_fcn(double *shpfc)
{
  SIERRA_FORTRAN(quad3d_shape_fcn)
    (&numIntPoints_,&intgLocShift_[0],shpfc);
}


//--------------------------------------------------------------------------
//-------- isInElement -----------------------------------------------------
//--------------------------------------------------------------------------
double
Quad3DSCS::isInElement(
  const double *elemNodalCoord,
  const double *pointCoord,
  double *isoParCoord )
{
  // square of the desired norm, 1.0e-8
  const double isInElemConverged = 1.0e-16;
  const int maxNonlinearIter = 20;

  // Translate element so that (x,y,z) coordinates of the first node are (0,0,0)

  double x[3] = { elemNodalCoord[1] - elemNodalCoord[0],
		elemNodalCoord[2] - elemNodalCoord[0],
		elemNodalCoord[3] - elemNodalCoord[0] };

  double y[3] = { elemNodalCoord[5] - elemNodalCoord[4],
		elemNodalCoord[6] - elemNodalCoord[4],
		elemNodalCoord[7] - elemNodalCoord[4] };

  double z[3] = { elemNodalCoord[9]  - elemNodalCoord[8],
		elemNodalCoord[10] - elemNodalCoord[8],
		elemNodalCoord[11] - elemNodalCoord[8] };

  // (xp,yp,zp) is the point at which we're searching for (xi,eta,d)
  // (must translate this also)
  // d = (scaled) distance in (x,y,z) space from point (xp,yp,zp) to the
  //     surface defined by the face element (the distance is scaled by
  //     the length of the non-unit normal vector; rescaling of d is done
  //     following the NR iteration below).

  double xp = pointCoord[0] - elemNodalCoord[0];
  double yp = pointCoord[1] - elemNodalCoord[4];
  double zp = pointCoord[2] - elemNodalCoord[8];


  // Newton-Raphson iteration for (xi,eta,d)

  double jdet;
  double j[9];
  double gn[3];
  double xcur[3];          // current (x,y,z) point on element surface
  double normal[3];        // (non-unit) normal computed at xcur

  // Solution vector solcur[3] = {xi,eta,d}
  double solcur[3] = {-0.5,-0.5,-0.5};     // initial guess
  double deltasol[] = {1.0,1.0, 1.0};

  int i = 0;
  do
  {
    // Update guess vector
    solcur[0] += deltasol[0];
    solcur[1] += deltasol[1];
    solcur[2] += deltasol[2];

    interpolatePoint(3,solcur,elemNodalCoord,xcur);

    // Translate xcur ((x,y,z) point corresponding
    // to current (xi,eta) guess)

    xcur[0] -= elemNodalCoord[0];
    xcur[1] -= elemNodalCoord[4];
    xcur[2] -= elemNodalCoord[8];

    non_unit_face_normal(solcur,elemNodalCoord,normal);

    gn[0] = xcur[0] - xp + solcur[2] * normal[0];
    gn[1] = xcur[1] - yp + solcur[2] * normal[1];
    gn[2] = xcur[2] - zp + solcur[2] * normal[2];

    // Mathematica-generated code for the jacobian

    j[0]=0.125*(-2.00*(-1.00+solcur[1])*x[0]
			    +(2.00*(1.00+solcur[1])*(x[1]-x[2])+solcur[2]
			      *(-(y[1]*z[0])+y[2]*z[0]+y[0]*z[1]-y[0]*z[2])));

    j[1]=0.125*(-2.00*(1.00+solcur[0])*x[0]
			    +2.00*(1.00+solcur[0])*x[1]-2.00
			    *(-1.00+solcur[0])*x[2]+(solcur[2]*(y[2]*(z[0]-z[1])+(-y[0]+y[1])*z[2])));

    j[2]= normal[0];

    j[3]=0.125*(-2.00*(-1.00+solcur[1])*y[0]
			    +(2.00*(1.00+solcur[1])*(y[1]-y[2])
			      +solcur[2]*(x[1]*z[0]-x[2]*z[0]-x[0]*z[1]+x[0]*z[2])));

    j[4]=0.125*(-2.00*(1.00+solcur[0])*y[0]
			    +2.00*(1.00+solcur[0])*y[1]
			    -2.00*(-1.00+solcur[0])*y[2]+(solcur[2]*(x[2]*(-z[0]+z[1])+(x[0]-x[1])*z[2])));

    j[5]= normal[1];

    j[6]=0.125*((solcur[2]*(-(x[1]*y[0])+x[2]*y[0]+x[0]*y[1]-x[0]*y[2]))
			    -2.00*((-1.00+solcur[1])*z[0]
					       -(1.00+solcur[1])*(z[1]-z[2])));

    j[7]=0.125*((solcur[2]*(x[2]*(y[0]-y[1])+(-x[0]+x[1])*y[2]))
			    -2.00*(1.00+solcur[0])*z[0]+2.00
			    *(1.00+solcur[0])*z[1]-2.00*(-1.00+solcur[0])*z[2]);
    
    j[8]= normal[2];
    
    jdet=-(j[2]*j[4]*j[6])+j[1]*j[5]*j[6]+j[2]*j[3]*j[7]-
      j[0]*j[5]*j[7]-j[1]*j[3]*j[8]+j[0]*j[4]*j[8];


    // Solve linear system (j*deltasol = -gn) for deltasol at step n+1
    
    deltasol[0] = (gn[2]*(j[2]*j[4]-j[1]*j[5])+gn[1]*(-(j[2]*j[7])+
		  j[1]*j[8])+gn[0]*(j[5]*j[7]-j[4]*j[8]))/jdet;
    deltasol[1] = (gn[2]*(-(j[2]*j[3])+j[0]*j[5])+gn[1]*(j[2]*j[6]-
		  j[0]*j[8])+gn[0]*(-(j[5]*j[6])+j[3]*j[8]))/jdet;
    deltasol[2] = (gn[2]*(j[1]*j[3]-j[0]*j[4])+gn[1]*(-(j[1]*j[6])+
		  j[0]*j[7])+gn[0]*(j[4]*j[6]-j[3]*j[7]))/jdet;

  } while ( !within_tolerance( vector_norm_sq(deltasol,3), isInElemConverged)
	    && ++i < maxNonlinearIter );

  // Fill in solution vector; only include the distance (in the third
  // solution vector slot) if npar_coord = 3 (this is how the user
  // requests it)

  isoParCoord[0] = std::numeric_limits<double>::max();
  isoParCoord[1] = std::numeric_limits<double>::max();
  isoParCoord[2] = std::numeric_limits<double>::max();
  double dist = std::numeric_limits<double>::max();

  if (i < maxNonlinearIter) {
    isoParCoord[0] = solcur[0] + deltasol[0];
    isoParCoord[1] = solcur[1] + deltasol[1];
    // Rescale the distance vector by the length of the (non-unit) normal vector,
    // which was used above in the NR iteration.
    const double area   = std::sqrt(vector_norm_sq(normal,3));
    const double length = std::sqrt(area);
    
    const double par_coor_2 = (solcur[2] + deltasol[2]) * length;
    //if ( npar_coord == 3 ) par_coor[2] = par_coor_2;
    isoParCoord[2] = par_coor_2;

    std::vector<double> xtmp(3);
    xtmp[0] = isoParCoord[0];
    xtmp[1] = isoParCoord[1];
    xtmp[2] = isoParCoord[2];
    dist = parametric_distance(xtmp);
  }
  return dist;
}

void
Quad3DSCS::non_unit_face_normal(
  const double * isoParCoord,            // (2)
  const double * elem_nodal_coor,        // (4,3)
	double * normal_vector )         // (3)
{
  double xi  = isoParCoord[0];
  double eta = isoParCoord[1];

  // Translate element so that node 0 is at (x,y,z) = (0,0,0)

  double x[3] = { elem_nodal_coor[1] - elem_nodal_coor[0],
                  elem_nodal_coor[2] - elem_nodal_coor[0],
                  elem_nodal_coor[3] - elem_nodal_coor[0] };
  
  double y[3] = { elem_nodal_coor[5] - elem_nodal_coor[4],
                  elem_nodal_coor[6] - elem_nodal_coor[4],
                  elem_nodal_coor[7] - elem_nodal_coor[4] };
  
  double z[3] = { elem_nodal_coor[9]  - elem_nodal_coor[8],
                  elem_nodal_coor[10] - elem_nodal_coor[8],
                  elem_nodal_coor[11] - elem_nodal_coor[8] };
  
  // Mathematica-generated and simplified code for the normal vector

  double n0 = 0.125*(xi*y[2]*z[0]+y[0]*z[1]+xi*y[0]*z[1]-y[2]*z[1]-
			       xi*y[0]*z[2]+y[1]*(-((1.00+xi)*z[0])+
	   (1.00+eta)*z[2])+eta*(y[2]*z[0]-y[2]*z[1]-y[0]*z[2]));

  double n1 = 0.125*(-(xi*x[2]*z[0])-x[0]*z[1]-xi*x[0]*z[1]+x[2]*z[1]+
				 xi*x[0]*z[2]+x[1]*((1.00+xi)*z[0]-
	   (1.00+eta)*z[2])+eta*(-(x[2]*z[0])+x[2]*z[1]+x[0]*z[2]));

  double n2 = 0.125*(xi*x[2]*y[0]+x[0]*y[1]+xi*x[0]*y[1]-x[2]*y[1]-
			       xi*x[0]*y[2]+x[1]*(-((1.00+xi)*y[0])+
	   (1.00+eta)*y[2])+eta*(x[2]*y[0]-x[2]*y[1]-x[0]*y[2]));

  normal_vector[0] = n0;
  normal_vector[1] = n1;
  normal_vector[2] = n2;

}

double 
Quad3DSCS::parametric_distance(const std::vector<double> &x)
{
  const int NCOORD   = 3;
  std::vector<double> y(NCOORD);
  
  for (int i=0; i<NCOORD; ++i) {
    y[i] = std::abs(x[i]);
  }

  double d = y[0];
  if (d < y[1]) d = y[1];
  if (elemThickness_ < y[2] && d < 1+y[2]) d = 1+y[2];
  return d;
}

//--------------------------------------------------------------------------
//-------- interpolatePoint ------------------------------------------------
//--------------------------------------------------------------------------
void
Quad3DSCS::interpolatePoint(
  const int &nComp,
  const double *isoParCoord,
  const double *field,
  double *result )
{
  // this is the same as the 2D implementation... Consider consolidation
  const double xi   = isoParCoord[0];
  const double eta  = isoParCoord[1];

  for ( int i = 0; i < nComp; i++ )
  {
    // Base 'field array' index for ith component
    int b = 4*i;

    result[i] = 0.250 * (
      (1.00-eta) * (1.00-xi ) * field[b+0] +
      (1.00-eta) * (1.00+xi ) * field[b+1] +
      (1.00+eta) * (1.00+xi ) * field[b+2] +
      (1.00+eta) * (1.00-xi ) * field[b+3] ) ;
  }
}

//--------------------------------------------------------------------------
//-------- general_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Quad3DSCS::general_shape_fcn(
  const int numIp,
  const double *isoParCoord,
  double *shpfc)
{
  // -1:1 isoparametric range
  const double npe = nodesPerElement_;
  for ( int ip = 0; ip < numIp; ++ip ) {
    
    const int rowIpc = 2*ip;
    const int rowSfc = npe*ip;
    
    const double s1 = isoParCoord[rowIpc];
    const double s2 = isoParCoord[rowIpc+1];
    shpfc[rowSfc  ] = 0.25*(1.0-s1)*(1.0-s2);
    shpfc[rowSfc+1] = 0.25*(1.0+s1)*(1.0-s2);
    shpfc[rowSfc+2] = 0.25*(1.0+s1)*(1.0+s2);
    shpfc[rowSfc+3] = 0.25*(1.0-s1)*(1.0+s2);
 
  }
}

//--------------------------------------------------------------------------
//-------- general_normal --------------------------------------------------
//--------------------------------------------------------------------------
void
Quad3DSCS::general_normal(
  const double *isoParCoord,
  const double *coords,
  double *normal)
{
  const int nDim = 3;

  const double psi0Xi = -0.25 * (1.0 - isoParCoord[1]);
  const double psi1Xi =  0.25 * (1.0 - isoParCoord[1]);
  const double psi2Xi =  0.25 * (1.0 + isoParCoord[1]);
  const double psi3Xi = -0.25 * (1.0 + isoParCoord[1]);
  
  const double psi0Eta =-0.25 * (1.0 - isoParCoord[0]);
  const double psi1Eta =-0.25 * (1.0 + isoParCoord[0]);
  const double psi2Eta = 0.25 * (1.0 + isoParCoord[0]);
  const double psi3Eta = 0.25 * (1.0 - isoParCoord[0]);
  
  const double DxDxi = coords[0*nDim+0]*psi0Xi +
    coords[1*nDim+0]*psi1Xi +
    coords[2*nDim+0]*psi2Xi +
    coords[3*nDim+0]*psi3Xi;  
      
  const double DyDxi = coords[0*nDim+1]*psi0Xi +
    coords[1*nDim+1]*psi1Xi +
    coords[2*nDim+1]*psi2Xi +
    coords[3*nDim+1]*psi3Xi;
  
  const double DzDxi = coords[0*nDim+2]*psi0Xi +
    coords[1*nDim+2]*psi1Xi +
    coords[2*nDim+2]*psi2Xi +
    coords[3*nDim+2]*psi3Xi;
  
  const double DxDeta = coords[0*nDim+0]*psi0Eta +
    coords[1*nDim+0]*psi1Eta +
    coords[2*nDim+0]*psi2Eta +
    coords[3*nDim+0]*psi3Eta;

  const double DyDeta = coords[0*nDim+1]*psi0Eta +
    coords[1*nDim+1]*psi1Eta +
    coords[2*nDim+1]*psi2Eta +
    coords[3*nDim+1]*psi3Eta;
  
  const double DzDeta = coords[0*nDim+2]*psi0Eta +
    coords[1*nDim+2]*psi1Eta +
    coords[2*nDim+2]*psi2Eta +
    coords[3*nDim+2]*psi3Eta;
  
  const double detXY =  DxDxi*DyDeta - DxDeta*DyDxi;
  const double detYZ =  DyDxi*DzDeta - DyDeta*DzDxi;
  const double detXZ = -DxDxi*DzDeta + DxDeta*DzDxi;
  
  const double det = std::sqrt( detXY*detXY + detYZ*detYZ + detXZ*detXZ );
    
  normal[0] = detYZ / det;
  normal[1] = detXZ / det;
  normal[2] = detXY / det;
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Quad93DSCS::Quad93DSCS()
  : HexahedralP2Element(),
    surfaceDimension_(2)
{
  // set up the one-dimensional quadrature rule
  set_quadrature_rule();

  // set up integration rule and relevant maps on scs
  set_interior_info();

  // compute and save shape functions and derivatives at ips
  eval_shape_functions_at_ips();
  eval_shape_derivs_at_ips();
  eval_shape_functions_at_shifted_ips();
  eval_shape_derivs_at_shifted_ips();
}

//--------------------------------------------------------------------------
//-------- set_interior_info -----------------------------------------------
//--------------------------------------------------------------------------
void
Quad93DSCS::set_interior_info()
{
  nodesPerElement_ = nodes1D_ * nodes1D_;

  std::vector<int> nodeMap = {
                               0, 4, 1,   // bottom row of nodes
                               7, 8, 5,   // middle row of nodes
                               3, 6, 2    // top row of nodes
                             };

  auto tensor_map_2D = [=] (int i, int j) { return nodeMap[i+nodes1D_*j]; };

  //1D integration rule per sub-control volume
   numIntPoints_ = (nodes1D_ * nodes1D_) * ( numQuad_ * numQuad_ ); // 36

   // define ip node mappings
   ipNodeMap_.resize(numIntPoints_);
   intgLoc_.resize(numIntPoints_*surfaceDimension_); // size = 72
   intgLocShift_.resize(numIntPoints_*surfaceDimension_); // size = 72
   ipWeight_.resize(numIntPoints_);

   // tensor product nodes (3x3) x tensor product quadrature (2x2)
   int vector_index_2D = 0; int scalar_index = 0;
   for (int l = 0; l < nodes1D_; ++l) {
     for (int k = 0; k < nodes1D_; ++k) {
       for (int j = 0; j < numQuad_; ++j) {
         for (int i = 0; i < numQuad_; ++i) {
           //integration point location
           intgLoc_[vector_index_2D]     = gauss_point_location(k,i);
           intgLoc_[vector_index_2D + 1] = gauss_point_location(l,j);

           intgLocShift_[vector_index_2D]     = shifted_gauss_point_location(k,i);
           intgLocShift_[vector_index_2D + 1] = shifted_gauss_point_location(l,j);

           //weight
           ipWeight_[scalar_index] = tensor_product_weight(k,l,i,j);

           //sub-control volume association
           ipNodeMap_[scalar_index] = tensor_map_2D(k,l);

           // increment indices
           ++scalar_index;
           vector_index_2D += surfaceDimension_;
         }
       }
     }
   }
}

//--------------------------------------------------------------------------
//-------- eval_shape_derivs_at_ips ----------------------------------------
//--------------------------------------------------------------------------
void
Quad93DSCS::eval_shape_functions_at_ips()
{
  shapeFunctions_.resize(numIntPoints_*nodesPerElement_);
  quad9_shape_fcn(numIntPoints_, intgLoc_.data(), shapeFunctions_.data());
}

//--------------------------------------------------------------------------
//-------- eval_shape_derivs_at_ips ----------------------------------------
//--------------------------------------------------------------------------
void
Quad93DSCS::eval_shape_derivs_at_ips()
{
  shapeDerivs_.resize(numIntPoints_*nodesPerElement_*surfaceDimension_);
  quad9_shape_deriv(numIntPoints_, intgLoc_.data(), shapeDerivs_.data());
}

//--------------------------------------------------------------------------
//-------- eval_shape_derivs_at_ips ----------------------------------------
//--------------------------------------------------------------------------
void
Quad93DSCS::eval_shape_functions_at_shifted_ips()
{
  shapeFunctionsShift_.resize(numIntPoints_*nodesPerElement_);
  quad9_shape_fcn(numIntPoints_, intgLocShift_.data(), shapeFunctionsShift_.data());
}

//--------------------------------------------------------------------------
//-------- eval_shape_derivs_at_ips ----------------------------------------
//--------------------------------------------------------------------------
void
Quad93DSCS::eval_shape_derivs_at_shifted_ips()
{
  shapeDerivsShift_.resize(numIntPoints_*nodesPerElement_*surfaceDimension_);
  quad9_shape_deriv(numIntPoints_, intgLocShift_.data(), shapeDerivsShift_.data());
}

//--------------------------------------------------------------------------
//-------- quad9_shape_fcn -------------------------------------------------
//--------------------------------------------------------------------------
void
Quad93DSCS::quad9_shape_fcn(
  int  numIntPoints,
  const double *intgLoc,
  double *shpfc) const
{
  for ( int ip = 0; ip < numIntPoints; ++ip ) {
    int nineIp = 9*ip; // nodes per element is always 9
    int k = 2*ip;
    const double s = intgLoc[k];
    const double t = intgLoc[k+1];

    const double one_m_s = 1.0 - s;
    const double one_p_s = 1.0 + s;
    const double one_m_t = 1.0 - t;
    const double one_p_t = 1.0 + t;

    const double one_m_ss = 1.0 - s * s;
    const double one_m_tt = 1.0 - t * t;

    shpfc[nineIp  ] =  0.25 * s * t *  one_m_s *  one_m_t;
    shpfc[nineIp+1] = -0.25 * s * t *  one_p_s *  one_m_t;
    shpfc[nineIp+2] =  0.25 * s * t *  one_p_s *  one_p_t;
    shpfc[nineIp+3] = -0.25 * s * t *  one_m_s *  one_p_t;
    shpfc[nineIp+4] = -0.50 *     t *  one_p_s *  one_m_s * one_m_t;
    shpfc[nineIp+5] =  0.50 * s     *  one_p_t *  one_m_t * one_p_s;
    shpfc[nineIp+6] =  0.50 *     t *  one_p_s *  one_m_s * one_p_t;
    shpfc[nineIp+7] = -0.50 * s     *  one_p_t *  one_m_t * one_m_s;
    shpfc[nineIp+8] =  one_m_ss * one_m_tt;
  }
}


//--------------------------------------------------------------------------
//-------- quad9_shape_deriv -----------------------------------------------
//--------------------------------------------------------------------------
void
Quad93DSCS::quad9_shape_deriv(
  int numIntPoints,
  const double *intgLoc,
  double *deriv) const
{
  for ( int ip = 0; ip < numIntPoints; ++ip ) {
    const int grad_offset = surfaceDimension_ * nodesPerElement_ * ip; // nodes per element is always 9
    const int vector_offset = surfaceDimension_ * ip;
    int node; int offset;

    const double s = intgLoc[vector_offset+0];
    const double t = intgLoc[vector_offset+1];

    const double s2 = s*s;
    const double t2 = t*t;

    node = 0;
    offset = grad_offset + surfaceDimension_ * node;
    deriv[offset+0] = 0.25 * (2.0 * s * t2 - 2.0 * s * t - t2 + t);
    deriv[offset+1] = 0.25 * (2.0 * s2 * t - 2.0 * s * t - s2 + s);

    node = 1;
    offset = grad_offset + surfaceDimension_ * node;
    deriv[offset+0] = 0.25 * (2.0 * s * t2 - 2.0 * s * t + t2 - t);
    deriv[offset+1] = 0.25 * (2.0 * s2 * t + 2.0 * s * t - s2 - s);

    node = 2;
    offset = grad_offset + surfaceDimension_ * node;
    deriv[offset+0] = 0.25 * (2.0 * s * t2 + 2.0 * s * t + t2 + t);
    deriv[offset+1] = 0.25 * (2.0 * s2 * t + 2.0 * s * t + s2 + s);

    node = 3;
    offset = grad_offset + surfaceDimension_ * node;
    deriv[offset+0] = 0.25 * (2.0 * s * t2 + 2.0 * s * t - t2 - t);
    deriv[offset+1] = 0.25 * (2.0 * s2 * t - 2.0 * s * t + s2 - s);

    node = 4;
    offset = grad_offset + surfaceDimension_ * node;
    deriv[offset+0] = -0.5 * (2.0 * s * t2 - 2.0 * s * t);
    deriv[offset+1] = -0.5 * (2.0 * s2 * t - s2 - 2.0 * t + 1.0);

    node = 5;
    offset = grad_offset + surfaceDimension_ * node;
    deriv[offset+0] = -0.5 * (2.0 * s * t2 + t2 - 2.0 * s - 1.0);
    deriv[offset+1] = -0.5 * (2.0 * s2 * t + 2.0 * s * t);

    node = 6;
    offset = grad_offset + surfaceDimension_ * node;
    deriv[offset+0] = -0.5 * (2.0 * s * t2 + 2.0 * s * t);
    deriv[offset+1] = -0.5 * (2.0 * s2 * t + s2 - 2.0 * t - 1.0);

    node = 7;
    offset = grad_offset + surfaceDimension_ * node;
    deriv[offset+0] = -0.5 * (2.0 * s * t2 - t2 - 2.0 * s + 1.0);
    deriv[offset+1] = -0.5 * (2.0 * s2 * t - 2.0 * s * t);

    node = 8;
    offset = grad_offset + surfaceDimension_ * node;
    deriv[offset+0] = 2.0 * s * t2 - 2.0 * s;
    deriv[offset+1] = 2.0 * s2 * t - 2.0 * t;
  }
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad93DSCS::ipNodeMap(
  int /*ordinal*/)
{
  // define ip->node mappings for each face (single ordinal); 
  return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void
Quad93DSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  std::array<double,3> areaVector;

  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;
    const int vector_elem_offset = nDim_ * numIntPoints_ * k;

    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int grad_offset = surfaceDimension_ * nodesPerElement_ * ip;
      const int offset = nDim_ * ip + vector_elem_offset;

      //compute area vector for this ip
      area_vector( &coords[coord_elem_offset],
                   &shapeDerivs_[grad_offset],
                   areaVector.data() );

      // apply quadrature weight and orientation (combined as weight)
      for (int j = 0; j < nDim_; ++j) {
        areav[offset+j]  = ipWeight_[ip] * areaVector[j];
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- isInElement -----------------------------------------------------
//--------------------------------------------------------------------------
double
Quad93DSCS::isInElement(
  const double *elemNodalCoord,
  const double *pointCoord,
  double *isoParCoord )
{
  const double isInElemConverged = 1.0e-16;
  // Translate element so that (x,y,z) coordinates of the first node are (0,0,0)

  double x[3] = { elemNodalCoord[1] - elemNodalCoord[0],
                  elemNodalCoord[2] - elemNodalCoord[0],
                  elemNodalCoord[3] - elemNodalCoord[0] };
  
  double y[3] = { elemNodalCoord[10] - elemNodalCoord[9],
                  elemNodalCoord[11] - elemNodalCoord[9],
                  elemNodalCoord[12] - elemNodalCoord[9] };
  
  double z[3] = { elemNodalCoord[19] - elemNodalCoord[18],
                  elemNodalCoord[20] - elemNodalCoord[18],
                  elemNodalCoord[21] - elemNodalCoord[18] };
  
  // (xp,yp,zp) is the point at which we're searching for (xi,eta,d)
  // (must translate this also)
  // d = (scaled) distance in (x,y,z) space from point (xp,yp,zp) to the
  //     surface defined by the face element (the distance is scaled by
  //     the length of the non-unit normal vector; rescaling of d is done
  //     following the NR iteration below).

  double xp = pointCoord[0] - elemNodalCoord[0];
  double yp = pointCoord[1] - elemNodalCoord[9];
  double zp = pointCoord[2] - elemNodalCoord[18];

  // Newton-Raphson iteration for (xi,eta,d)

  double jdet;
  double j[9];
  double gn[3];
  double xcur[3];          // current (x,y,z) point on element surface
  double normal[3];        // (non-unit) normal computed at xcur

  // Solution vector solcur[3] = {xi,eta,d}
  double solcur[3] = {-0.5,-0.5,-0.5};     // initial guess
  double deltasol[] = {1.0,1.0, 1.0};

  unsigned i = 0;
  const unsigned MAX_NR_ITER = 100;

  do
  {
    // Update guess vector
    solcur[0] += deltasol[0];
    solcur[1] += deltasol[1];
    solcur[2] += deltasol[2];

    interpolatePoint(3,solcur,elemNodalCoord,xcur);

    // Translate xcur ((x,y,z) point corresponding
    // to current (xi,eta) guess)

    xcur[0] -= elemNodalCoord[0];
    xcur[1] -= elemNodalCoord[9];
    xcur[2] -= elemNodalCoord[18];

    non_unit_face_normal(solcur,elemNodalCoord,normal);

    gn[0] = xcur[0] - xp + solcur[2] * normal[0];
    gn[1] = xcur[1] - yp + solcur[2] * normal[1];
    gn[2] = xcur[2] - zp + solcur[2] * normal[2];

    // Mathematica-generated code for the jacobian

    j[0]=0.125000000000000*(-2.00000000000000*(-1.00000000000000+solcur[1])*x[0]+(2.00000000000000*(1.00000000000000+solcur[1])*(x[1]-x[2])+solcur[2]*(-(y[1]*z[0])+y[2]*z[0]+y[0]*z[1]-y[0]*z[2])));
    
    j[1]=0.125000000000000*(-2.00000000000000*(1.00000000000000+solcur[0])*x[0]+2.00000000000000*(1.00000000000000+solcur[0])*x[1]-2.00000000000000*(-1.00000000000000+solcur[0])*x[2]+(solcur[2]*(y[2]*(z[0]-z[1])+(-y[0]+y[1])*z[2])));
    
    j[2]= normal[0];

    j[3]=0.125000000000000*(-2.00000000000000*(-1.00000000000000+solcur[1])*y[0]+(2.00000000000000*(1.00000000000000+solcur[1])*(y[1]-y[2])+solcur[2]*(x[1]*z[0]-x[2]*z[0]-x[0]*z[1]+x[0]*z[2])));

    j[4]=0.125000000000000*(-2.00000000000000*(1.00000000000000+solcur[0])*y[0]+2.00000000000000*(1.00000000000000+solcur[0])*y[1]-2.00000000000000*(-1.00000000000000+solcur[0])*y[2]+(solcur[2]*(x[2]*(-z[0]+z[1])+(x[0]-x[1])*z[2])));

    j[5]= normal[1];

    j[6]=0.125000000000000*((solcur[2]*(-(x[1]*y[0])+x[2]*y[0]+x[0]*y[1]-x[0]*y[2]))-2.00000000000000*((-1.00000000000000+solcur[1])*z[0]-(1.00000000000000+solcur[1])*(z[1]-z[2])));

    j[7]=0.125000000000000*((solcur[2]*(x[2]*(y[0]-y[1])+(-x[0]+x[1])*y[2]))-2.00000000000000*(1.00000000000000+solcur[0])*z[0]+2.00000000000000*(1.00000000000000+solcur[0])*z[1]-2.00000000000000*(-1.00000000000000+solcur[0])*z[2]);

    j[8]= normal[2];

    jdet=-(j[2]*j[4]*j[6])+j[1]*j[5]*j[6]+j[2]*j[3]*j[7]-
	   j[0]*j[5]*j[7]-j[1]*j[3]*j[8]+j[0]*j[4]*j[8];

    // Solve linear system (j*deltasol = -gn) for deltasol at step n+1

    deltasol[0] = (gn[2]*(j[2]*j[4]-j[1]*j[5])+gn[1]*(-(j[2]*j[7])+
                                                      j[1]*j[8])+gn[0]*(j[5]*j[7]-j[4]*j[8]))/jdet;
    deltasol[1] = (gn[2]*(-(j[2]*j[3])+j[0]*j[5])+gn[1]*(j[2]*j[6]-
                                                         j[0]*j[8])+gn[0]*(-(j[5]*j[6])+j[3]*j[8]))/jdet;
    deltasol[2] = (gn[2]*(j[1]*j[3]-j[0]*j[4])+gn[1]*(-(j[1]*j[6])+
                                                      j[0]*j[7])+gn[0]*(j[4]*j[6]-j[3]*j[7]))/jdet;
    
  } while ( !within_tolerance( vector_norm_sq(deltasol,3), isInElemConverged) &&
	    ++i < MAX_NR_ITER );

  // Fill in solution vector; only include the distance (in the third
  // solution vector slot) if npar_coord = 3 (this is how the user
  // requests it)

  isoParCoord[0] = isoParCoord[1] = isoParCoord[2] = std::numeric_limits<double>::max();
  double dist = std::numeric_limits<double>::max();

  if ( i < MAX_NR_ITER ) {
    isoParCoord[0] = solcur[0] + deltasol[0];
    isoParCoord[1] = solcur[1] + deltasol[1];
    // Rescale the distance vector by the length of the (non-unit) normal vector,
    // which was used above in the NR iteration.
    const double area   = std::sqrt(vector_norm_sq(normal,3));
    const double length = std::sqrt(area);
    const double isoParCoord_2  = (solcur[2] + deltasol[2]) * length;
    isoParCoord[2] = isoParCoord_2;

    std::vector<double> xtmp = { isoParCoord[0], isoParCoord[1], isoParCoord_2 };
    dist = parametric_distance(xtmp);
  }
  return dist;
}

//--------------------------------------------------------------------------
//-------- interpolatePoint ------------------------------------------------
//--------------------------------------------------------------------------
void
Quad93DSCS::interpolatePoint(
  const int &nComp,
  const double *isoParCoord,
  const double *field,
  double *result )
{
  double one = 1.0;
  double half   = one / 2.0;
  double qtr    = one / 4.0;

  double s = isoParCoord[0];
  double t = isoParCoord[1];

  double one_m_s = one - s;
  double one_p_s = one + s;
  double one_m_t = one - t;
  double one_p_t = one + t;

  double one_m_ss = one - s * s;
  double one_m_tt = one - t * t;

  for ( int i = 0; i < nComp; i++ ) {
    int b = 9*i;       // Base 'field array' index for ith component

    result[i] =   qtr * s * t *  one_m_s * one_m_t * field[b+ 0]+
      -qtr * s * t *  one_p_s *  one_m_t * field[b+ 1]+
      qtr * s * t *  one_p_s *  one_p_t * field[b+ 2]+
      -qtr * s * t *  one_m_s *  one_p_t * field[b+ 3]+
      -half * t * one_p_s * one_m_s * one_m_t * field[b+ 4]+
      half * s * one_p_t * one_m_t * one_p_s * field[b+ 5]+
      half * t * one_p_s * one_m_s * one_p_t * field[b+ 6]+
      -half * s * one_p_t * one_m_t * one_m_s * field[b+ 7]+
      one_m_ss * one_m_tt * field[b+ 8];
  }
}

//--------------------------------------------------------------------------
//-------- general_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Quad93DSCS::general_shape_fcn(
  const int numIp,
  const double *isoParCoord,
  double *shpfc)
{
  quad9_shape_fcn(numIp, isoParCoord, shpfc);
}

//--------------------------------------------------------------------------
//-------- general_normal --------------------------------------------------
//--------------------------------------------------------------------------
void
Quad93DSCS::general_normal(
  const double *isoParCoord,
  const double *coords,
  double *normal)
{
  // coords(3,9)
  const int nDim = 3;

  const double s = isoParCoord[0];
  const double t = isoParCoord[1];
  
  const double t2 = t*t;
  const double s2 = s*s;
  
  const double psi0Xi  =  0.25 * (2.0 * s  * t2 - 2.0*s*t-t2+t);
  const double psi1Xi  =  0.25 * (2.0 * s  * t2 - 2.0*s*t+t2-t);
  const double psi2Xi  =  0.25 * (2.0 * s  * t2 + 2.0*s*t+t2+t);
  const double psi3Xi  =  0.25 * (2.0 * s  * t2 + 2.0*s*t-t2-t);
  const double psi4Xi  =  -0.5 * (2.0 * s  * t2 - 2.0*s*t);
  const double psi5Xi  =  -0.5 * (2.0 * s  * t2 + t2 - 2.0*s-1.0);
  const double psi6Xi  =  -0.5 * (2.0 * s  * t2 + 2.0*s*t);
  const double psi7Xi  =  -0.5 * (2.0 * s  * t2 - t2 - 2.0*s+1.0);
  const double psi8Xi  =          2.0 * s  * t2      - 2.0*s;
  
  const double psi0Eta = 0.25 * (2.0 * s2 * t  - 2.0*s*t-s2+s);
  const double psi1Eta = 0.25 * (2.0 * s2 * t  + 2.0*s*t-s2-s);
  const double psi2Eta = 0.25 * (2.0 * s2 * t  + 2.0*s*t+s2+s);
  const double psi3Eta = 0.25 * (2.0 * s2 * t  - 2.0*s*t+s2-s);
  const double psi4Eta = -0.5 * (2.0 * s2 * t  - s2 - 2.0*t+1.0);
  const double psi5Eta = -0.5 * (2.0 * s2 * t  + 2.0*s*t);
  const double psi6Eta = -0.5 * (2.0 * s2 * t  + s2 - 2.0*t-1.0);
  const double psi7Eta = -0.5 * (2.0 * s2 * t  - 2.0*s*t);
  const double psi8Eta =         2.0 * s2 * t       - 2.0*t;
  
  const double DxDxi = coords[0*nDim+0]*psi0Xi +
    coords[1*nDim+0]*psi1Xi +
    coords[2*nDim+0]*psi2Xi +
    coords[3*nDim+0]*psi3Xi +
    coords[4*nDim+0]*psi4Xi +
    coords[5*nDim+0]*psi5Xi +
    coords[6*nDim+0]*psi6Xi +
    coords[7*nDim+0]*psi7Xi +
    coords[8*nDim+0]*psi8Xi;
  
  const double DyDxi = coords[0*nDim+1]*psi0Xi +
    coords[1*nDim+1]*psi1Xi +
    coords[2*nDim+1]*psi2Xi +
    coords[3*nDim+1]*psi3Xi +
    coords[4*nDim+1]*psi4Xi +
    coords[5*nDim+1]*psi5Xi +
    coords[6*nDim+1]*psi6Xi +
    coords[7*nDim+1]*psi7Xi +
    coords[8*nDim+1]*psi8Xi;
  
  const double DzDxi = coords[0*nDim+2]*psi0Xi +
    coords[1*nDim+2]*psi1Xi +
    coords[2*nDim+2]*psi2Xi +
    coords[3*nDim+2]*psi3Xi +
    coords[4*nDim+2]*psi4Xi +
    coords[5*nDim+2]*psi5Xi +
    coords[6*nDim+2]*psi6Xi +
    coords[7*nDim+2]*psi7Xi +
    coords[8*nDim+2]*psi8Xi;
  
  const double DxDeta = coords[0*nDim+0]*psi0Eta +
    coords[1*nDim+0]*psi1Eta +
    coords[2*nDim+0]*psi2Eta +
    coords[3*nDim+0]*psi3Eta +
    coords[4*nDim+0]*psi4Eta +
    coords[5*nDim+0]*psi5Eta +
    coords[6*nDim+0]*psi6Eta +
    coords[7*nDim+0]*psi7Eta +
    coords[8*nDim+0]*psi8Eta;
  
  const double DyDeta = coords[0*nDim+1]*psi0Eta +
    coords[1*nDim+1]*psi1Eta +
    coords[2*nDim+1]*psi2Eta +
    coords[3*nDim+1]*psi3Eta +
    coords[4*nDim+1]*psi4Eta +
    coords[5*nDim+1]*psi5Eta +
    coords[6*nDim+1]*psi6Eta +
    coords[7*nDim+1]*psi7Eta +
    coords[8*nDim+1]*psi8Eta;

  const double DzDeta = coords[0*nDim+2]*psi0Eta +
    coords[1*nDim+2]*psi1Eta +
    coords[2*nDim+2]*psi2Eta +
    coords[3*nDim+2]*psi3Eta +
    coords[4*nDim+2]*psi4Eta +
    coords[5*nDim+2]*psi5Eta +
    coords[6*nDim+2]*psi6Eta +
    coords[7*nDim+2]*psi7Eta +
    coords[8*nDim+2]*psi8Eta;
  
  const double detXY =  DxDxi*DyDeta - DxDeta*DyDxi;
  const double detYZ =  DyDxi*DzDeta - DyDeta*DzDxi;
  const double detXZ = -DxDxi*DzDeta + DxDeta*DzDxi;
  
  const double det = std::sqrt( detXY*detXY + detYZ*detYZ + detXZ*detXZ );
  
  normal[0] = detYZ / det;
  normal[1] = detXZ / det;
  normal[2] = detXY / det;
}

//--------------------------------------------------------------------------
//-------- non_unit_face_normal --------------------------------------------
//--------------------------------------------------------------------------
void
Quad93DSCS::non_unit_face_normal(	
  const double * isoParCoord, 
  const double * elemNodalCoord,
  double * normalVector )
{
  double xi  = isoParCoord[0];
  double eta = isoParCoord[1];

  // Translate element so that node 0 is at (x,y,z) = (0,0,0)

  double x[3] = { elemNodalCoord[1] - elemNodalCoord[0],
                  elemNodalCoord[2] - elemNodalCoord[0],
                  elemNodalCoord[3] - elemNodalCoord[0] };
  
  double y[3] = { elemNodalCoord[10] - elemNodalCoord[9],
                  elemNodalCoord[11] - elemNodalCoord[9],
                  elemNodalCoord[12] - elemNodalCoord[9] };
  
  double z[3] = { elemNodalCoord[19] - elemNodalCoord[18],
                  elemNodalCoord[20] - elemNodalCoord[18],
                  elemNodalCoord[21] - elemNodalCoord[18] };
  
  // Mathematica-generated and simplified code for the normal vector

  const double n0 = 0.125000000000000*(xi*y[2]*z[0]+y[0]*z[1]+xi*y[0]*z[1]-y[2]*z[1]-
                                       xi*y[0]*z[2]+y[1]*(-((1.00000000000000+xi)*z[0])+
                                                          (1.00000000000000+eta)*z[2])+eta*(y[2]*z[0]-y[2]*z[1]-y[0]*z[2]));

  const double n1 = 0.125000000000000*(-(xi*x[2]*z[0])-x[0]*z[1]-xi*x[0]*z[1]+x[2]*z[1]+
                                       xi*x[0]*z[2]+x[1]*((1.00000000000000+xi)*z[0]-
                                                          (1.00000000000000+eta)*z[2])+eta*(-(x[2]*z[0])+x[2]*z[1]+x[0]*z[2]));
  
  const double n2 = 0.125000000000000*(xi*x[2]*y[0]+x[0]*y[1]+xi*x[0]*y[1]-x[2]*y[1]-
                                       xi*x[0]*y[2]+x[1]*(-((1.00000000000000+xi)*y[0])+
                                                          (1.00000000000000+eta)*y[2])+eta*(x[2]*y[0]-x[2]*y[1]-x[0]*y[2]));
  
  normalVector[0] = n0;
  normalVector[1] = n1;
  normalVector[2] = n2;
}

//--------------------------------------------------------------------------
//-------- parametric_distance ---------------------------------------------
//--------------------------------------------------------------------------
double Quad93DSCS::parametric_distance(const std::vector<double> &x)
{
  const double ELEM_THICK  = 0.01;
  std::vector<double> y = { std::fabs(x[0]), std::fabs(x[1]), std::fabs(x[2]) };
  double d = y[0];
  if (d < y[1]) d = y[1];
  if (ELEM_THICK < y[2] && d < 1+y[2]) d = 1+y[2];
  return d;
}

//--------------------------------------------------------------------------
//-------- area_vector -----------------------------------------------------
//--------------------------------------------------------------------------
void
Quad93DSCS::area_vector(
  const double *elemNodalCoords,
  const double *shapeDeriv,
  double *areaVector) const
{
   // return the normal area vector given shape derivatives dnds OR dndt
   double dx_ds1 = 0.0; double dy_ds1 = 0.0; double dz_ds1 = 0.0;
   double dx_ds2 = 0.0; double dy_ds2 = 0.0; double dz_ds2 = 0.0;

   for (int node = 0; node < nodesPerElement_; ++node) {
     const int vector_offset = nDim_ * node;
     const int surface_vector_offset = surfaceDimension_ * node;

     const double xCoord = elemNodalCoords[vector_offset+0];
     const double yCoord = elemNodalCoords[vector_offset+1];
     const double zCoord = elemNodalCoords[vector_offset+2];

     const double dn_ds1 = shapeDeriv[surface_vector_offset+0];
     const double dn_ds2 = shapeDeriv[surface_vector_offset+1];

     dx_ds1 += dn_ds1 * xCoord;
     dx_ds2 += dn_ds2 * xCoord;

     dy_ds1 += dn_ds1 * yCoord;
     dy_ds2 += dn_ds2 * yCoord;

     dz_ds1 += dn_ds1 * zCoord;
     dz_ds2 += dn_ds2 * zCoord;
   }

   //cross product
   areaVector[0] = dy_ds1 * dz_ds2 - dz_ds1 * dy_ds2;
   areaVector[1] = dz_ds1 * dx_ds2 - dx_ds1 * dz_ds2;
   areaVector[2] = dx_ds1 * dy_ds2 - dy_ds1 * dx_ds2;
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Tri3DSCS::Tri3DSCS()
  : MasterElement()
{
  nDim_ = 3;
  nodesPerElement_ = 3;
  numIntPoints_ = 3;

  // define ip node mappings; ordinal size = 1
  ipNodeMap_.resize(3);
  ipNodeMap_[0] = 0;
  ipNodeMap_[1] = 1;
  ipNodeMap_[2] = 2;

  // standard integration location
  intgLoc_.resize(6);
  const double five24th = 5.0/24.0;
  const double seven12th = 7.0/12.0;
  intgLoc_[0]  = five24th;  intgLoc_[1] = five24th;  // surf 1
  intgLoc_[2]  = seven12th; intgLoc_[3] = five24th;  // surf 2
  intgLoc_[4]  = five24th;  intgLoc_[5] = seven12th; // surf 3

  // shifted
  intgLocShift_.resize(6);
  intgLocShift_[0]  =  0.00; intgLocShift_[1] =  0.00; // surf 1
  intgLocShift_[2]  =  1.00; intgLocShift_[3] =  0.00; // surf 2
  intgLocShift_[4]  =  0.00; intgLocShift_[5] =  1.00; // surf 3
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Tri3DSCS::~Tri3DSCS()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Tri3DSCS::ipNodeMap(
  int /*ordinal*/)
{
  // define ip->node mappings for each face (single ordinal); 
  return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void Tri3DSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(tri3d_scs_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_,
      coords, areav );

  // fake check
  *error = (lerr == 0) ? 0.0 : 1.0;
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Tri3DSCS::shape_fcn(double *shpfc)
{
  tri_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Tri3DSCS::shifted_shape_fcn(double *shpfc)
{
  tri_shape_fcn(numIntPoints_, &intgLocShift_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- tri_shape_fcn ---------------------------------------------------
//--------------------------------------------------------------------------
void
Tri3DSCS::tri_shape_fcn(
  const int  &npts,
  const double *isoParCoord,
  double *shape_fcn)
{
  for (int j = 0; j < npts; ++j ) {
    const int threej = 3*j;
    const int k = 2*j;
    const double xi = isoParCoord[k];
    const double eta = isoParCoord[k+1];
    shape_fcn[    threej] = 1.0 - xi - eta;
    shape_fcn[1 + threej] = xi;
    shape_fcn[2 + threej] = eta;
  }
}

//--------------------------------------------------------------------------
//-------- isInElement -----------------------------------------------------
//--------------------------------------------------------------------------
double
Tri3DSCS::isInElement(
    const double * elem_nodal_coor,
    const double * point_coor,
	  double * par_coor ) 
{
  // always intended for 3D...
  const int npar_coord = 3;
  // Translate element so that (x,y,z) coordinates of the
  // first node
  double x[2] = { elem_nodal_coor[1] - elem_nodal_coor[0],
                  elem_nodal_coor[2] - elem_nodal_coor[0] };
  double y[2] = { elem_nodal_coor[4] - elem_nodal_coor[3],
                  elem_nodal_coor[5] - elem_nodal_coor[3] };
  double z[2] = { elem_nodal_coor[7] - elem_nodal_coor[6],
                  elem_nodal_coor[8] - elem_nodal_coor[6] };

  // Translate position vector of point in same manner

  double xp = point_coor[0] - elem_nodal_coor[0];
  double yp = point_coor[1] - elem_nodal_coor[3];
  double zp = point_coor[2] - elem_nodal_coor[6];

  // Set new nodal coordinates with Node 1 at origin and with new
  // x and y axes lying in the plane of the element
  double len12 = std::sqrt( x[0]*x[0] + y[0]*y[0] + z[0] *z[0] );
  double len13 = std::sqrt( x[1]*x[1] + y[1]*y[1] + z[1] *z[1] );

  double xnew[3];
  double ynew[3];
  double znew[3];

  // Use cross-product of 12 and 13 to find enclosed angle and
  // direction of new z-axis

  znew[0] = y[0]*z[1] - y[1]*z[0];
  znew[1] = x[1]*z[0] - x[0]*z[1];
  znew[2] = x[0]*y[1] - x[1]*y[0];

  double Area2 = std::sqrt( znew[0]*znew[0] + znew[1]*znew[1] +
                            znew[2]*znew[2] );

  // find sin of angle
  double sin_theta = Area2 / ( len12 * len13 ) ;

  // find cosine of angle
  double cos_theta = (x[0]*x[1] + y[0]*y[1] + z[0]*z[1])/(len12 * len13);

  // nodal coordinates of nodes 2 and 3 in new system
  // (coordinates of node 1 are identically 0.0)
  double x_nod_new[2] = { len12, len13*cos_theta};
  double y_nod_new[2] = {  0.0, len13*sin_theta};

  // find direction cosines transform position vector of
  // point to be checked into new coordinate system

  // direction cosines of new x axis along side 12

  xnew[0] = x[0]/len12;
  xnew[1] = y[0]/len12;
  xnew[2] = z[0]/len12;

  // direction cosines of new z axis
  znew[0] = znew[0]/Area2;
  znew[1] = znew[1]/Area2;
  znew[2] = znew[2]/Area2;

  // direction cosines of new y-axis (cross-product of znew and xnew)
  ynew[0] = znew[1]*xnew[2] - xnew[1]*znew[2];
  ynew[1] = xnew[0]*znew[2] - znew[0]*xnew[2];
  ynew[2] = znew[0]*xnew[1] - xnew[0]*znew[1];

  // compute transformed coordinates of point
  // (coordinates in xnew,ynew,znew)
  double xpnew = xnew[0]*xp + xnew[1]*yp + xnew[2]*zp;
  double ypnew = ynew[0]*xp + ynew[1]*yp + ynew[2]*zp;
  double zpnew = znew[0]*xp + znew[1]*yp + znew[2]*zp;

  // Find parametric coordinates of point and check that
  // it lies in the element
  par_coor[0] = 1. - xpnew / x_nod_new[0] +
		 ypnew*( x_nod_new[1] - x_nod_new[0] ) / Area2;
  par_coor[1] = ( xpnew*y_nod_new[1] - ypnew*x_nod_new[1] ) / Area2;

  if (3 == npar_coord) par_coor[2] = zpnew/std::sqrt(Area2);

  std::vector<double> w = { par_coor[0], par_coor[1], zpnew/std::sqrt(Area2) };

  par_coor[0] = w[1];
  par_coor[1] = 1.0-w[0]-w[1];

  const double dist = parametric_distance(w);

  return dist;
}

//--------------------------------------------------------------------------
//-------- parametric_distance ---------------------------------------------
//--------------------------------------------------------------------------
double 
Tri3DSCS::parametric_distance(
  const std::vector<double> &x)
{
  const double ELEM_THICK = 0.01;
  const double X=x[0] - 1./3.;
  const double Y=x[1] - 1./3.;
  const double dist0 = -3*X;
  const double dist1 = -3*Y;
  const double dist2 =  3*(X+Y);
  double dist = std::max(std::max(dist0,dist1),dist2);
  const double y = std::fabs(x[2]);
  if (ELEM_THICK < y && dist < 1+y) dist = 1+y;
  return dist;
}

//--------------------------------------------------------------------------
//-------- interpolatePoint ------------------------------------------------
//--------------------------------------------------------------------------
void
Tri3DSCS::interpolatePoint(
  const int  & ncomp_field,
  const double * isoParCoord,
  const double * field,
  double * result )
{
  const double r = isoParCoord[0];
  const double s = isoParCoord[1];
  const double t = 1.0 - r - s;

  for ( int i = 0; i < ncomp_field; i++ ) {
    int b = 3*i;  //Base 'field array' index for ith component
    result[i] = t*field[b] + r*field[b+1] + s*field[b+2];
  }
}

//--------------------------------------------------------------------------
//-------- general_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Tri3DSCS::general_shape_fcn(
  const int numIp,
  const double *isoParCoord,
  double *shpfc)
{
  tri_shape_fcn(numIp, isoParCoord, shpfc);
}


//--------------------------------------------------------------------------
//-------- general_normal --------------------------------------------------
//--------------------------------------------------------------------------
void
Tri3DSCS::general_normal(
  const double */*isoParCoord*/,
  const double *coords,
  double *normal)
{
  // can be only linear
  const double ax  = coords[3] - coords[0];
  const double ay  = coords[4] - coords[1];
  const double az  = coords[5] - coords[2];
  const double bx  = coords[6] - coords[0];
  const double by  = coords[7] - coords[1];
  const double bz  = coords[8] - coords[2];

  normal[0] = ( ay*bz - az*by );
  normal[1] = ( az*bx - ax*bz );
  normal[2] = ( ax*by - ay*bx );

  const double mag = std::sqrt( normal[0]*normal[0] +
                                normal[1]*normal[1] +
                                normal[2]*normal[2] );
  normal[0] /= mag;
  normal[1] /= mag;
  normal[2] /= mag;
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Edge2DSCS::Edge2DSCS()
  : MasterElement(),
    elemThickness_(0.01)
{
  nDim_ = 2;
  nodesPerElement_ = 2;
  numIntPoints_ = 2;
  scaleToStandardIsoFac_ = 2.0;

  // define ip node mappings; ordinal size = 1
  ipNodeMap_.resize(2);
  ipNodeMap_[0] = 0;
  ipNodeMap_[1] = 1;

  intgLoc_.resize(2);
  intgLoc_[0]  =  -0.25; intgLoc_[1]  = 0.25;
 
  intgLocShift_.resize(2);
  intgLocShift_[0]  =  -0.50; intgLocShift_[1]  = 0.50; 
  
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Edge2DSCS::~Edge2DSCS()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Edge2DSCS::ipNodeMap(
  int /*ordinal*/)
{
  // define ip->node mappings for each face (single ordinal); 
  return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void Edge2DSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  int lerr = 0;

  SIERRA_FORTRAN(edge2d_scs_det)
    ( &nelem, &nodesPerElement_, &numIntPoints_,
      coords, areav );

  // fake check
  *error = (lerr == 0) ? 0.0 : 1.0;
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Edge2DSCS::shape_fcn(double *shpfc)
{
  for ( int i =0; i < nodesPerElement_; ++i ) {
    int j = 2*i;
    shpfc[j  ] = 0.5-intgLoc_[i];
    shpfc[j+1] = 0.5+intgLoc_[i];
  }
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Edge2DSCS::shifted_shape_fcn(double *shpfc)
{
  for ( int i =0; i< nodesPerElement_; ++i ) {
    int j = 2*i;
    shpfc[j  ] = 0.5-intgLocShift_[i];
    shpfc[j+1] = 0.5+intgLocShift_[i];
  }
}

//--------------------------------------------------------------------------
//-------- isInElement -----------------------------------------------------
//--------------------------------------------------------------------------
double
Edge2DSCS::isInElement(
    const double * elem_nodal_coor,     // (2,2)
    const double * point_coor,          // (2)
	  double * par_coor ) 
{
  // elem_nodal_coor has the endpoints of the line
  // segment defining this element.  Set the first
  // endpoint to zero.  This means subtrace the
  // first endpoint from the second.
  const double X1 = elem_nodal_coor[1]-elem_nodal_coor[0];
  const double X2 = elem_nodal_coor[3]-elem_nodal_coor[2];

  // Now subtract the first endpoint from the target point
  const double P1 = point_coor[0] - elem_nodal_coor[0];
  const double P2 = point_coor[1] - elem_nodal_coor[2];

  // Now find the projection along the line of the point
  // This is the parametric coordinate in range (0,1)
  const double norm2 = X1*X1 + X2*X2;
  
  const double xi = (P1*X1 + P2*X2) / norm2;
  // rescale to (-1,1)
  par_coor[0] = 2*xi - 1;

  // Now find the projection from the point to a perpenducular
  // line.  This gives the distance from the point to the element.
  const double alpha = std::abs(P1*X2 - P2*X1) / norm2;
  if (2 == nDim_) 
    par_coor[1] = alpha;

  std::vector<double> x(2);
  x[0] = par_coor[0];
  x[1] = alpha;
  const double dist = parametric_distance(x);

  return dist;
}

//--------------------------------------------------------------------------
//-------- parametric_distance ---------------------------------------------
//--------------------------------------------------------------------------
double
Edge2DSCS::parametric_distance(const std::vector<double> &x)
{
  double dist = std::fabs(x[0]);
  if (elemThickness_ < x[1] && dist < 1.0+x[1]) 
    dist = 1+x[1];
  return dist;
}

//--------------------------------------------------------------------------
//-------- interpolatePoint ------------------------------------------------
//--------------------------------------------------------------------------
void
Edge2DSCS::interpolatePoint(
  const int &nComp,
  const double *isoParCoord,
  const double *field,
  double *result )
{
  double xi = isoParCoord[0]; 
  for ( int i = 0; i < nComp; i++ ) {
    // Base 'field array' index for ith component
    int b = 2*i;
    result[i] = 0.5*(1.0-xi) * field[b+0] +
      0.5*(1.0+xi) * field[b+1];
  }
}

//--------------------------------------------------------------------------
//-------- general_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Edge2DSCS::general_shape_fcn(
  const int numIp,
  const double *isoParCoord,
  double *shpfc)
{
  const double npe = nodesPerElement_;
  for ( int ip = 0; ip < numIp; ++ip ) {
    int j = npe*ip;
    shpfc[j  ] = 0.5*(1.0-isoParCoord[ip]);
    shpfc[j+1] = 0.5*(1.0+isoParCoord[ip]);
  }
}

//--------------------------------------------------------------------------
//-------- general_normal --------------------------------------------------
//--------------------------------------------------------------------------
void
Edge2DSCS::general_normal(
  const double */*isoParCoord*/,
  const double *coords,
  double *normal)
{
  // can be only linear
  const double dx  = coords[2] - coords[0];
  const double dy  = coords[3] - coords[1];
  const double mag = std::sqrt(dx*dx + dy*dy);

  normal[0] =  dy/mag;
  normal[1] = -dx/mag;
}

//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Edge32DSCS::Edge32DSCS()
  : QuadrilateralP2Element()
{
  nodesPerElement_ = nodes1D_;

  // set up the one-dimensional quadrature rule
  set_quadrature_rule();

  numIntPoints_ = numQuad_ * nodes1D_;

  ipNodeMap_.resize(numIntPoints_);
  intgLoc_.resize(numIntPoints_);

  intgLocShift_.resize(6);
  intgLocShift_[0] = -1.00; intgLocShift_[1]  = -1.00;
  intgLocShift_[2] =  0.00; intgLocShift_[3]  =  0.00;
  intgLocShift_[4] =  1.00; intgLocShift_[5]  =  1.00;

  ipWeight_.resize(numIntPoints_);

  std::vector<int> stk1DNodeMap = { 0, 2, 1 };

  int scalar_index = 0;
  for (int k = 0; k < nodes1D_; ++k) {
    for (int i = 0; i < numQuad_; ++i) {
      ipNodeMap_[scalar_index] = stk1DNodeMap[k];
      intgLoc_[scalar_index] = gauss_point_location(k,i);
      ipWeight_[scalar_index] = tensor_product_weight(k,i);

      ++scalar_index;
    }
  }
}

//--------------------------------------------------------------------------
//-------- ipNodeMap -------------------------------------------------------
//--------------------------------------------------------------------------
const int *
Edge32DSCS::ipNodeMap(
  int /*ordinal*/)
{
  // define ip->node mappings for each face (single ordinal); 
  return &ipNodeMap_[0];
}

//--------------------------------------------------------------------------
//-------- determinant -----------------------------------------------------
//--------------------------------------------------------------------------
void Edge32DSCS::determinant(
  const int nelem,
  const double *coords,
  double *areav,
  double *error)
{
  std::array<double,2> areaVector;

  for (int k = 0; k < nelem; ++k) {
    const int coord_elem_offset = nDim_ * nodesPerElement_ * k;

    for (int ip = 0; ip < numIntPoints_; ++ip) {
      const int offset = nDim_ * ip + coord_elem_offset;

      // calculate the area vector
      area_vector( &coords[coord_elem_offset],
                   intgLoc_[ip],
                   areaVector.data() );

      // weight the area vector with the Gauss-quadrature weight for this IP
      areav[offset + 0] = ipWeight_[ip] * areaVector[0];
      areav[offset + 1] = ipWeight_[ip] * areaVector[1];
    }
  }

  // check
  *error = 0.0;
}

//--------------------------------------------------------------------------
//-------- shape_fcn -------------------------------------------------------
//--------------------------------------------------------------------------
void
Edge32DSCS::shape_fcn(double *shpfc)
{
  for ( int i =0; i< numIntPoints_; ++i ) {
    int j = 3*i;
    const double s = intgLoc_[i];
    shpfc[j  ] = -s*(1.0-s)*0.5;
    shpfc[j+1] = s*(1.0+s)*0.5;
    shpfc[j+2] = (1.0-s)*(1.0+s);
  }
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
Edge32DSCS::shifted_shape_fcn(double *shpfc)
{
  for ( int i =0; i< numIntPoints_; ++i ) {
    int j = 3*i;
    const double s = intgLocShift_[i];
    shpfc[j  ] = -s*(1.0-s)*0.5;
    shpfc[j+1] = s*(1.0+s)*0.5;
    shpfc[j+2] = (1.0-s)*(1.0+s);
  }
}

//--------------------------------------------------------------------------
//-------- area_vector -----------------------------------------------------
//--------------------------------------------------------------------------
void
Edge32DSCS::area_vector(
  const double *coords,
  const double s,
  double *areaVector) const
{
  // returns the normal area vector (dyds,-dxds) evaluated at s

  // create a parameterization of the curve
  // r(s) = (x(s),y(s)) s.t. r(-1) = (x0,y0); r(0) = (x2,y2); r(1) = (x1,y1);
  // x(s) = x2 + 0.5 (x1-x0) s + 0.5 (x1 - 2 x2 + x0) s^2,
  // y(s) = y2 + 0.5 (y1-y0) s + 0.5 (y1 - 2 y2 + y0) s^2
  // could equivalently use the shape function derivatives . . .

  // coordinate names
  const double x0 = coords[0]; const double y0 = coords[1];
  const double x1 = coords[2]; const double y1 = coords[3];
  const double x2 = coords[4]; const double y2 = coords[5];

  const double dxds = 0.5 * (x1 - x0) + (x1 - 2.0 * x2 + x0) * s;
  const double dyds = 0.5 * (y1 - y0) + (y1 - 2.0 * y2 + y0) * s;

  areaVector[0] =  dyds;
  areaVector[1] = -dxds;
}

} // namespace nalu
} // namespace sierra
