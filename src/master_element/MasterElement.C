/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <master_element/MasterElement.h>
#include <FORTRAN_Proto.h>

#include <stk_topology/topology.hpp>

#include <iostream>

#include <cmath>
#include <limits>

namespace sierra{
namespace nalu{

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
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
HexSCV::HexSCV()
  : MasterElement()
{
  nDim_ = 3;
  nodesPerElement_ = 8;
  numIntPoints_ = 8;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
HexSCV::~HexSCV()
{
  // does nothing
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

  // mapping between exposed face and extruded element's overlapping face   
  faceNodeOnExtrudedElem_.resize(24);
  faceNodeOnExtrudedElem_[0]  = 1; faceNodeOnExtrudedElem_[1]  = 0; faceNodeOnExtrudedElem_[2]  = 4; faceNodeOnExtrudedElem_[3]  = 5;
  faceNodeOnExtrudedElem_[4]  = 1; faceNodeOnExtrudedElem_[5]  = 0; faceNodeOnExtrudedElem_[6]  = 4; faceNodeOnExtrudedElem_[7]  = 5;
  faceNodeOnExtrudedElem_[8]  = 1; faceNodeOnExtrudedElem_[9]  = 0; faceNodeOnExtrudedElem_[10] = 4; faceNodeOnExtrudedElem_[11] = 5;
  faceNodeOnExtrudedElem_[12] = 0; faceNodeOnExtrudedElem_[13] = 4; faceNodeOnExtrudedElem_[14] = 5; faceNodeOnExtrudedElem_[15] = 1;
  faceNodeOnExtrudedElem_[16] = 0; faceNodeOnExtrudedElem_[17] = 4; faceNodeOnExtrudedElem_[18] = 5; faceNodeOnExtrudedElem_[19] = 1;
  faceNodeOnExtrudedElem_[20] = 1; faceNodeOnExtrudedElem_[21] = 0; faceNodeOnExtrudedElem_[22] = 4; faceNodeOnExtrudedElem_[23] = 5;

  // mapping between exposed face and extruded element's opposing face
  opposingNodeOnExtrudedElem_.resize(24);
  opposingNodeOnExtrudedElem_[0]  = 2; opposingNodeOnExtrudedElem_[1]  = 3; opposingNodeOnExtrudedElem_[2]  = 7; opposingNodeOnExtrudedElem_[3]  = 6;
  opposingNodeOnExtrudedElem_[4]  = 2; opposingNodeOnExtrudedElem_[5]  = 3; opposingNodeOnExtrudedElem_[6]  = 7; opposingNodeOnExtrudedElem_[7]  = 6;
  opposingNodeOnExtrudedElem_[8]  = 2; opposingNodeOnExtrudedElem_[9]  = 3; opposingNodeOnExtrudedElem_[10] = 7; opposingNodeOnExtrudedElem_[11] = 6;
  opposingNodeOnExtrudedElem_[12] = 3; opposingNodeOnExtrudedElem_[13] = 7; opposingNodeOnExtrudedElem_[14] = 6; opposingNodeOnExtrudedElem_[15] = 2;
  opposingNodeOnExtrudedElem_[16] = 3; opposingNodeOnExtrudedElem_[17] = 7; opposingNodeOnExtrudedElem_[18] = 6; opposingNodeOnExtrudedElem_[19] = 2;
  opposingNodeOnExtrudedElem_[20] = 2; opposingNodeOnExtrudedElem_[21] = 3; opposingNodeOnExtrudedElem_[22] = 7; opposingNodeOnExtrudedElem_[23] = 6;

  // mapping between exposed face scs ips and halo edge
  faceScsIpOnExtrudedElem_.resize(24);
  faceScsIpOnExtrudedElem_[0]  = 1; faceScsIpOnExtrudedElem_[1]  = 3; faceScsIpOnExtrudedElem_[2]  = 7; faceScsIpOnExtrudedElem_[3]  = 5;
  faceScsIpOnExtrudedElem_[4]  = 1; faceScsIpOnExtrudedElem_[5]  = 3; faceScsIpOnExtrudedElem_[6]  = 7; faceScsIpOnExtrudedElem_[7]  = 5;
  faceScsIpOnExtrudedElem_[8]  = 1; faceScsIpOnExtrudedElem_[9]  = 3; faceScsIpOnExtrudedElem_[10] = 7; faceScsIpOnExtrudedElem_[11] = 5;
  faceScsIpOnExtrudedElem_[12] = 3; faceScsIpOnExtrudedElem_[13] = 7; faceScsIpOnExtrudedElem_[14] = 5; faceScsIpOnExtrudedElem_[15] = 1;
  faceScsIpOnExtrudedElem_[16] = 3; faceScsIpOnExtrudedElem_[17] = 7; faceScsIpOnExtrudedElem_[18] = 5; faceScsIpOnExtrudedElem_[19] = 1;
  faceScsIpOnExtrudedElem_[20] = 1; faceScsIpOnExtrudedElem_[21] = 3; faceScsIpOnExtrudedElem_[22] = 7; faceScsIpOnExtrudedElem_[23] = 5;

  // mapping between exposed face scs ips and exposed face edge
  faceScsIpOnFaceEdges_.resize(24);
  faceScsIpOnFaceEdges_[0]  = 0; faceScsIpOnFaceEdges_[1]  = 8; faceScsIpOnFaceEdges_[2]  = 4; faceScsIpOnFaceEdges_[3]  = 9;
  faceScsIpOnFaceEdges_[4]  = 0; faceScsIpOnFaceEdges_[5]  = 8; faceScsIpOnFaceEdges_[6]  = 4; faceScsIpOnFaceEdges_[7]  = 9;
  faceScsIpOnFaceEdges_[8]  = 0; faceScsIpOnFaceEdges_[9]  = 8; faceScsIpOnFaceEdges_[10] = 4; faceScsIpOnFaceEdges_[11] = 9;
  faceScsIpOnFaceEdges_[12] = 8; faceScsIpOnFaceEdges_[13] = 4; faceScsIpOnFaceEdges_[14] = 9; faceScsIpOnFaceEdges_[15] = 0;
  faceScsIpOnFaceEdges_[16] = 8; faceScsIpOnFaceEdges_[17] = 4; faceScsIpOnFaceEdges_[18] = 9; faceScsIpOnFaceEdges_[19] = 0;
  faceScsIpOnFaceEdges_[20] = 0; faceScsIpOnFaceEdges_[21] = 8; faceScsIpOnFaceEdges_[22] = 4; faceScsIpOnFaceEdges_[23] = 9;
  
  // alignment of face:edge ordering and scsip area vector
  edgeAlignedArea_.resize(24);
  edgeAlignedArea_[0]  = -1.0; edgeAlignedArea_[1]  = +1.0; edgeAlignedArea_[2]  = +1.0; edgeAlignedArea_[3]  = -1.0;
  edgeAlignedArea_[4]  = -1.0; edgeAlignedArea_[5]  = +1.0; edgeAlignedArea_[6]  = +1.0; edgeAlignedArea_[7]  = -1.0;
  edgeAlignedArea_[8]  = -1.0; edgeAlignedArea_[9]  = +1.0; edgeAlignedArea_[10] = +1.0; edgeAlignedArea_[11] = -1.0;
  edgeAlignedArea_[12] = +1.0; edgeAlignedArea_[13] = +1.0; edgeAlignedArea_[14] = -1.0; edgeAlignedArea_[15] = -1.0;
  edgeAlignedArea_[16] = +1.0; edgeAlignedArea_[17] = +1.0; edgeAlignedArea_[18] = -1.0; edgeAlignedArea_[19] = -1.0;
  edgeAlignedArea_[20] = -1.0; edgeAlignedArea_[21] = +1.0; edgeAlignedArea_[22] = +1.0; edgeAlignedArea_[23] = -1.0;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
HexSCS::~HexSCS()
{
  // does nothing
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
  while ( !within_tol( vector_norm(xidiff,3), isInElemConverged) && ++i < maxNonlinearIter);

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
  const int ndim = 3;
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
//-------- faceNodeOnExtrudedElem ------------------------------------------
//--------------------------------------------------------------------------
const int *
HexSCS::faceNodeOnExtrudedElem()
{
  return &faceNodeOnExtrudedElem_[0];
}

//--------------------------------------------------------------------------
//-------- opposingNodeOnExtrudedElem --------------------------------------
//--------------------------------------------------------------------------
const int *
HexSCS::opposingNodeOnExtrudedElem()
{
  return &opposingNodeOnExtrudedElem_[0];
}

//--------------------------------------------------------------------------
//-------- faceScsIpOnExtrudedElem -----------------------------------------
//--------------------------------------------------------------------------
const int *
HexSCS::faceScsIpOnExtrudedElem()
{
  return &faceScsIpOnExtrudedElem_[0];
}

//--------------------------------------------------------------------------
//-------- faceScsIpOnFaceEdges --------------------------------------------
//--------------------------------------------------------------------------
const int *
HexSCS::faceScsIpOnFaceEdges()
{
  return &faceScsIpOnFaceEdges_[0];
}

//--------------------------------------------------------------------------
//-------- edgeAlignedArea -------------------------------------------------
//--------------------------------------------------------------------------
const double *
HexSCS::edgeAlignedArea()
{
  return &edgeAlignedArea_[0];
}

//--------------------------------------------------------------------------
//-------- vector_norm -----------------------------------------------------
//--------------------------------------------------------------------------
double HexSCS::vector_norm( const double * vect, int len )
{
  double norm_sq = 0.0;
  for (int i=0; i<len; i++)
  {
    norm_sq += vect[i]*vect[i];
  }
  return norm_sq;
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
//-------- within_tol ------------------------------------------------------
//--------------------------------------------------------------------------
bool HexSCS::within_tol( const double & val, const double & tol )
{
  return (std::abs(val)<tol);
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
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
TetSCV::~TetSCV()
{
  // does nothing
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

}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
TetSCS::~TetSCS()
{
  // does nothing
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
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
PyrSCV::PyrSCV()
  : MasterElement()
{
  nDim_ = 3;
  nodesPerElement_ = 5;
  numIntPoints_ = 5; 
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
PyrSCV::~PyrSCV()
{
  // does nothing
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
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
WedSCV::~WedSCV()
{
  // does nothing
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
  intgExpFace_[0] = 0.25;       intgExpFace_[1]  = 0.0;       intgExpFace_[2] = -0.5;  // surf 0,  nodes 0,1,4,3
  intgExpFace_[3] = 0.75;       intgExpFace_[4]  = 0.0;       intgExpFace_[5] = -0.5;  // face 0, surf 1
  intgExpFace_[6] = 0.75;       intgExpFace_[7]  = 0.0;       intgExpFace_[8] =  0.5;  // face 0, surf 2
  intgExpFace_[9] = 0.25;       intgExpFace_[10] = 0.0;       intgExpFace_[11] = 0.5;  // face 0, surf 3
  intgExpFace_[12] = 0.75;      intgExpFace_[13] = 0.25;      intgExpFace_[14] = -0.5; // surf 0,  nodes 1,2,5,4
  intgExpFace_[15] = 0.25;      intgExpFace_[16] = 0.75;      intgExpFace_[17] = -0.5; // face 1, surf 1
  intgExpFace_[18] = 0.25;      intgExpFace_[19] = 0.75;      intgExpFace_[20] =  0.5; // face 1, surf 2
  intgExpFace_[21] = 0.75;      intgExpFace_[22] = 0.25;      intgExpFace_[23] =  0.5; // face 1, surf 3
  intgExpFace_[24] = 0.0;       intgExpFace_[25] = 0.25;      intgExpFace_[26] = -0.5; // surf 0   nodes 0,3,5,2
  intgExpFace_[27] = 0.0;       intgExpFace_[28] = 0.25;      intgExpFace_[29] =  0.5; // face 2, surf 1
  intgExpFace_[30] = 0.0;       intgExpFace_[31] = 0.75;      intgExpFace_[32] =  0.5; // face 2, surf 2
  intgExpFace_[33] = 0.0;       intgExpFace_[34] = 0.75;      intgExpFace_[35] = -0.5; // face 2, surf 3
  intgExpFace_[36] = five24th;  intgExpFace_[37] = five24th;  intgExpFace_[38] = -1.0; // surf 0   nodes 0,2,1
  intgExpFace_[39] = five24th;  intgExpFace_[40] = seven12th; intgExpFace_[41] = -1.0; // face 3, surf 1
  intgExpFace_[42] = seven12th; intgExpFace_[43] = five24th;  intgExpFace_[44] = -1.0; // face 3, surf 2
  intgExpFace_[45] = 0.0;       intgExpFace_[46] = 0.0;       intgExpFace_[47] =  0.0; // (blank)
  intgExpFace_[48] = five24th;  intgExpFace_[49] = five24th;  intgExpFace_[50] = 1.0;  // surf 0   nodes 3,4,5
  intgExpFace_[51] = seven12th; intgExpFace_[52] = five24th;  intgExpFace_[53] = 1.0;  // face 4, surf 1
  intgExpFace_[54] = five24th;  intgExpFace_[55] = seven12th; intgExpFace_[56] = 1.0;  // face 4, surf 2
  intgExpFace_[57] = 0.0;       intgExpFace_[58] = 0.0;       intgExpFace_[59] = 0.0;  // (blank)
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
WedSCS::~WedSCS()
{
  // does nothing
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
void WedSCS::face_grad_op(
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
  wed_shape_fcn(numIntPoints_, &intgLoc_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- shifted_shape_fcn -----------------------------------------------
//--------------------------------------------------------------------------
void
WedSCS::shifted_shape_fcn(double *shpfc)
{
  wed_shape_fcn(numIntPoints_, &intgLocShift_[0], shpfc);
}

//--------------------------------------------------------------------------
//-------- wed_shape_fcn ---------------------------------------------------
//--------------------------------------------------------------------------
void
WedSCS::wed_shape_fcn(
  const int  &npts,
  const double *par_coord, 
  double *shape_fcn)
{
  for (int j = 0; j < npts; ++j ) {
    int sixj = 6 * j;
    int k    = 3 * j;
    double r    = par_coord[k];
    double s    = par_coord[k + 1];
    double t    = 1.0 - r - s;
    double xi   = par_coord[k + 2];
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
Quad2DSCV::Quad2DSCV()
  : MasterElement()
{
  nDim_ = 2;
  nodesPerElement_ = 4;
  numIntPoints_ = 4;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Quad2DSCV::~Quad2DSCV()
{
  // does nothing
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

  // mapping between exposed face and extruded element's overlapping face   
  faceNodeOnExtrudedElem_.resize(8);
  faceNodeOnExtrudedElem_[0] = 1; faceNodeOnExtrudedElem_[1] = 0;
  faceNodeOnExtrudedElem_[2] = 1; faceNodeOnExtrudedElem_[3] = 0;
  faceNodeOnExtrudedElem_[4] = 1; faceNodeOnExtrudedElem_[5] = 0;
  faceNodeOnExtrudedElem_[6] = 1; faceNodeOnExtrudedElem_[7] = 0;

  // mapping between exposed face and extruded element's opposing face
  opposingNodeOnExtrudedElem_.resize(8);
  opposingNodeOnExtrudedElem_[0] = 2; opposingNodeOnExtrudedElem_[1] = 3;
  opposingNodeOnExtrudedElem_[2] = 2; opposingNodeOnExtrudedElem_[3] = 3;
  opposingNodeOnExtrudedElem_[4] = 2; opposingNodeOnExtrudedElem_[5] = 3;
  opposingNodeOnExtrudedElem_[6] = 2; opposingNodeOnExtrudedElem_[7] = 3;

  // mapping between exposed face scs ips and halo edge
  faceScsIpOnExtrudedElem_.resize(8);
  faceScsIpOnExtrudedElem_[0] = 1; faceScsIpOnExtrudedElem_[1] = 3;
  faceScsIpOnExtrudedElem_[2] = 1; faceScsIpOnExtrudedElem_[3] = 3;
  faceScsIpOnExtrudedElem_[4] = 1; faceScsIpOnExtrudedElem_[5] = 3;
  faceScsIpOnExtrudedElem_[6] = 1; faceScsIpOnExtrudedElem_[7] = 3;

  // mapping between exposed face scs ips and exposed face edge
  faceScsIpOnFaceEdges_.resize(4);
  faceScsIpOnFaceEdges_[0] = 0; faceScsIpOnFaceEdges_[1] = 0;
  faceScsIpOnFaceEdges_[2] = 0; faceScsIpOnFaceEdges_[3] = 0;
  
  // alignment of face:edge ordering and scsip area vector
  edgeAlignedArea_.resize(4);
  edgeAlignedArea_[0] = -1.0; edgeAlignedArea_[1] = -1.0;
  edgeAlignedArea_[2] = -1.0; edgeAlignedArea_[3] = -1.0;

}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Quad2DSCS::~Quad2DSCS()
{
  // does nothing
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
  const double *par_coord, 
  double *shape_fcn)
{
  for (int j = 0; j < npts; ++j ) {
    const int fourj = 4*j;
    const int k = 2*j;
    const double s1 = par_coord[k];
    const double s2 = par_coord[k+1];
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
  const int ndim = 2;
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
//-------- faceNodeOnExtrudedElem ------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad2DSCS::faceNodeOnExtrudedElem()
{
  return &faceNodeOnExtrudedElem_[0];
}

//--------------------------------------------------------------------------
//-------- opposingNodeOnExtrudedElem --------------------------------------
//--------------------------------------------------------------------------
const int *
Quad2DSCS::opposingNodeOnExtrudedElem()
{
  return &opposingNodeOnExtrudedElem_[0];
}

//--------------------------------------------------------------------------
//-------- faceScsIpOnExtrudedElem -----------------------------------------
//--------------------------------------------------------------------------
const int *
Quad2DSCS::faceScsIpOnExtrudedElem()
{
  return &faceScsIpOnExtrudedElem_[0];
}

//--------------------------------------------------------------------------
//-------- faceScsIpOnFaceEdges --------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad2DSCS::faceScsIpOnFaceEdges()
{
  return &faceScsIpOnFaceEdges_[0];
}

//--------------------------------------------------------------------------
//-------- edgeAlignedArea -------------------------------------------------
//--------------------------------------------------------------------------
const double *
Quad2DSCS::edgeAlignedArea()
{
  return &edgeAlignedArea_[0];
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
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Tri2DSCV::~Tri2DSCV()
{
  // does nothing
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
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Tri2DSCS::~Tri2DSCS()
{
  // does nothing
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
  const double *par_coord, 
  double *shape_fcn)
{
  for (int j = 0; j < npts; ++j ) {
    const int threej = 3*j;
    const int k = 2*j;
    const double xi = par_coord[k];
    const double eta = par_coord[k+1];
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

  } while ( !within_tol( vector_norm2(deltasol,3), isInElemConverged)
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
    const double area   = std::sqrt(vector_norm2(normal,3));
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

bool 
Quad3DSCS::within_tol( const double & val, const double & tol )
{
  return (std::abs(val)<tol);
}

double 
Quad3DSCS::vector_norm2( const double * vect, int len )
{
  double norm_sq = 0.0;
  for (int i=0; i<len; i++)
  {
    norm_sq += vect[i]*vect[i];
  }
  return norm_sq;
}

void
Quad3DSCS::non_unit_face_normal(
  const double * par_coord,              // (2)
  const double * elem_nodal_coor,        // (4,3)
	double * normal_vector )         // (3)
{
  double xi  = par_coord[0];
  double eta = par_coord[1];

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
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Tri3DSCS::Tri3DSCS()
  : MasterElement()
{
  nDim_ = 3;
  nodesPerElement_ = 3;
  numIntPoints_ = 3;

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
  const double *par_coord,
  double *shape_fcn)
{
  for (int j = 0; j < npts; ++j ) {
    const int threej = 3*j;
    const int k = 2*j;
    const double xi = par_coord[k];
    const double eta = par_coord[k+1];
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
  const double * par_coord,
  const double * field,
  double * result )
{
  const double r = par_coord[0];
  const double s = par_coord[1];
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
  tri_shape_fcn(numIp, &isoParCoord[0], shpfc);
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

} // namespace nalu
} // namespace sierra
