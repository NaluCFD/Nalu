/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <master_element/MasterElement.h>
#include <master_element/MasterElementFunctions.h>

#include <master_element/MasterElementHO.h>
#include <master_element/MasterElementUtils.h>

#include <element_promotion/LagrangeBasis.h>
#include <element_promotion/TensorProductQuadratureRule.h>
#include <element_promotion/QuadratureRule.h>
#include <AlgTraits.h>

#include <NaluEnv.h>
#include <FORTRAN_Proto.h>

#include <stk_util/environment/ReportHandler.hpp>
#include <stk_topology/topology.hpp>

#include <iostream>

#include <cmath>
#include <limits>
#include <array>
#include <map>
#include <memory>

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
//-------- constructor------------------------------------------------------
//--------------------------------------------------------------------------
QuadrilateralP2Element::QuadrilateralP2Element()
  : MasterElement(),
    scsDist_(std::sqrt(3.0)/3.0),
    nodes1D_(3),
    numQuad_(2)
{
  nDim_ = 2;
  nodesPerElement_ = nodes1D_ * nodes1D_;

  // map the standard stk (refinement consistent) node numbering
  // to a tensor-product style node numbering (i.e. node (m,l,k) -> m+npe*l+npe^2*k)
  stkNodeMap_ = {
                  0, 4, 1, // bottom row of nodes
                  7, 8, 5, // middle row of nodes
                  3, 6, 2  // top row of nodes
                };

  sideNodeOrdinals_ = {
      0, 1, 4,
      1, 2, 5,
      2, 3, 6,
      3, 0, 7
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
  std::tie(gaussAbscissae_, gaussWeight_) = gauss_legendre_rule(numQuad_);
  for (unsigned j = 0; j < gaussWeight_.size(); ++j) {
    gaussWeight_[j] *= 0.5;
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
   return isoparametric_mapping( scsEndLoc_[nodeOrdinal+1],
     scsEndLoc_[nodeOrdinal],
     gaussAbscissae_[gaussPointOrdinal] );
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
  const double Ls1 = scsEndLoc_[s1Node+1]-scsEndLoc_[s1Node];
  const double Ls2 = scsEndLoc_[s2Node+1]-scsEndLoc_[s2Node];
  const double isoparametricArea = Ls1 * Ls2;

  const double weight = isoparametricArea * gaussWeight_[s1Ip] * gaussWeight_[s2Ip];

  return weight;
}

//--------------------------------------------------------------------------
//-------- tensor_product_weight -------------------------------------------
//--------------------------------------------------------------------------
double
QuadrilateralP2Element::tensor_product_weight(int s1Node, int s1Ip) const
{
  //line integration
  const double isoparametricLength = scsEndLoc_[s1Node+1]-scsEndLoc_[s1Node];
  const double weight = isoparametricLength * gaussWeight_[s1Ip];

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
//-------- parametric_distance ---------------------------------------------
//--------------------------------------------------------------------------
double QuadrilateralP2Element::parametric_distance(const std::array<double, 2>& x)
{
  double absXi  = std::abs(x[0]);
  double absEta = std::abs(x[1]);
  return (absXi > absEta) ? absXi : absEta;
}

//--------------------------------------------------------------------------
//-------- interpolatePoint ------------------------------------------------
//--------------------------------------------------------------------------
void
QuadrilateralP2Element::interpolatePoint(
  const int &nComp,
  const double *isoParCoord,
  const double *field,
  double *result )
{
  constexpr int nNodes = 9;
  std::array<double, nNodes> shapefct;
  quad9_shape_fcn(1, isoParCoord, shapefct.data());

  for (int i = 0; i < nComp; i++) {
    result[i] = ddot(shapefct.data(), field + nNodes * i, nNodes);
  }
}

//--------------------------------------------------------------------------
//-------- isInElement -----------------------------------------------------
//--------------------------------------------------------------------------
double QuadrilateralP2Element::isInElement(
  const double *elemNodalCoord,
  const double *pointCoord,
  double *isoParCoord)
{
  // control the interation
  double isInElemConverged = 1.0e-16; // NOTE: the square of the tolerance on the distance
  int N_MAX_ITER = 100;

  constexpr int dim = 2;
  std::array<double, dim> guess = { { 0.0, 0.0 } };
  std::array<double, dim> delta;
  int iter = 0;

  do {
    // interpolate coordinate at guess
    constexpr int nNodes = 9;
    std::array<double, nNodes> weights;
    quad9_shape_fcn(1, guess.data(), weights.data());

    // compute difference between coordinates interpolated to the guessed isoParametric coordinates
    // and the actual point's coordinates
    std::array<double, dim> error_vec;
    error_vec[0] = pointCoord[0] - ddot(weights.data(), elemNodalCoord + 0 * nNodes, nNodes);
    error_vec[1] = pointCoord[1] - ddot(weights.data(), elemNodalCoord + 1 * nNodes, nNodes);

    // update guess along gradient of mapping from physical-to-reference coordinates
    // transpose of the jacobian of the forward mapping
    constexpr int deriv_size = nNodes * dim;
    std::array<double, deriv_size> deriv;
    quad9_shape_deriv(1, guess.data(), deriv.data());

    std::array<double, dim * dim> jact{};
    for(int j = 0; j < nNodes; ++j) {
      jact[0] += deriv[0 + j * dim] * elemNodalCoord[j + 0 * nNodes];
      jact[1] += deriv[1 + j * dim] * elemNodalCoord[j + 0 * nNodes];
      jact[2] += deriv[0 + j * dim] * elemNodalCoord[j + 1 * nNodes];
      jact[3] += deriv[1 + j * dim] * elemNodalCoord[j + 1 * nNodes];
    }

    // apply its inverse on the error vector
    solve22(jact.data(), error_vec.data(), delta.data());

    // update guess
    guess[0] += delta[0];
    guess[1] += delta[1];

    //continue to iterate if update was larger than the set tolerance until max iterations are reached
  } while(!within_tolerance(vector_norm_sq(delta.data(), 2), isInElemConverged) && (++iter < N_MAX_ITER));

  // output if failed:
  isoParCoord[0] = std::numeric_limits<double>::max();
  isoParCoord[1] = std::numeric_limits<double>::max();
  double dist = std::numeric_limits<double>::max();

  if (iter < N_MAX_ITER) {
    // output if succeeded:
    isoParCoord[0] = guess[0];
    isoParCoord[1] = guess[1];
    dist = parametric_distance(guess);
  }
  return dist;
}

void
QuadrilateralP2Element::sidePcoords_to_elemPcoords(
  const int & side_ordinal,
  const int & npoints,
  const double *side_pcoords,
  double *elem_pcoords)
{
  switch (side_ordinal) {
  case 0:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*2+0] = side_pcoords[i];
      elem_pcoords[i*2+1] = -1;
    }
    break;
  case 1:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*2+0] = 1;
      elem_pcoords[i*2+1] = side_pcoords[i];
    }
    break;
  case 2:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*2+0] = -side_pcoords[i];
      elem_pcoords[i*2+1] = 1;
    }
    break;
  case 3:
    for (int i=0; i<npoints; i++) {
      elem_pcoords[i*2+0] = -1;
      elem_pcoords[i*2+1] = -side_pcoords[i];
    }
    break;
  default:
    throw std::runtime_error("QuadrilateralP2Element::sideMap invalid ordinal");
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
    for (int ip = 0; ip < Traits::numScvIp_; ++ip) {
      const int grad_offset = nDim_ * nodesPerElement_ * ip;

      //weighted jacobian determinant
      const double det_j = jacobian_determinant(coords, &shapeDerivs_[grad_offset]);

      //apply weight and store to volume
      volume[ip] = ipWeight_[ip] * det_j;

      //flag error
      if (det_j < tiny_positive_value()) {
        *error = 1.0;
      }
    }

}

//--------------------------------------------------------------------------
//-------- jacobian_determinant --------------------------------------------
//--------------------------------------------------------------------------
double
Quad92DSCV::jacobian_determinant(
  const double *POINTER_RESTRICT elemNodalCoords,
  const double *POINTER_RESTRICT shapeDerivs) const
{
  double dx_ds1 = 0.0;  double dx_ds2 = 0.0;
  double dy_ds1 = 0.0;  double dy_ds2 = 0.0;

  for (int node = 0; node < Traits::nodesPerElement_; ++node) {
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
//-------- side_node_ordinals ----------------------------------------------
//--------------------------------------------------------------------------
const int *
Quad92DSCS::side_node_ordinals(
  int ordinal)
{
  // define face_ordinal->node_ordinal mappings for each face (ordinal);
  return &sideNodeOrdinals_[ordinal*3];
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

  ThrowRequireMsg(nelem == 1, "P2 elements are processed one-at-a-time");

  constexpr int dim = Traits::nDim_;
  constexpr int ipsPerDirection = Traits::numScsIp_ / dim;
  static_assert ( ipsPerDirection * dim == Traits::numScsIp_, "Number of ips incorrect");

  constexpr int deriv_increment = dim * Traits::nodesPerElement_;

  int index = 0;

   //returns the normal vector x_u x x_s for constant t surfaces
  for (int ip = 0; ip < ipsPerDirection; ++ip) {
    ThrowAssert(ipInfo_[index].direction == Jacobian::T_DIRECTION);
    area_vector<Jacobian::T_DIRECTION>(coords, &shapeDerivs_[deriv_increment * index], &areav[index*dim]);
    ++index;
  }

  //returns the normal vector x_t x x_u for constant s curves
  for (int ip = 0; ip < ipsPerDirection; ++ip) {
    ThrowAssert(ipInfo_[index].direction == Jacobian::S_DIRECTION);
    area_vector<Jacobian::S_DIRECTION>(coords, &shapeDerivs_[deriv_increment * index], &areav[index*dim]);
    ++index;
  }

  // Multiply with the integration point weighting
  for (int ip = 0; ip < Traits::numScsIp_; ++ip) {
    double weight = ipInfo_[ip].weight;
    areav[ip * dim + 0] *= weight;
    areav[ip * dim + 1] *= weight;
  }

  *error = 0; // no error checking available
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

  constexpr int numShapeDerivs = Traits::numScsIp_*Traits::nodesPerElement_*Traits::nDim_;
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

  constexpr int numShapeDerivs = Traits::numScsIp_*Traits::nodesPerElement_*Traits::nDim_;
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
  ThrowRequireMsg(nelem == 1, "P2 elements are processed one-at-a-time");

  int lerr = 0;

  const int nface = 1;
  const int face_offset =  nDim_ * ipsPerFace_ * nodesPerElement_ * face_ordinal;
  double* offsetFaceDerivs = &expFaceShapeDerivs_[face_offset];

  for (int ip = 0; ip < ipsPerFace_; ++ip) {
    const int grad_offset = nDim_ * nodesPerElement_ * ip;

    SIERRA_FORTRAN(quad_gradient_operator)
    ( & nface,
        &nodesPerElement_,
        &nface,
        &offsetFaceDerivs[grad_offset],
        coords,
        &gradop[grad_offset],
        &det_j[ip],
        error,
        &lerr
    );

    if (det_j[ip] < tiny_positive_value() || lerr != 0) {
      *error = 1.0;
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
template <Jacobian::Direction direction> void
Quad92DSCS::area_vector(
  const double *POINTER_RESTRICT elemNodalCoords,
  double *POINTER_RESTRICT shapeDeriv,
  double *POINTER_RESTRICT normalVec ) const
{
  constexpr int s1Component = (direction == Jacobian::S_DIRECTION) ?
      Jacobian::T_DIRECTION : Jacobian::S_DIRECTION;

  double dxdr = 0.0;  double dydr = 0.0;
  for (int node = 0; node < Traits::nodesPerElement_; ++node) {
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


  sideNodeOrdinals_ = {
      0, 1,  // ordinal 0
      1, 2,  // ordinal 1
      2, 0   // ordinal 2
  };

  std::vector<std::vector<double>> nodeLocations =
  {
      {0.0,0.0}, {1.0,0}, {0.0,1.0}
  };
  intgExpFaceShift_.resize(12);
  int index = 0;
  stk::topology topo = stk::topology::TRIANGLE_3_2D;
  for (unsigned k = 0; k < topo.num_sides(); ++k) {
    stk::topology side_topo = topo.side_topology(k);
    const int* ordinals = side_node_ordinals(k);
    for (unsigned n = 0; n < side_topo.num_nodes(); ++n) {
      intgExpFaceShift_[2*index + 0] = nodeLocations[ordinals[n]][0];
      intgExpFaceShift_[2*index + 1] = nodeLocations[ordinals[n]][1];
      ++index;
    }
  }
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
//-------- side_node_ordinals ----------------------------------------------
//--------------------------------------------------------------------------
const int *
Tri2DSCS::side_node_ordinals(
  int ordinal)
{
  // define face_ordinal->node_ordinal mappings for each face (ordinal);
  return &sideNodeOrdinals_[ordinal*2];
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
//-------- shifted_face_grad_op --------------------------------------------
//--------------------------------------------------------------------------
void Tri2DSCS::shifted_face_grad_op(
  const int nelem,
  const int /*face_ordinal*/,
  const double *coords,
  double *gradop,
  double *det_j,
  double *error)
{
  // same as regular face_grad_op

  int lerr = 0;
  int npf = 2;

  const int nface = 1;
  double dpsi[6];

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
          &coords[12*n], &gradop[k*nelem*6+n*6], &det_j[npf*n+k], error, &lerr );

      if ( lerr )
        std::cout << "sorry, issue with face_grad_op.." << std::endl;

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
  
  // derivatives are constant
  SIERRA_FORTRAN(tri_derivative)
    ( &nface, dpsi );
      
  SIERRA_FORTRAN(tri_gradient_operator)
    ( &nface,
      &nodesPerElement_,
      &nface,
      dpsi,
      &coords[0], &gradop[0], &det_j[0], error, &lerr );
      
  if ( lerr )
    std::cout << "sorry, issue with face_grad_op.." << std::endl;
  
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
//-------- interpolate_point -----------------------------------------------
//--------------------------------------------------------------------------
void
Edge32DSCS::interpolatePoint(
  const int &nComp,
  const double *isoParCoord,
  const double *field,
  double *result)
{
  constexpr int nNodes = 3;

  double s = isoParCoord[0];
  std::array<double, nNodes> shapefct = {{-0.5*s*(1-s), +0.5*s*(1+s), (1-s)*(1+s)}};

  for ( int i =0; i< nComp; ++i ) {
    result[i] = shapefct[0] * field[3*i+0] + shapefct[1] * field[3*i+1] + shapefct[2] * field[3*i+2];
  }
}

//--------------------------------------------------------------------------
//-------- area_vector -----------------------------------------------------
//--------------------------------------------------------------------------
void
Edge32DSCS::area_vector(
  const double *POINTER_RESTRICT coords,
  const double s,
  double *POINTER_RESTRICT areaVector) const
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
