/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MasterElement_h
#define MasterElement_h

#include <vector>
#include <cstdlib>
#include <stdexcept>

namespace stk {
struct topology;
}

namespace sierra{
namespace nalu{

namespace Jacobian{
enum Direction
{
  S_DIRECTION = 0,
  T_DIRECTION = 1,
  U_DIRECTION = 2
};
}

class MasterElement
{
public:

  MasterElement();
  virtual ~MasterElement();

  virtual void determinant(
    const int nelem,
    const double *coords,
    double *volume,
    double * error ) {
    throw std::runtime_error("determinant not implemented");}

  virtual void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error ) {
    throw std::runtime_error("grad_op not implemented");}

  virtual void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error ) {
    throw std::runtime_error("grad_op not implemented");}

  virtual void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv) {
    throw std::runtime_error("gij not implemented");}

  virtual void nodal_grad_op(
    const int nelem,
    double *deriv,
    double * error ) {
    throw std::runtime_error("nodal_grad_op not implemented");}

  virtual void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error ) {
    throw std::runtime_error("face_grad_op not implemented; avoid this element type at open bcs, walls and symms");}
  
  virtual const int * adjacentNodes() {
    throw std::runtime_error("adjacentNodes not implementedunknown bc");
    return NULL;}

  virtual const int * ipNodeMap(int ordinal = 0) {
      throw std::runtime_error("ipNodeMap not implemented");
      return NULL;}

  virtual void shape_fcn(
    double *shpfc) {
    throw std::runtime_error("shape_fcn not implemented"); }

  virtual void shifted_shape_fcn(
    double *shpfc) {
    throw std::runtime_error("shifted_shape_fcn not implemented"); }

  virtual int opposingNodes(
    const int ordinal, const int node) {
    throw std::runtime_error("adjacentNodes not implemented"); }

  virtual int opposingFace(
    const int ordinal, const int node) {
    throw std::runtime_error("opposingFace not implemented"); 
    return 0; }

  virtual double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord) {
    throw std::runtime_error("isInElement not implemented"); 
    return 1.0e6; }

  virtual void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result) {
    throw std::runtime_error("interpolatePoint not implemented"); }
  
  virtual void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc) {
    throw std::runtime_error("general_shape_fcn not implement"); }

  virtual void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error ) {
    throw std::runtime_error("general_face_grad_op not implemented");}

  virtual void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords) {
    throw std::runtime_error("sidePcoords_to_elemPcoords");}

  virtual const int * faceNodeOnExtrudedElem() {
    throw std::runtime_error("faceNodeOnExtrudedElem not implement"); }

  virtual const int * opposingNodeOnExtrudedElem() {
    throw std::runtime_error("opposingNodeOnExtrudedElem not implement"); }

  virtual const int * faceScsIpOnExtrudedElem() {
    throw std::runtime_error("faceScsIpOnExtrudedElem not implement"); }

  virtual const int * faceScsIpOnFaceEdges() {
    throw std::runtime_error("faceScsIpOnFaceEdges not implement"); }

  virtual const double * edgeAlignedArea() {
    throw std::runtime_error("edgeAlignedArea not implement"); }

  double isoparametric_mapping(const double b, const double a, const double xi) const;

  int nDim_;
  int nodesPerElement_;
  int numIntPoints_;
  double scaleToStandardIsoFac_;

  std::vector<int> lrscv_;
  std::vector<int> ipNodeMap_;
  std::vector<int> oppNode_;
  std::vector<int> oppFace_;
  std::vector<double> intgLoc_;
  std::vector<double> intgLocShift_;
  std::vector<double> intgExpFace_;
  std::vector<double> nodeLoc_;
  // extrusion-based scheme
  std::vector<int> faceNodeOnExtrudedElem_;
  std::vector<int> opposingNodeOnExtrudedElem_;
  std::vector<int> faceScsIpOnExtrudedElem_;
  std::vector<int> faceScsIpOnFaceEdges_;
  std::vector<double> edgeAlignedArea_;
  
};

// Hex 8 subcontrol volume
class HexSCV : public MasterElement
{
public:

  HexSCV();
  virtual ~HexSCV();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *volume,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shape_fcn(
    double *shpfc);
};

// Hex 8 subcontrol surface
class HexSCS : public MasterElement
{
public:

  HexSCS();
  virtual ~HexSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv);

  const int * adjacentNodes();

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);

  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);

  // extrusion data structure access
  const int * faceNodeOnExtrudedElem();
  const int * opposingNodeOnExtrudedElem();
  const int * faceScsIpOnExtrudedElem();
  const int * faceScsIpOnFaceEdges();
  const double * edgeAlignedArea();
  
  // helper
  double vector_norm( const double * vect, int len );
  double parametric_distance(const std::vector<double> &x);
  bool within_tol( const double & val, const double & tol );
};

class HexahedralP2Element : public MasterElement
{
public:
  HexahedralP2Element();
  virtual ~HexahedralP2Element() {}

  void shape_fcn(double *shpfc);
  void shifted_shape_fcn(double *shpfc);

protected:
  struct ContourData {
    Jacobian::Direction direction;
    double weight;
  };

  int tensor_product_node_map(int i, int j, int k) const;

  double gauss_point_location(
    int nodeOrdinal,
    int gaussPointOrdinal) const;

  double shifted_gauss_point_location(
    int nodeOrdinal,
    int gaussPointOrdinal) const;

  double tensor_product_weight(
    int s1Node, int s2Node, int s3Node,
    int s1Ip, int s2Ip, int s3Ip) const;

  double tensor_product_weight(
    int s1Node, int s2Node,
    int s1Ip, int s2Ip) const;

  virtual void eval_shape_functions_at_ips();
  virtual void eval_shape_functions_at_shifted_ips();

  virtual void eval_shape_derivs_at_ips();
  virtual void eval_shape_derivs_at_shifted_ips();

  void eval_shape_derivs_at_face_ips();

  void set_quadrature_rule();
  void GLLGLL_quadrature_weights();

  const double scsDist_;
  const bool useGLLGLL_;
  const int nodes1D_;
  const int numQuad_;

  // quadrature info
  std::vector<double> gaussAbscissae1D_;
  std::vector<double> gaussAbscissae_;
  std::vector<double> gaussAbscissaeShift_;
  std::vector<double> gaussWeight_;
  std::vector<double> scsEndLoc_;

  std::vector<int> stkNodeMap_;

  std::vector<double> shapeFunctions_;
  std::vector<double> shapeFunctionsShift_;
  std::vector<double> shapeDerivs_;
  std::vector<double> shapeDerivsShift_;
  std::vector<double> expFaceShapeDerivs_;

private:
  void hex27_shape_fcn(
    int npts,
    const double *par_coord,
    double* shape_fcn
  ) const;

  void hex27_shape_deriv(
    int npts,
    const double *par_coord,
    double* shape_fcn
  ) const;
};

// 3D Quad 27 subcontrol volume
class Hex27SCV : public HexahedralP2Element
{
public:
  Hex27SCV();
  virtual ~Hex27SCV() {}

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

private:
  void set_interior_info();

  double jacobian_determinant(
    const double *elemNodalCoords,
    const double *shapeDerivs ) const;

  std::vector<double> ipWeight_;
};

// 3D Hex 27 subcontrol surface
class Hex27SCS : public HexahedralP2Element
{
public:
  Hex27SCS();
  virtual ~Hex27SCS() {}

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv);

  const int * adjacentNodes();

  const int * ipNodeMap(int ordinal = 0);

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

private:
  void set_interior_info();
  void set_boundary_info();

  void area_vector(
    const Jacobian::Direction direction,
    const double *elemNodalCoords,
    double *shapeDeriv,
    double *areaVector ) const;

  void gradient(
    const double* elemNodalCoords,
    const double* shapeDeriv,
    double* grad,
    double* det_j ) const;

  std::vector<ContourData> ipInfo_;
  int ipsPerFace_;
};

// Tet 4 subcontrol volume
class TetSCV : public MasterElement
{
public:

  TetSCV();
  virtual ~TetSCV();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );
};

// Tet 4 subcontrol surface
class TetSCS : public MasterElement
{
public:

  TetSCS();
  virtual ~TetSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv);

  const int * adjacentNodes();

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void tet_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);

  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );
  
  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);

  // helper
  double parametric_distance(const std::vector<double> &x);
};

// Pyramid 5 subcontrol volume
class PyrSCV : public MasterElement
{
public:

  PyrSCV();
  virtual ~PyrSCV();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );
};

// Pyramid 5 subcontrol surface
class PyrSCS : public MasterElement
{
public:

  PyrSCS();
  virtual ~PyrSCS();

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void pyr_derivative(
    const int npts,
    const double *intLoc,
    double *deriv);

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv);

  const int * adjacentNodes();

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void pyr_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  int opposingNodes(
    const int ordinal, const int node);
};

// Wedge 6 subcontrol volume
class WedSCV : public MasterElement
{
public:
  WedSCV();
  virtual ~WedSCV();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );
};

// Wedge 6 subcontrol surface
class WedSCS : public MasterElement
{
public:
  WedSCS();
  virtual ~WedSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void wedge_derivative(
    const int npts,
    const double *intLoc,
    double *deriv);

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv);

  const int * adjacentNodes();

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);
  
  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);

  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);
  
  void wedge_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  // helper functions to isInElement
  bool within_tolerance( const double & val, const double & tol );
  double vector_norm_sq( const double *theVector );
  double parametric_distance( const double X, const double Y);
  double parametric_distance( const std::vector<double> &x);
};

// 2D Quad 4 subcontrol volume
class Quad2DSCV : public MasterElement
{
public:
  Quad2DSCV();
  virtual ~Quad2DSCV();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void quad_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);
};

// 2D Quad 4 subcontrol surface
class Quad2DSCS : public MasterElement
{
public:
  Quad2DSCS();
  virtual ~Quad2DSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void gij(
     const double *coords,
     double *gupperij,
     double *gij,
     double *deriv);

  const int * adjacentNodes();

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void quad_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);
  
  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);
  
  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);

  // extrusion data structure access
  const int * faceNodeOnExtrudedElem();
  const int * opposingNodeOnExtrudedElem();
  const int * faceScsIpOnExtrudedElem();
  const int * faceScsIpOnFaceEdges();
  const double * edgeAlignedArea();
};

class QuadrilateralP2Element : public MasterElement
{
public:
  QuadrilateralP2Element();
  virtual ~QuadrilateralP2Element() {}

  void shape_fcn(double *shpfc);
  void shifted_shape_fcn(double *shpfc);
protected:
  struct ContourData {
    Jacobian::Direction direction;
    double weight;
  };

  void set_quadrature_rule();
  void GLLGLL_quadrature_weights();

  int tensor_product_node_map(int i, int j) const;

  double gauss_point_location(
    int nodeOrdinal,
    int gaussPointOrdinal) const;

  double shifted_gauss_point_location(
    int nodeOrdinal,
    int gaussPointOrdinal) const;

  double tensor_product_weight(
    int s1Node, int s2Node,
    int s1Ip, int s2Ip) const;

  double tensor_product_weight(int s1Node, int s1Ip) const;

  void eval_shape_functions_at_ips();
  void eval_shape_functions_at_shifted_ips();

  void eval_shape_derivs_at_ips();
  void eval_shape_derivs_at_shifted_ips();

  void eval_shape_derivs_at_face_ips();

  const double scsDist_;
  bool useGLLGLL_;
  const int nodes1D_;
  int numQuad_;

  //quadrature info
  std::vector<double> gaussAbscissae1D_;
  std::vector<double> gaussAbscissae_;
  std::vector<double> gaussAbscissaeShift_;
  std::vector<double> gaussWeight_;

  std::vector<int> stkNodeMap_;
  std::vector<double> scsEndLoc_;

  std::vector<double> shapeFunctions_;
  std::vector<double> shapeFunctionsShift_;
  std::vector<double> shapeDerivs_;
  std::vector<double> shapeDerivsShift_;
  std::vector<double> expFaceShapeDerivs_;
private:
  void quad9_shape_fcn(
    int npts,
    const double *par_coord,
    double* shape_fcn
  ) const;

  void quad9_shape_deriv(
    int npts,
    const double *par_coord,
    double* shape_fcn
  ) const;
};

// 3D Quad 27 subcontrol volume
class Quad92DSCV : public QuadrilateralP2Element
{
public:
  Quad92DSCV();
  virtual ~Quad92DSCV() {}

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

private:
  void set_interior_info();

  double jacobian_determinant(
    const double *elemNodalCoords,
    const double *shapeDerivs ) const;

  std::vector<double> ipWeight_;
};

// 3D Hex 27 subcontrol surface
class Quad92DSCS : public QuadrilateralP2Element
{
public:
  Quad92DSCS();
  virtual ~Quad92DSCS() {}

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv);

  const int * adjacentNodes();

  const int * ipNodeMap(int ordinal = 0);

  int opposingNodes(
    const int ordinal, const int node);

  int opposingFace(
    const int ordinal, const int node);

private:
  void set_interior_info();
  void set_boundary_info();

  void area_vector(
    const Jacobian::Direction direction,
    const double *elemNodalCoords,
    double *shapeDeriv,
    double *areaVector ) const;

  std::vector<ContourData> ipInfo_;
  int ipsPerFace_;
};

// 2D Tri 3 subcontrol volume
class Tri2DSCV : public MasterElement
{
public:
  Tri2DSCV();
  virtual ~Tri2DSCV();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

};

// 2D Tri 3 subcontrol surface
class Tri2DSCS : public MasterElement
{
public:
  Tri2DSCS();
  virtual ~Tri2DSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void shifted_grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error );

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv);

  const int * adjacentNodes();

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);
  
  void tri_shape_fcn(
    const int &npts,
    const double *par_coord, 
    double* shape_fcn);

  int opposingNodes(
    const int ordinal, const int node);
  
  int opposingFace(
    const int ordinal, const int node);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);
  
  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  double tri_parametric_distance(
    const std::vector<double> &x);
  
  void general_face_grad_op(
    const int face_ordinal,
    const double *isoParCoord,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error );

  void sidePcoords_to_elemPcoords(
    const int & side_ordinal,
    const int & npoints,
    const double *side_pcoords,
    double *elem_pcoords);

};

// 3D Quad 4
class Quad3DSCS : public MasterElement
{
public:

  Quad3DSCS();
  virtual ~Quad3DSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);
  
  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  bool within_tol( const double & val, const double & tol );
  
  double vector_norm2( const double * vect, int len );

  void non_unit_face_normal(
    const double * par_coord,
    const double * elem_nodal_coor,
    double * normal_vector );
  
  double parametric_distance(const std::vector<double> &x);

  const double elemThickness_;
};

// 3D Quad 9
class Quad93DSCS : public HexahedralP2Element
{
public:
  Quad93DSCS();
  virtual ~Quad93DSCS() {}

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

private:
  void set_interior_info();
  void eval_shape_functions_at_ips() final;
  void eval_shape_derivs_at_ips() final;

  void eval_shape_functions_at_shifted_ips() final;
  void eval_shape_derivs_at_shifted_ips() final;

  void area_vector(
    const double *coords,
    const double *shapeDerivs,
    double *areaVector) const;

  void quad9_shape_fcn(
    int npts,
    const double *par_coord,
    double* shape_fcn
  ) const;

  void quad9_shape_deriv(
    int npts,
    const double *par_coord,
    double* shape_fcn
  ) const;

  std::vector<double> ipWeight_;
  const int surfaceDimension_;
};

// 3D Tri 3
class Tri3DSCS : public MasterElement
{
public:

  Tri3DSCS();
  virtual ~Tri3DSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
     double *shpfc);

   void shifted_shape_fcn(
     double *shpfc);

   void tri_shape_fcn(
     const int &npts,
     const double *par_coord,
     double* shape_fcn);

   double isInElement(
     const double *elemNodalCoord,
     const double *pointCoord,
     double *isoParCoord);

  double parametric_distance(
    const std::vector<double> &x);

  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

};

// edge 2d
class Edge2DSCS : public MasterElement
{
public:
  Edge2DSCS();
  virtual ~Edge2DSCS();

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

  double isInElement(
    const double *elemNodalCoord,
    const double *pointCoord,
    double *isoParCoord);
  
  void interpolatePoint(
    const int &nComp,
    const double *isoParCoord,
    const double *field,
    double *result);

  void general_shape_fcn(
    const int numIp,
    const double *isoParCoord,
    double *shpfc);

  double parametric_distance(const std::vector<double> &x);

  const double elemThickness_;  
};

// edge 2d
class Edge32DSCS : public QuadrilateralP2Element
{
public:
  Edge32DSCS();
  virtual ~Edge32DSCS() {}

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  void shape_fcn(
    double *shpfc);

  void shifted_shape_fcn(
    double *shpfc);

private:
  void area_vector(
    const double *coords,
    const double s,
    double *areaVector) const;

  std::vector<double> ipWeight_;
};

} // namespace nalu
} // namespace Sierra

#endif
