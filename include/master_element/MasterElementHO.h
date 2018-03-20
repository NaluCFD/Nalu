/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef MasterElementHO_h
#define MasterElementHO_h

#include <master_element/MasterElement.h>
#include <element_promotion/TensorProductQuadratureRule.h>
#include <element_promotion/LagrangeBasis.h>

#include <element_promotion/ElementDescription.h>
#include <element_promotion/HexNElementDescription.h>
#include <element_promotion/QuadNElementDescription.h>

#include <vector>
#include <array>

namespace sierra{
namespace nalu{

  struct ContourData {
    Jacobian::Direction direction;
    double weight;
  };

struct ElementDescription;
struct HexNElementDescription;

class LagrangeBasis;
class TensorProductQuadratureRule;

class HigherOrderHexSCV final: public MasterElement
{
public:
  HigherOrderHexSCV(
    ElementDescription elem,
    LagrangeBasis basis,
    TensorProductQuadratureRule quadrature);

  virtual ~HigherOrderHexSCV() {}

  void shape_fcn(double *shpfc) final;
  const int * ipNodeMap(int ordinal = 0) final;

  void determinant(
    const int nelem,
    const double *coords,
    double *volume,
    double * error ) final;

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error) final;

  std::vector<double> shape_functions() {
    return shapeFunctionVals_;
  }

  std::vector<double> shape_function_derivatives() {
    return shapeDerivs_;
  }

  std::vector<double> ip_weights() {
    return ipWeights_;
  }


private:
  void set_interior_info();

  double jacobian_determinant(
    const double* POINTER_RESTRICT elemNodalCoords,
    const double* POINTER_RESTRICT shapeDerivs ) const;

  const ElementDescription elem_;
  LagrangeBasis basis_;
  const TensorProductQuadratureRule quadrature_;

  std::vector<double> shapeFunctionVals_;
  std::vector<double> shapeDerivs_;
  std::vector<double> ipWeights_;
};

// 3D Hex 27 subcontrol surface
class HigherOrderHexSCS final: public MasterElement
{
public:
  HigherOrderHexSCS(
    ElementDescription elem,
    LagrangeBasis basis,
    TensorProductQuadratureRule quadrature);
  virtual ~HigherOrderHexSCS() {}

  void shape_fcn(double *shpfc) final;

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error) final;

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error) final;

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error) final;

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv) final;

  double isInElement(
      const double *elemNodalCoord,
      const double *pointCoord,
      double *isoParCoord) final;

  void interpolatePoint(
      const int &nComp,
      const double *isoParCoord,
      const double *field,
      double *result) final;

  const int * adjacentNodes() final;

  const int * ipNodeMap(int ordinal = 0) final;

  const int * side_node_ordinals(int ordinal = 0) final;

  int opposingNodes(
    const int ordinal, const int node) final;

  int opposingFace(
    const int ordinal, const int node) final;

  std::vector<double> shape_functions() {
    return shapeFunctionVals_;
  }

  std::vector<double> shape_function_derivatives() {
    return shapeDerivs_;
  }


private:
  void set_interior_info();
  void set_boundary_info();

  template <Jacobian::Direction direction> void
  area_vector(
    const double *POINTER_RESTRICT elemNodalCoords,
    double *POINTER_RESTRICT shapeDeriv,
    double *POINTER_RESTRICT areaVector) const;

  const ElementDescription elem_;
  LagrangeBasis basis_;
  const TensorProductQuadratureRule quadrature_;

  std::vector<double> shapeFunctionVals_;
  std::vector<double> shapeDerivs_;
  std::vector<double> expFaceShapeDerivs_;
  std::vector<ContourData> ipInfo_;
  int ipsPerFace_;
};

// 3D Quad 9
class HigherOrderQuad3DSCS final: public MasterElement
{
public:
  HigherOrderQuad3DSCS(
    ElementDescription elem,
    LagrangeBasis basis,
    TensorProductQuadratureRule quadrature);

  virtual ~HigherOrderQuad3DSCS() {}

  void shape_fcn(double *shpfc) final;

  const int * ipNodeMap(int ordinal = 0);

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error );

  std::vector<double> shape_functions() {
    return shapeFunctionVals_;
  }

  std::vector<double> shape_function_derivatives() {
    return shapeDerivs_;
  }

  std::vector<double> ip_weights() {
    return ipWeights_;
  }

private:
  void set_interior_info();
  void eval_shape_functions_at_ips();
  void eval_shape_derivs_at_ips();

  void area_vector(
    const double* POINTER_RESTRICT elemNodalCoords,
    const double* POINTER_RESTRICT shapeDeriv,
    std::array<double,3>& areaVector) const;

  const ElementDescription elem_;
  LagrangeBasis basis_;
  const TensorProductQuadratureRule quadrature_;

  std::vector<double> shapeFunctionVals_;
  std::vector<double> shapeDerivs_;
  std::vector<double> ipWeights_;
  int surfaceDimension_;
};

class HigherOrderQuad2DSCV final: public MasterElement
{
public:
  HigherOrderQuad2DSCV(
    ElementDescription elem,
    LagrangeBasis basis,
    TensorProductQuadratureRule quadrature);
  virtual ~HigherOrderQuad2DSCV() {}

  void shape_fcn(double *shpfc) final;

  const int * ipNodeMap(int ordinal = 0) final;

  void determinant(
    const int nelem,
    const double *coords,
    double *volume,
    double * error ) final;

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error) final;

  std::vector<double> shape_functions() {
    return shapeFunctionVals_;
  }

  std::vector<double> shape_function_derivatives() {
    return shapeDerivs_;
  }

  std::vector<double> ip_weights() {
    return ipWeights_;
  }

private:
  void set_interior_info();

  double jacobian_determinant(
    const double* POINTER_RESTRICT elemNodalCoords,
    const double* POINTER_RESTRICT shapeDerivs ) const;

  const ElementDescription elem_;
  LagrangeBasis basis_;
  const TensorProductQuadratureRule quadrature_;

  std::vector<double> shapeFunctionVals_;
  std::vector<double> shapeDerivs_;
  std::vector<double> ipWeights_;
};
class HigherOrderQuad2DSCS final: public MasterElement
{
public:
  HigherOrderQuad2DSCS(
    ElementDescription elem,
    LagrangeBasis basis,
    TensorProductQuadratureRule quadrature);
  virtual ~HigherOrderQuad2DSCS() {}

  void shape_fcn(double *shpfc) final;

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error) final;

  void grad_op(
    const int nelem,
    const double *coords,
    double *gradop,
    double *deriv,
    double *det_j,
    double * error) final;

  void face_grad_op(
    const int nelem,
    const int face_ordinal,
    const double *coords,
    double *gradop,
    double *det_j,
    double * error) final;

  void gij(
    const double *coords,
    double *gupperij,
    double *glowerij,
    double *deriv) final;

  double isInElement(
      const double *elemNodalCoord,
      const double *pointCoord,
      double *isoParCoord) final;

  void interpolatePoint(
      const int &nComp,
      const double *isoParCoord,
      const double *field,
      double *result) final;

  const int * adjacentNodes() final;

  const int * ipNodeMap(int ordinal = 0) final;

  int opposingNodes(
    const int ordinal, const int node) final;

  int opposingFace(
    const int ordinal, const int node) final;

  const int * side_node_ordinals(int ordinal = 0) final;

  std::vector<double> shape_functions() {
    return shapeFunctionVals_;
  }

  std::vector<double> shape_function_derivatives() {
    return shapeDerivs_;
  }

private:
  void set_interior_info();
  void set_boundary_info();

  template <Jacobian::Direction direction> void
  area_vector(
    const double *POINTER_RESTRICT elemNodalCoords,
    double *POINTER_RESTRICT shapeDeriv,
    double *POINTER_RESTRICT normalVec ) const;

  const ElementDescription elem_;
  LagrangeBasis basis_;
  const TensorProductQuadratureRule quadrature_;

  std::vector<double> shapeFunctionVals_;
  std::vector<double> shapeDerivs_;
  std::vector<ContourData> ipInfo_;
  int ipsPerFace_;
  std::vector<double> expFaceShapeDerivs_;
};

class HigherOrderEdge2DSCS final: public MasterElement
{
public:
  explicit HigherOrderEdge2DSCS(
    ElementDescription elem,
    LagrangeBasis basis,
    TensorProductQuadratureRule quadrature);
  virtual ~HigherOrderEdge2DSCS() {}

  const int * ipNodeMap(int ordinal = 0) final;

  void determinant(
    const int nelem,
    const double *coords,
    double *areav,
    double * error ) final;

  void shape_fcn(
    double *shpfc) final;

  std::vector<double> shape_functions() {
    return shapeFunctionVals_;
  }

  std::vector<double> shape_function_derivatives() {
    return shapeDerivs_;
  }

  std::vector<double> ip_weights() {
    return ipWeights_;
  }

private:
  void area_vector(
    const double* POINTER_RESTRICT elemNodalCoords,
    const double* POINTER_RESTRICT shapeDeriv,
    std::array<double,2>& areaVector) const;

  const ElementDescription elem_;
  LagrangeBasis basis_;
  const TensorProductQuadratureRule quadrature_;

  std::vector<double> shapeFunctionVals_;
  std::vector<double> shapeDerivs_;
  std::vector<double> ipWeights_;
};

} // namespace nalu
} // namespace Sierra

#endif
