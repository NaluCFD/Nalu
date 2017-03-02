#include <gtest/gtest.h>
#include <limits>
#include <random>
#include <stdexcept>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>

#include <element_promotion/ElementDescription.h>
#include <element_promotion/MasterElementHO.h>
#include <element_promotion/QuadratureRule.h>
#include <element_promotion/LagrangeBasis.h>
#include <element_promotion/TensorProductQuadratureRule.h>
#include <element_promotion/QuadratureKernels.h>

#include "UnitTestUtils.h"

namespace {

double poly_val(std::vector<double> coeffs, double x)
{
  double val = 0.0;
  for (unsigned j = 0; j < coeffs.size(); ++j) {
    val += coeffs[j]*std::pow(x,j);
  }
  return val;
}
//--------------------------------------------------------------------------
double poly_der(std::vector<double> coeffs, double x)
{
  double val = 0.0;
  for (unsigned j = 1; j < coeffs.size(); ++j) {
    val += coeffs[j]*std::pow(x,j-1)*j;
  }
  return val;
}
//--------------------------------------------------------------------------
double poly_int(std::vector<double> coeffs,
  double xlower, double xupper)
{
  double upper = 0.0; double lower = 0.0;
  for (unsigned j = 0; j < coeffs.size(); ++j) {
    upper += coeffs[j]*std::pow(xupper,j+1)/(j+1.0);
    lower += coeffs[j]*std::pow(xlower,j+1)/(j+1.0);
  }
  return (upper-lower);
}
//--------------------------------------------------------------------------
void
check_interpolation_quad(int polyOrder, int numIps, double tol)
{
  // create a (-1,1) x (-1,1) element filled with polynomial values
  // and interpolate to a vector of random points
  // (between -1.05 and 1.05 to check edges)

  auto elemDesc = sierra::nalu::ElementDescription::create(2, polyOrder);
  auto basis = sierra::nalu::LagrangeBasis(elemDesc->inverseNodeMap, elemDesc->nodeLocs1D);

  std::mt19937 rng;
  rng.seed(0);
  std::uniform_real_distribution<double> coeff(-1.0, 1.0);
  std::uniform_real_distribution<double> loc(-1.05, 1.05);

  int nodesPerElement = elemDesc->nodesPerElement;
  std::vector<double> intgLoc(numIps * elemDesc->dimension);
  std::vector<double> coeffsX(elemDesc->polyOrder+1);
  std::vector<double> coeffsY(elemDesc->polyOrder+1);
  std::vector<double> nodalValues(nodesPerElement);
  std::vector<double> interpWeights(numIps * nodesPerElement);
  std::vector<double> exactInterp(numIps);
  std::vector<double> approxInterp(numIps);

  for (int ip = 0; ip < numIps; ++ip) {
    int offset = ip * elemDesc->dimension;
    intgLoc[offset + 0] = loc(rng);
    intgLoc[offset + 1] = loc(rng);
  }

  for (int k = 0; k < elemDesc->polyOrder+1; ++k) {
    coeffsX[k] = coeff(rng);
    coeffsY[k] = coeff(rng);
  }

  interpWeights = basis.eval_basis_weights(intgLoc);

  for (int i = 0; i < elemDesc->nodes1D; ++i) {
    for (int j = 0; j < elemDesc->nodes1D; ++j) {
      nodalValues[elemDesc->node_map(i, j)] =
            poly_val(coeffsX, elemDesc->nodeLocs1D[i])
          * poly_val(coeffsY, elemDesc->nodeLocs1D[j]);
    }
  }

  int offset_1 = 0;
  for (int ip = 0; ip < numIps; ++ip) {
    exactInterp[ip] =
          poly_val(coeffsX, intgLoc[offset_1 + 0])
        * poly_val(coeffsY, intgLoc[offset_1 + 1]);
    offset_1 += 2;
  }

  for (int ip = 0; ip < numIps; ++ip) {
    int ip_offset = ip * nodesPerElement;
    double val = 0.0;
    for (int nodeNumber = 0; nodeNumber < nodesPerElement; ++nodeNumber) {
      int offset = ip_offset + nodeNumber;
      val += interpWeights[offset] * nodalValues[nodeNumber];
    }
    approxInterp[ip] = val;
  }

  for (unsigned j = 0; j < exactInterp.size(); ++j) {
    EXPECT_NEAR(approxInterp[j], exactInterp[j], tol);
  }
}
//--------------------------------------------------------------------------
void
check_interpolation_hex(int polyOrder, int numIps, double tol)
{
  // create a (-1,1) x (-1,1) element filled with polynomial values
  // and interpolate to a vector of random points
  // (between -1.05 and 1.05 to check edges)
  auto elemDesc = sierra::nalu::ElementDescription::create(3, polyOrder);
  auto basis = sierra::nalu::LagrangeBasis(elemDesc->inverseNodeMap, elemDesc->nodeLocs1D);

  std::mt19937 rng;
  rng.seed(0);
  std::uniform_real_distribution<double> coeff(-1.0, 1.0);
  std::uniform_real_distribution<double> loc(-1.05, 1.05);

  int nodesPerElement = elemDesc->nodesPerElement;
  std::vector<double> intgLoc(numIps * elemDesc->dimension);
  std::vector<double> coeffsX(elemDesc->polyOrder+1);
  std::vector<double> coeffsY(elemDesc->polyOrder+1);
  std::vector<double> coeffsZ(elemDesc->polyOrder+1);
  std::vector<double> nodalValues(nodesPerElement);
  std::vector<double> interpWeights(numIps * nodesPerElement);
  std::vector<double> exactInterp(numIps);
  std::vector<double> approxInterp(numIps);

  for (int ip = 0; ip < numIps; ++ip) {
    int offset = ip * elemDesc->dimension;
    intgLoc[offset + 0] = loc(rng);
    intgLoc[offset + 1] = loc(rng);
    intgLoc[offset + 2] = loc(rng);
  }

  for (int k = 0; k < elemDesc->polyOrder+1; ++k) {
    coeffsX[k] = coeff(rng);
    coeffsY[k] = coeff(rng);
    coeffsZ[k] = coeff(rng);
  }

  interpWeights = basis.eval_basis_weights(intgLoc);

  for (int i = 0; i < elemDesc->nodes1D; ++i) {
    for (int j = 0; j < elemDesc->nodes1D; ++j) {
      for (int k = 0; k < elemDesc->nodes1D; ++k) {
        nodalValues[elemDesc->node_map(i, j, k)] =
            poly_val(coeffsX, elemDesc->nodeLocs1D[i])
            * poly_val(coeffsY, elemDesc->nodeLocs1D[j])
            * poly_val(coeffsZ, elemDesc->nodeLocs1D[k]);
      }
    }
  }

  for (int ip = 0; ip < numIps; ++ip) {
    int offset = ip*elemDesc->dimension;
    exactInterp[ip] =
          poly_val(coeffsX, intgLoc[offset + 0])
        * poly_val(coeffsY, intgLoc[offset + 1])
        * poly_val(coeffsZ, intgLoc[offset + 2]);
  }

  for (int ip = 0; ip < numIps; ++ip) {
    int ip_offset = ip * nodesPerElement;
    double val = 0.0;
    for (int nodeNumber = 0; nodeNumber < nodesPerElement; ++nodeNumber) {
      int offset = ip_offset + nodeNumber;
      val += interpWeights[offset] * nodalValues[nodeNumber];
    }
    approxInterp[ip] = val;
  }

  for (unsigned j = 0; j < exactInterp.size(); ++j) {
    EXPECT_NEAR(approxInterp[j], exactInterp[j], tol);
  }
}
//--------------------------------------------------------------------------
void
check_derivative_quad(int polyOrder, int numIps, double tol)
{
  // create a (-1,1) x (-1,1) element filled with polynomial values
  // and compute derivatives at a vector of random points
  // (between -1.05 and 1.05 to check edges)
  auto elemDesc = sierra::nalu::ElementDescription::create(2, polyOrder);
  auto basis = sierra::nalu::LagrangeBasis(elemDesc->inverseNodeMap, elemDesc->nodeLocs1D);

  std::mt19937 rng;
  rng.seed(0);
  std::uniform_real_distribution<double> coeff(-1.0,1.0);
  std::uniform_real_distribution<double> loc(-1.05,1.05);

  std::vector<double> intgLoc(numIps*elemDesc->dimension);
  std::vector<double> coeffsX(elemDesc->polyOrder+1);
  std::vector<double> coeffsY(elemDesc->polyOrder+1);
  std::vector<double> exactDeriv(numIps*elemDesc->dimension);
  std::vector<double> approxDeriv(numIps*elemDesc->dimension);

  for (int ip = 0; ip < numIps; ++ip) {
    int offset = ip*elemDesc->dimension;
    intgLoc[offset+0] = loc(rng);
    intgLoc[offset+1] = loc(rng);
  }

  std::vector<double> diffWeights = basis.eval_deriv_weights(intgLoc);

  for (int k = 0; k < elemDesc->polyOrder+1; ++k) {
    coeffsX[k] = coeff(rng);
    coeffsY[k] = coeff(rng);
  }

  // create a (-1,1) x (-1,1) element and fill it with polynomial values
  // expect exact values to floating-point precision
  std::vector<double> nodalValues(elemDesc->nodesPerElement);
  for (int i = 0; i < elemDesc->nodes1D; ++i) {
    for(int j = 0; j < elemDesc->nodes1D; ++j) {
      nodalValues[elemDesc->node_map(i,j)] =
          poly_val(coeffsX,elemDesc->nodeLocs1D[i]) * poly_val(coeffsY,elemDesc->nodeLocs1D[j]);
    }
  }

  for (int ip = 0; ip < numIps; ++ip) {
    int offset = ip*elemDesc->dimension;
    exactDeriv[offset + 0] =
          poly_der(coeffsX, intgLoc[offset])
        * poly_val(coeffsY, intgLoc[offset + 1]);
    exactDeriv[offset + 1] =
          poly_val(coeffsX, intgLoc[offset])
        * poly_der(coeffsY, intgLoc[offset + 1]);
  }

  for (int ip = 0; ip < numIps; ++ip) {
    double dndx = 0.0;
    double dndy = 0.0;
    for (int node = 0; node < elemDesc->nodesPerElement; ++node) {
      int deriv_offset = (ip*elemDesc->nodesPerElement+node)*elemDesc->dimension;
      dndx += diffWeights[deriv_offset + 0] * nodalValues[node];
      dndy += diffWeights[deriv_offset + 1] * nodalValues[node];
    }
    approxDeriv[ip*elemDesc->dimension+0] = dndx;
    approxDeriv[ip*elemDesc->dimension+1] = dndy;
  }

  for (unsigned j = 0; j < exactDeriv.size(); ++j) {
    EXPECT_NEAR(approxDeriv[j], exactDeriv[j], tol);
  }
}
//--------------------------------------------------------------------------
void
check_derivative_hex(int polyOrder, int numIps, double tol)
{
  // create a (-1,1) x (-1,1) element filled with polynomial values
  // and compute derivatives at a vector of random points
  // (between -1.05 and 1.05 to check edges)
  auto elemDesc = sierra::nalu::ElementDescription::create(3, polyOrder);
  auto basis = sierra::nalu::LagrangeBasis(elemDesc->inverseNodeMap, elemDesc->nodeLocs1D);

  std::mt19937 rng;
  rng.seed(0);
  std::uniform_real_distribution<double> coeff(-1.0,1.0);
  std::uniform_real_distribution<double> loc(-1.05,1.05);

  std::vector<double> intgLoc(numIps*elemDesc->dimension);
  std::vector<double> coeffsX(elemDesc->polyOrder+1);
  std::vector<double> coeffsY(elemDesc->polyOrder+1);
  std::vector<double> coeffsZ(elemDesc->polyOrder+1);
  std::vector<double> exactDeriv(numIps*elemDesc->dimension);
  std::vector<double> approxDeriv(numIps*elemDesc->dimension);

  for (int ip = 0; ip < numIps; ++ip) {
    int offset = ip*elemDesc->dimension;
    intgLoc[offset+0] = loc(rng);
    intgLoc[offset+1] = loc(rng);
    intgLoc[offset+2] = loc(rng);
  }

  std::vector<double> diffWeights = basis.eval_deriv_weights(intgLoc);

  for (int k = 0; k < elemDesc->polyOrder+1; ++k) {
    coeffsX[k] = coeff(rng);
    coeffsY[k] = coeff(rng);
    coeffsZ[k] = coeff(rng);
  }

  std::vector<double> nodalValues(elemDesc->nodesPerElement);
  for (int i = 0; i < elemDesc->nodes1D; ++i) {
    for(int j = 0; j < elemDesc->nodes1D; ++j) {
      for(int k = 0; k < elemDesc->nodes1D; ++k) {
        nodalValues[elemDesc->node_map(i,j,k)] =
              poly_val(coeffsX,elemDesc->nodeLocs1D[i])
            * poly_val(coeffsY,elemDesc->nodeLocs1D[j])
            * poly_val(coeffsZ,elemDesc->nodeLocs1D[k]);
      }
    }
  }

  for (int ip = 0; ip < numIps; ++ip) {
    int offset = ip*elemDesc->dimension;
    exactDeriv[offset + 0] =
          poly_der(coeffsX, intgLoc[offset])
        * poly_val(coeffsY, intgLoc[offset + 1])
        * poly_val(coeffsZ, intgLoc[offset + 2]);
    exactDeriv[offset + 1] =
          poly_val(coeffsX, intgLoc[offset])
        * poly_der(coeffsY, intgLoc[offset + 1])
        * poly_val(coeffsZ, intgLoc[offset + 2]);
    exactDeriv[offset + 2] =
          poly_val(coeffsX, intgLoc[offset])
        * poly_val(coeffsY, intgLoc[offset + 1])
        * poly_der(coeffsZ, intgLoc[offset + 2]);
  }

  for (int ip = 0; ip < numIps; ++ip) {
    double dndx = 0.0;
    double dndy = 0.0;
    double dndz = 0.0;
    for (int node = 0; node < elemDesc->nodesPerElement; ++node) {
      int deriv_offset = (ip*elemDesc->nodesPerElement+node)*elemDesc->dimension;
      dndx += diffWeights[deriv_offset + 0] * nodalValues[node];
      dndy += diffWeights[deriv_offset + 1] * nodalValues[node];
      dndz += diffWeights[deriv_offset + 2] * nodalValues[node];
    }
    approxDeriv[ip*elemDesc->dimension+0] = dndx;
    approxDeriv[ip*elemDesc->dimension+1] = dndy;
    approxDeriv[ip*elemDesc->dimension+2] = dndz;
  }

  for (unsigned j = 0; j < exactDeriv.size(); ++j) {
    EXPECT_NEAR(approxDeriv[j], exactDeriv[j], tol);
  }
}
//--------------------------------------------------------------------------
void
check_volume_quadrature_quad(int polyOrder, double tol)
{
  // create a (-1,1) x (-1,1) element filled with polynomial values
  // and integrate the polynomial over the dual nodal volumes

  std::mt19937 rng;
  rng.seed(0);
  std::uniform_real_distribution<double> coeff(-1.0, 1.0);

  auto elemDesc = sierra::nalu::ElementDescription::create(2, polyOrder);
  auto basis = sierra::nalu::LagrangeBasis(elemDesc->inverseNodeMap, elemDesc->nodeLocs1D);
  auto quad = sierra::nalu::TensorProductQuadratureRule("GaussLegendre", elemDesc->polyOrder);
  auto masterElement = sierra::nalu::HigherOrderQuad2DSCV(*elemDesc, basis, quad);

  const auto& interpWeights  = masterElement.shape_functions();
  const auto& ipWeights = masterElement.ip_weights();
  const auto* ipNodeMap = masterElement.ipNodeMap();
  const auto& scsEndLoc = quad.scs_end_loc();
  std::vector<double> approxInt(elemDesc->nodesPerElement);
  std::vector<double> coeffsX(elemDesc->polyOrder+1);
  std::vector<double> coeffsY(elemDesc->polyOrder+1);
  std::vector<double> exactInt(elemDesc->nodesPerElement);
  std::vector<double> nodalValues(elemDesc->nodesPerElement);

  for (int k = 0; k < elemDesc->polyOrder+1; ++k) {
    coeffsX[k] = coeff(rng);
    coeffsY[k] = coeff(rng);
  }

  for (int i = 0; i < elemDesc->nodes1D; ++i) {
    for (int j = 0; j < elemDesc->nodes1D; ++j) {
      nodalValues[elemDesc->node_map(i, j)] =
            poly_val(coeffsX, elemDesc->nodeLocs1D[i])
          * poly_val(coeffsY, elemDesc->nodeLocs1D[j]);

      exactInt[elemDesc->node_map(i, j)] =
            poly_int(coeffsX, scsEndLoc[i], scsEndLoc[i + 1])
          * poly_int(coeffsY, scsEndLoc[j], scsEndLoc[j + 1]);
    }
  }

  approxInt.assign(approxInt.size(), 0.0);
  for (int ip = 0; ip < masterElement.numIntPoints_; ++ip) {
    double interpValue = 0.0;
    for (int nodeNumber = 0; nodeNumber < elemDesc->nodesPerElement; ++nodeNumber) {
      interpValue += interpWeights[ip*elemDesc->nodesPerElement+nodeNumber] * nodalValues[nodeNumber];
    }
    approxInt[ipNodeMap[ip]] += ipWeights[ip] * interpValue; //ipweights -> ws_scv_volume if not square
  }

  for (unsigned j = 0; j < exactInt.size(); ++j) {
    EXPECT_NEAR(approxInt[j], exactInt[j], tol);
  }
}
  //--------------------------------------------------------------------------
void
check_volume_quadrature_hex(int polyOrder, double tol)
{
  // create a (-1,1) x (-1,1) element filled with polynomial values
  // and integrate the polynomial over the dual nodal volumes

  std::mt19937 rng;
  rng.seed(0);
  std::uniform_real_distribution<double> coeff(-10.0, 10.0);

  auto elemDesc = sierra::nalu::ElementDescription::create(3, polyOrder);
  auto basis = sierra::nalu::LagrangeBasis(elemDesc->inverseNodeMap, elemDesc->nodeLocs1D);
  auto quad = sierra::nalu::TensorProductQuadratureRule("GaussLegendre", elemDesc->polyOrder);
  auto masterElement = sierra::nalu::HigherOrderHexSCV(*elemDesc, basis, quad);

  const auto& interpWeights = masterElement.shape_functions();
  const auto& ipWeights = masterElement.ip_weights();
  const auto* ipNodeMap = masterElement.ipNodeMap();
  const auto& scsEndLoc = quad.scs_end_loc();
  std::vector<double> approxInt(elemDesc->nodesPerElement);
  std::vector<double> coeffsX(elemDesc->polyOrder+1);
  std::vector<double> coeffsY(elemDesc->polyOrder+1);
  std::vector<double> coeffsZ(elemDesc->polyOrder+1);
  std::vector<double> exactInt(elemDesc->nodesPerElement);
  std::vector<double> nodalValues(elemDesc->nodesPerElement);

  for (int k = 0; k < elemDesc->polyOrder+1; ++k) {
    coeffsX[k] = coeff(rng);
    coeffsY[k] = coeff(rng);
    coeffsZ[k] = coeff(rng);
  }

  for (int i = 0; i < elemDesc->nodes1D; ++i) {
    for (int j = 0; j < elemDesc->nodes1D; ++j) {
      for (int k = 0; k < elemDesc->nodes1D; ++k) {
        nodalValues[elemDesc->node_map(i, j, k)] =
            poly_val(coeffsX, elemDesc->nodeLocs1D[i])
            * poly_val(coeffsY, elemDesc->nodeLocs1D[j])
            * poly_val(coeffsZ, elemDesc->nodeLocs1D[k]);

        exactInt[elemDesc->node_map(i, j, k)] =
              poly_int(coeffsX, scsEndLoc[i], scsEndLoc[i + 1])
            * poly_int(coeffsY, scsEndLoc[j], scsEndLoc[j + 1])
            * poly_int(coeffsZ, scsEndLoc[k], scsEndLoc[k + 1]);
      }
    }
  }

  approxInt.assign(approxInt.size(), 0.0);
  for (int ip = 0; ip < masterElement.numIntPoints_; ++ip) {
    double interpValue = 0.0;
    for (int nodeNumber = 0; nodeNumber < elemDesc->nodesPerElement; ++nodeNumber) {
      interpValue += interpWeights[ip*elemDesc->nodesPerElement+nodeNumber] * nodalValues[nodeNumber];
    }
    approxInt[ipNodeMap[ip]] +=  ipWeights[ip]*interpValue;
  }

  for (unsigned j = 0; j < exactInt.size(); ++j) {
    EXPECT_NEAR(approxInt[j], exactInt[j], tol);
  }
}
//--------------------------------------------------------------------------
void
check_volume_quadrature_quad_SGL(int polyOrder, double tol)
{
  // create a (-1,1) x (-1,1) element filled with polynomial values
  // and integrate the polynomial over the dual nodal volumes
  std::mt19937 rng;
  rng.seed(0);
  std::uniform_real_distribution<double> coeff(-10.0, 10.0);

  auto elemDesc = sierra::nalu::ElementDescription::create(2, polyOrder);
  auto basis = sierra::nalu::LagrangeBasis(elemDesc->inverseNodeMap, elemDesc->nodeLocs1D);
  auto quad = sierra::nalu::TensorProductQuadratureRule("SGL", elemDesc->polyOrder);
  auto masterElement = sierra::nalu::HigherOrderQuad2DSCV(*elemDesc, basis, quad);

  const auto& scsEndLoc = quad.scs_end_loc();
  std::vector<double> approxInt(elemDesc->nodesPerElement);
  std::vector<double> coeffsX(elemDesc->polyOrder+1);
  std::vector<double> coeffsY(elemDesc->polyOrder+1);
  std::vector<double> exactInt(elemDesc->nodesPerElement);
  std::vector<double> nodalValues(elemDesc->nodesPerElement);
  std::vector<double> nodalValuesTensor(elemDesc->nodesPerElement);
  std::vector<double> approxIntTensor(elemDesc->nodesPerElement, 0.0);

  auto quadOp = sierra::nalu::SGLQuadratureOps(*elemDesc);
  int nodes1D = elemDesc->nodes1D;

  std::vector<double> temp(elemDesc->nodesPerElement,0.0);

  // get a random polyOrder-degree polynomial
  for (int k = 0; k < elemDesc->polyOrder+1; ++k) {
    coeffsX[k] = coeff(rng);
    coeffsY[k] = coeff(rng);
  }

  // exact solution
  for (int i = 0; i < elemDesc->nodes1D; ++i) {
    for (int j = 0; j < elemDesc->nodes1D; ++j) {
      nodalValues[elemDesc->node_map(i, j)] =
            poly_val(coeffsX, elemDesc->nodeLocs1D[i])
          * poly_val(coeffsY, elemDesc->nodeLocs1D[j]);

      exactInt[elemDesc->node_map(i, j)] =
            poly_int(coeffsX, scsEndLoc[i], scsEndLoc[i + 1])
          * poly_int(coeffsY, scsEndLoc[j], scsEndLoc[j + 1]);
    }
  }

  approxIntTensor.assign(nodes1D*nodes1D,0.0);

  for (int i = 0; i < elemDesc->nodes1D; ++i) {
    for (int j = 0; j < elemDesc->nodes1D; ++j) {
      // this algorithm requires the nodes to be ordered like a tensor
      nodalValuesTensor[i+nodes1D*j] =  nodalValues[elemDesc->node_map(i, j)];

      //multiply by det(J)_ij here if not a square domain
    }
  }

  quadOp.volume_2D(nodalValuesTensor.data(),approxIntTensor.data());

  // convert back to tensor-product form
  for (int i = 0; i < elemDesc->nodes1D; ++i) {
    for (int j = 0; j < elemDesc->nodes1D; ++j) {
      approxInt[elemDesc->node_map(i,j)] = approxIntTensor[i+nodes1D*j];
    }
  }

  for (unsigned j = 0; j < exactInt.size(); ++j) {
    EXPECT_NEAR(approxInt[j], exactInt[j], tol);
  }
}
//--------------------------------------------------------------------------
void
check_volume_quadrature_hex_SGL(int polyOrder, double tol)
{
  // create a (-1,1) x (-1,1) element filled with polynomial values
  // and integrate the polynomial over the dual nodal volumes

  std::mt19937 rng;
  rng.seed(0);
  std::uniform_real_distribution<double> coeff(-10.0, 10.0);

  auto elemDesc = sierra::nalu::ElementDescription::create(3, polyOrder);
  auto basis = sierra::nalu::LagrangeBasis(elemDesc->inverseNodeMap, elemDesc->nodeLocs1D);
  auto quad = sierra::nalu::TensorProductQuadratureRule("SGL", elemDesc->polyOrder);
  auto masterElement = sierra::nalu::HigherOrderHexSCV(*elemDesc, basis, quad);

  const auto& scsEndLoc = quad.scs_end_loc();
  std::vector<double> approxInt(elemDesc->nodesPerElement);
  std::vector<double> coeffsX(elemDesc->polyOrder+1);
  std::vector<double> coeffsY(elemDesc->polyOrder+1);
  std::vector<double> coeffsZ(elemDesc->polyOrder+1);
  std::vector<double> exactInt(elemDesc->nodesPerElement);
  std::vector<double> nodalValues(elemDesc->nodesPerElement);
  std::vector<double> nodalValuesTensor(elemDesc->nodesPerElement);
  std::vector<double> approxIntTensor(elemDesc->nodesPerElement, 0.0);
  int nodes1D = elemDesc->nodes1D;
  int nodes2D = nodes1D*nodes1D;

  auto quadOp = sierra::nalu::SGLQuadratureOps(*elemDesc);
  std::vector<double> temp1(nodes2D, 0.0);
  std::vector<double> temp2(nodes2D, 0.0);

  // get a random polyOrder-degree polynomial
  for (int k = 0; k < elemDesc->polyOrder+1; ++k) {
    coeffsX[k] = coeff(rng);
    coeffsY[k] = coeff(rng);
    coeffsZ[k] = coeff(rng);
  }

  // exact solution
  for (int i = 0; i < elemDesc->nodes1D; ++i) {
    for (int j = 0; j < elemDesc->nodes1D; ++j) {
      for (int k = 0; k < elemDesc->nodes1D; ++k) {
        nodalValues[elemDesc->node_map(i, j, k)] =
              poly_val(coeffsX, elemDesc->nodeLocs1D[i])
            * poly_val(coeffsY, elemDesc->nodeLocs1D[j])
            * poly_val(coeffsZ, elemDesc->nodeLocs1D[k]);

        exactInt[elemDesc->node_map(i, j, k)] =
              poly_int(coeffsX, scsEndLoc[i], scsEndLoc[i + 1])
            * poly_int(coeffsY, scsEndLoc[j], scsEndLoc[j + 1])
            * poly_int(coeffsZ, scsEndLoc[k], scsEndLoc[k + 1]);
      }
    }
  }

  approxIntTensor.assign(nodes1D * nodes1D * nodes1D, 0.0);

  for (int i = 0; i < elemDesc->nodes1D; ++i) {
    for (int j = 0; j < elemDesc->nodes1D; ++j) {
      for (int k = 0; k < elemDesc->nodes1D; ++k) {
        // this algorithm requires the nodes to be ordered like a tensor
        nodalValuesTensor[i + nodes1D * (j + nodes1D * k)] =
            nodalValues[elemDesc->node_map(i, j, k)];

        //multiply by det(J)_ijk here if not a square domain
      }
    }
  }

  quadOp.volume_3D(nodalValues.data(), approxInt.data());

  for (unsigned j = 0; j < exactInt.size(); ++j) {
    EXPECT_NEAR(approxInt[j], exactInt[j], tol);
  }
}


}//namespace

constexpr int MAX_POLY_ORDER = 6;

TEST(HOMasterElements, interp_quad)
{
  double tol = 1.0e-10;
  int numIps = 10;

  for (int p = 1; p < MAX_POLY_ORDER; ++p) {
    check_interpolation_quad(p, numIps, tol);
  }
}

TEST(HOMasterElements, interp_hex)
{
  double tol = 1.0e-10;
  int numIps = 10;

  for (int p = 1; p < MAX_POLY_ORDER; ++p) {
    check_interpolation_hex(p, numIps, tol);
  }
}

TEST(HOMasterElements, deriv_quad)
{
  double tol = 1.0e-10;
  int numIps = 10;

  for (int p = 1; p < MAX_POLY_ORDER; ++p) {
    check_derivative_quad(p, numIps, tol);
  }
}

TEST(HOMasterElements, deriv_hex)
{
  double tol = 1.0e-10;
  int numIps = 10;

  for (int p = 1; p < MAX_POLY_ORDER; ++p) {
    check_derivative_hex(p, numIps, tol);
  }
}

TEST(HOMasterElements, integ_GL_quad)
{
  double tol = 1.0e-10;

  for (int p = 1; p < MAX_POLY_ORDER; ++p) {
    check_volume_quadrature_quad(p,  tol);
  }
}

TEST(HOMasterElements, integ_GL_hex)
{
  double tol = 1.0e-10;

  for (int p = 1; p < MAX_POLY_ORDER; ++p) {
    check_volume_quadrature_hex(p,  tol);
  }
}

TEST(HOMasterElements, integ_SGL_quad)
{
  double tol = 1.0e-10;

  for (int p = 1; p < MAX_POLY_ORDER; ++p) {
    check_volume_quadrature_quad_SGL(p,  tol);
  }
}

TEST(HOMasterElements, integ_SGL_hex)
{
  double tol = 1.0e-10;

  for (int p = 1; p < MAX_POLY_ORDER; ++p) {
    check_volume_quadrature_hex_SGL(p,  tol);
  }
}


