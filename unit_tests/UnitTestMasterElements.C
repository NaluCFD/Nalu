#include <gtest/gtest.h>
#include <limits>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <master_element/MasterElement.h>

#include <memory>
#include <random>

#include "UnitTestUtils.h"

namespace {

using VectorFieldType = stk::mesh::Field<double, stk::mesh::Cartesian>;
//-------------------------------------------------------------------------
double linear_scalar_value(int dim, double a, const double* b, const double* x)
{
  if (dim == 2u) {
    return (a + b[0] * x[0] + b[1] * x[1]);
  }
  return (a + b[0] * x[0] + b[1] * x[1] + b[2] * x[2]);
}
//-------------------------------------------------------------------------
void check_interpolation_at_ips(
  const stk::mesh::Entity* node_rels,
  const VectorFieldType& coordField,
  sierra::nalu::MasterElement& me)
{
  // Check that we can interpolate a random 3D polynomial
  // to the integration points
  int dim = me.nDim_;

  std::mt19937 rng;
  rng.seed(0);
  std::uniform_real_distribution<double> coeff(-1.0, 1.0);
  std::vector<double> coeffs(dim);

  double a = coeff(rng);
  for (int j = 0; j < dim; ++j) {
    coeffs[j] = coeff(rng);
  }

  const auto& intgLoc = me.intgLoc_;
  std::vector<double> polyResult(me.numIntPoints_);
  for (int j = 0; j < me.numIntPoints_; ++j) {
    polyResult[j] = linear_scalar_value(dim, a, coeffs.data(), &intgLoc[j*dim]);
  }

  std::vector<double> ws_field(me.nodesPerElement_);
  for (int j = 0; j < me.nodesPerElement_; ++j) {
    ws_field[j] = linear_scalar_value(dim, a, coeffs.data(), stk::mesh::field_data(coordField, node_rels[j]));
  }

  std::vector<double> meResult(me.numIntPoints_, 0.0);

  std::vector<double> meShapeFunctions(me.nodesPerElement_ * me.numIntPoints_);
  me.shape_fcn(meShapeFunctions.data());

  for (int j = 0; j < me.numIntPoints_; ++j) {
    for (int i = 0; i < me.nodesPerElement_; ++i) {
      meResult[j] += meShapeFunctions[j*me.nodesPerElement_+i] * ws_field[i];
    }
  }

  for (unsigned j = 0 ; j < meResult.size(); ++j) {
    EXPECT_NEAR(meResult[j], polyResult[j], tol);
  }
}
//-------------------------------------------------------------------------
void check_derivatives_at_ips(
  const stk::mesh::Entity* node_rels,
  const VectorFieldType& coordField,
  sierra::nalu::MasterElement& me)
{
  // Check that we can interpolate a random 3D linear field
  int dim = me.nDim_;

  std::mt19937 rng;
  rng.seed(0);
  std::uniform_real_distribution<double> coeff(-1.0, 1.0);
  std::vector<double> coeffs(dim);

  double a = coeff(rng);
  for (int j = 0; j < dim; ++j) {
    coeffs[j] = coeff(rng);
  }

  std::vector<double> polyResult(me.numIntPoints_ * dim);
  for (int j = 0; j < me.numIntPoints_; ++j) {
    for (int d = 0; d < dim; ++d) {
      polyResult[j*dim+d] = coeffs[d];
    }
  }

  std::vector<double> ws_field(me.nodesPerElement_);
  std::vector<double> ws_coords(me.nodesPerElement_ * dim);
  for (int j = 0; j < me.nodesPerElement_; ++j) {
    const double* coords = stk::mesh::field_data(coordField, node_rels[j]);
    for (int d = 0; d < dim; ++d) {
      ws_coords[j*dim+d] = coords[d];
    }
    ws_field[j] = linear_scalar_value(dim, a, coeffs.data(), coords);
  }

  std::vector<double> meResult(me.numIntPoints_ * dim, 0.0);
  std::vector<double> meGrad(me.numIntPoints_ * me.nodesPerElement_ * dim);
  std::vector<double> meDeriv(me.numIntPoints_ * me.nodesPerElement_ * dim);
  std::vector<double> meDetj(me.numIntPoints_);

  double error = 0.0;
  me.grad_op(1, ws_coords.data(), meGrad.data(), meDeriv.data(), meDetj.data(), &error);
  EXPECT_EQ(error, 0.0);

  for (int j = 0; j < me.numIntPoints_; ++j) {
    for (int i = 0; i < me.nodesPerElement_; ++i) {
      for (int d = 0; d < dim; ++d) {
        meResult[j*dim+d] += meGrad[j*me.nodesPerElement_*dim + i * dim + d] * ws_field[i];
      }
    }
  }

  // detj should be unity to floating point error
  for (unsigned j = 0 ; j < meDetj.size(); ++j) {
    EXPECT_NEAR(1, meDetj[j], tol) ;
  }

  // derivative should be exact to floating point error
  for (unsigned j = 0 ; j < meResult.size(); ++j) {
    EXPECT_NEAR(meResult[j], polyResult[j], tol);
  }
}
//-------------------------------------------------------------------------
void check_scv_shifted_ips_are_nodal(
  const stk::mesh::Entity* node_rels,
  const VectorFieldType& coordField,
  sierra::nalu::MasterElement& meSV)
{
  int dim = meSV.nDim_;
  std::vector<double> ws_coords(meSV.nodesPerElement_ * dim);

  for (int j = 0; j < meSV.nodesPerElement_; ++j) {
    const double* coords = stk::mesh::field_data(coordField, node_rels[j]);
    for (int d = 0; d < dim; ++d) {
      ws_coords[j*dim+d] = coords[d];
    }
  }

  const auto& shiftedIps = meSV.intgLocShift_;
  EXPECT_EQ(ws_coords.size(), shiftedIps.size()) << "P1 test";
  for (unsigned j = 0; j < shiftedIps.size(); ++j) {
    EXPECT_NEAR(ws_coords[j], shiftedIps[j], tol);
  }



}
//-------------------------------------------------------------------------
void check_is_in_element(
  const stk::mesh::Entity* node_rels,
  const VectorFieldType& coordField,
  sierra::nalu::MasterElement& me)
{
  int dim = me.nDim_;

  std::mt19937 rng;
  rng.seed(0);

  // randomly select a point within (boxmin, boxmax)^3 \subset reference element domain
  const double boxmin = 0.125;
  const double boxmax = 0.25;
  std::uniform_real_distribution<double> coeff(boxmin, boxmax);
  std::vector<double> random_pt(dim);
  for (int j = 0; j < dim; ++j) {
    random_pt[j] = coeff(rng);
  }

  // is in element uses a different stride for the coordinate data
  // compared to the gradient computation
  std::vector<double> ws_field(me.nodesPerElement_);
  std::vector<double> ws_coords(me.nodesPerElement_ * dim);

  // Hex8/Quad4's is_in_element and interpolatePoint routines use a different,
  // but self-consistent reference element compared to the core shape functions
  // and derivatives

  bool isHexSCS = dynamic_cast<sierra::nalu::HexSCS*>(&me) != nullptr;
  bool isQuadSCS = dynamic_cast<sierra::nalu::Quad2DSCS*>(&me) != nullptr;
  double fac = (isHexSCS || isQuadSCS) ? 2.0 : 1.0;

  for (int j = 0; j < me.nodesPerElement_; ++j) {
    const double* coords = stk::mesh::field_data(coordField, node_rels[j]);
    for (int d = 0; d < dim; ++d) {
      ws_coords[d * me.nodesPerElement_ + j] = fac*coords[d];
    }
  }

  std::vector<double> mePt(dim);
  auto dist = me.isInElement(ws_coords.data(), random_pt.data(), mePt.data());

  EXPECT_LT(dist, 1.0+tol);
  for (int d = 0; d < dim; ++d) {
    EXPECT_NEAR(random_pt[d], mePt[d], tol);
  }
}
//-------------------------------------------------------------------------
void check_is_not_in_element(
  const stk::mesh::Entity* node_rels,
  const VectorFieldType& coordField,
  sierra::nalu::MasterElement& me)
{
  int dim = me.nDim_;

  // choose a point not in the element
  std::vector<double> exterior_pt = { 100, 100, 100 };

  std::vector<double> ws_field(me.nodesPerElement_);
  std::vector<double> ws_coords(me.nodesPerElement_ * dim);

  for (int j = 0; j < me.nodesPerElement_; ++j) {
    const double* coords = stk::mesh::field_data(coordField, node_rels[j]);
    for (int d = 0; d < dim; ++d) {
      ws_coords[d * me.nodesPerElement_ + j] = coords[d];
    }
  }

  std::vector<double> mePt(dim);
  double dist = me.isInElement(ws_coords.data(), exterior_pt.data(), mePt.data());
  EXPECT_GT(dist, 1 + tol);
}
//-------------------------------------------------------------------------
void check_particle_interp(
  const stk::mesh::Entity* node_rels,
  const VectorFieldType& coordField,
  sierra::nalu::MasterElement& me)
{
  int dim = me.nDim_;

  std::mt19937 rng;
  rng.seed(0);

  // randomly select a point within (boxmin, boxmax)^3 \subset reference element domain
  const double boxmin = 0.125;
  const double boxmax = 0.25;
  std::uniform_real_distribution<double> coeff(boxmin, boxmax);
  std::vector<double> random_pt(dim);

  double const_value = coeff(rng);
  for (int j = 0; j < dim; ++j) {
    random_pt[j] = coeff(rng);
  }

  std::vector<double> coeffs(dim);
  for (int j = 0; j < dim; ++j) {
    coeffs[j] = coeff(rng);
  }

  // randomly perturb each of the coordinates of by a factor of delta
  // the element still needs to actually contain the box, (boxmin, boxmax)^3
  const double delta = 0.25;
  std::uniform_real_distribution<double> coord_perturb(-delta/2, delta/2);

  // is in element uses a different stride for the coordinate data
  // compared to the gradient computation
  std::vector<double> ws_field(me.nodesPerElement_);
  std::vector<double> ws_coords(me.nodesPerElement_ * dim);
  std::vector<double> perturbed_coords(dim);

  for (int j = 0; j < me.nodesPerElement_; ++j) {
    const double* coords = stk::mesh::field_data(coordField, node_rels[j]);
    for (int d = 0; d < dim; ++d) {
      perturbed_coords[d] = coords[d] + coord_perturb(rng);
      ws_coords[d * me.nodesPerElement_ + j] = perturbed_coords[d];
    }
    ws_field[j] = linear_scalar_value(dim, const_value, coeffs.data(), perturbed_coords.data());
  }

  std::vector<double> mePt(dim);
  double dist = me.isInElement(ws_coords.data(), random_pt.data(), mePt.data());
  EXPECT_LT(dist, 1.0+tol);

  double meInterp = 0.0;
  me.interpolatePoint(1, mePt.data(), ws_field.data(), &meInterp);
  double exactVal = linear_scalar_value(dim, const_value, coeffs.data(), random_pt.data());
  EXPECT_NEAR(meInterp, exactVal, tol);
}

}

class MasterElement : public ::testing::Test
{
protected:
    MasterElement() : comm(MPI_COMM_WORLD) {}

    void choose_topo(stk::topology topo)
    {
      meta = std::unique_ptr<stk::mesh::MetaData>(new stk::mesh::MetaData(topo.dimension()));
      bulk = std::unique_ptr<stk::mesh::BulkData>(new stk::mesh::BulkData(*meta, comm));
      elem = unit_test_utils::create_one_reference_element(*bulk, topo);
      meSS = unit_test_utils::get_surface_master_element(topo);
      meSV = unit_test_utils::get_volume_master_element(topo);
    }

    void scs_interpolation(stk::topology topo) {
      choose_topo(topo);
      check_interpolation_at_ips(bulk->begin_nodes(elem), coordinate_field(), *meSS);
    }

    void scv_interpolation(stk::topology topo) {
      choose_topo(topo);
      check_interpolation_at_ips(bulk->begin_nodes(elem), coordinate_field(), *meSV);
    }

    void scs_derivative(stk::topology topo) {
      choose_topo(topo);
      check_derivatives_at_ips(bulk->begin_nodes(elem), coordinate_field(), *meSS);
    }

    void is_not_in_element(stk::topology topo) {
      choose_topo(topo);
      check_is_not_in_element(bulk->begin_nodes(elem), coordinate_field(), *meSS);
    }

    void scv_shifted_ips_are_nodal(stk::topology topo) {
      choose_topo(topo);
      check_scv_shifted_ips_are_nodal(bulk->begin_nodes(elem), coordinate_field(), *meSV);
    }

    void is_in_element(stk::topology topo) {
      choose_topo(topo);
      check_is_in_element(bulk->begin_nodes(elem), coordinate_field(), *meSS);
    }

    void particle_interpolation(stk::topology topo) {
      choose_topo(topo);
      check_particle_interp(bulk->begin_nodes(elem), coordinate_field(), *meSS);
    }

    const VectorFieldType& coordinate_field() const {
      return *static_cast<const VectorFieldType*>(meta->coordinate_field());
    }

    stk::ParallelMachine comm;
    std::unique_ptr<stk::mesh::MetaData> meta;
    std::unique_ptr<stk::mesh::BulkData> bulk;
    stk::mesh::Entity elem;
    sierra::nalu::MasterElement* meSS;
    sierra::nalu::MasterElement* meSV;
};

#define TEST_F_ALL_TOPOS_NO_PYR(x, y) \
    TEST_F(x, tri##_##y)   { y(stk::topology::TRI_3_2D); }   \
    TEST_F(x, quad4##_##y)  { y(stk::topology::QUAD_4_2D); } \
    TEST_F(x, quad9##_##y)  { y(stk::topology::QUAD_9_2D); } \
    TEST_F(x, tet##_##y)   { y(stk::topology::TET_4); }      \
    TEST_F(x, wedge##_##y) { y(stk::topology::WEDGE_6); }    \
    TEST_F(x, hex8##_##y)   { y(stk::topology::HEX_8); }     \
    TEST_F(x, hex27##_##y)  { y(stk::topology::HEX_27); }

#define TEST_F_ALL_P1_TOPOS_NO_PYR(x, y) \
    TEST_F(x, tri##_##y)   { y(stk::topology::TRI_3_2D); }   \
    TEST_F(x, quad4##_##y)  { y(stk::topology::QUAD_4_2D); } \
    TEST_F(x, tet##_##y)   { y(stk::topology::TET_4); }      \
    TEST_F(x, wedge##_##y) { y(stk::topology::WEDGE_6); }    \
    TEST_F(x, hex8##_##y)   { y(stk::topology::HEX_8); }

// Patch tests: pyramids fail
TEST_F_ALL_TOPOS_NO_PYR(MasterElement, scs_interpolation);
TEST_F_ALL_TOPOS_NO_PYR(MasterElement, scs_derivative);
TEST_F_ALL_TOPOS_NO_PYR(MasterElement, scv_interpolation);

// Test needs to be modified substantially for HO.  Pyramid works.
TEST_F_ALL_P1_TOPOS_NO_PYR(MasterElement, scv_shifted_ips_are_nodal);
TEST_F(MasterElement, pyr_scv_shifted_ips_are_nodal) { scv_shifted_ips_are_nodal(stk::topology::PYRAMID_5); }

// HO elems / pyramids are missing these functions
TEST_F_ALL_P1_TOPOS_NO_PYR(MasterElement, is_not_in_element);
TEST_F_ALL_P1_TOPOS_NO_PYR(MasterElement, is_in_element);
TEST_F_ALL_P1_TOPOS_NO_PYR(MasterElement, particle_interpolation);
