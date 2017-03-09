#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <random>
#include <tuple>
#include <ostream>
#include <memory>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <master_element/MasterElement.h>
#include <NaluEnv.h>

#include "UnitTestUtils.h"

namespace {

void matvec22(const double* A, const double* x, double* b)
{
  b[0] = A[0] * x[0] + A[1] * x[1];
  b[1] = A[2] * x[0] + A[3] * x[1];
}

void mxm22(const double* A, const double* B, double* C)
{
  C[0] = A[0] * B[0] + A[1] * B[2];
  C[1] = A[0] * B[1] + A[1] * B[3];
  C[2] = A[2] * B[0] + A[3] * B[2];
  C[3] = A[2] * B[1] + A[3] * B[3];
}

void mxm33(const double* A, const double* B, double* C)
{
  enum { XX = 0, XY = 1, XZ = 2, YX = 3, YY = 4, YZ = 5, ZX = 6, ZY = 7, ZZ = 8 };
  C[XX] = A[XX] * B[XX] + A[XY] * B[YX] + A[XZ] * B[ZX];
  C[XY] = A[XX] * B[XY] + A[XY] * B[YY] + A[XZ] * B[ZY];
  C[XZ] = A[XX] * B[XZ] + A[XY] * B[YZ] + A[XZ] * B[ZZ];
  C[YX] = A[YX] * B[XX] + A[YY] * B[YX] + A[YZ] * B[ZX];
  C[YY] = A[YX] * B[XY] + A[YY] * B[YY] + A[YZ] * B[ZY];
  C[YZ] = A[YX] * B[XZ] + A[YY] * B[YZ] + A[YZ] * B[ZZ];
  C[ZX] = A[ZX] * B[XX] + A[ZY] * B[YX] + A[ZZ] * B[ZX];
  C[ZY] = A[ZX] * B[XY] + A[ZY] * B[YY] + A[ZZ] * B[ZY];
  C[ZZ] = A[ZX] * B[XZ] + A[ZY] * B[YZ] + A[ZZ] * B[ZZ];
}


void matvec33(const double* A, const double* x, double* b)
{
  b[0] = A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
  b[1] = A[3] * x[0] + A[4] * x[1] + A[5] * x[2];
  b[2] = A[6] * x[0] + A[7] * x[1] + A[8] * x[2];
}

void transpose22(const double* A, double* At)
{
  At[0] = A[0];
  At[1] = A[2];
  At[2] = A[1];
  At[3] = A[3];
}

void transpose33(const double* A, double* At)
{
  enum { XX = 0, XY = 1, XZ = 2, YX = 3, YY = 4, YZ = 5, ZX = 6, ZY = 7, ZZ = 8 };
  At[XX] = A[XX];
  At[YY] = A[YY];
  At[ZZ] = A[ZZ];

  At[XY] = A[YX];
  At[XZ] = A[ZX];
  At[YX] = A[XY];
  At[YZ] = A[ZY];
  At[ZX] = A[XZ];
  At[ZY] = A[YZ];
}

std::pair<std::vector<double>, std::vector<double>>
calculate_metric_tensor(sierra::nalu::MasterElement& me, const std::vector<double>& ws_coords)
{
  double scs_error = 0.0;
  int gradSize = me.numIntPoints_ * me.nodesPerElement_ * me.nDim_;
  std::vector<double> ws_dndx(gradSize);
  std::vector<double> ws_deriv(gradSize);
  std::vector<double> ws_det_j(me.numIntPoints_);
  me.grad_op(1, ws_coords.data(), ws_dndx.data(), ws_deriv.data(), ws_det_j.data(), &scs_error);

  int metricSize = me.nDim_ * me.nDim_ * me.numIntPoints_;
  std::vector<double> contravariant_metric_tensor(metricSize);
  std::vector<double> covariant_metric_tensor(metricSize);
  me.gij(ws_coords.data(), contravariant_metric_tensor.data(), covariant_metric_tensor.data(), ws_deriv.data());

  return {contravariant_metric_tensor, covariant_metric_tensor};
}

using VectorFieldType = stk::mesh::Field<double, stk::mesh::Cartesian>;

void test_metric_for_topo_2D(stk::topology topo, double tol) {
  int dim = topo.dimension();
  ASSERT_EQ(dim, 2);

  stk::mesh::MetaData meta(dim);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  stk::mesh::Entity elem = unit_test_utils::create_one_reference_element(bulk, topo);

  auto* mescs = unit_test_utils::get_surface_master_element(topo);

  // apply some arbitrary linear map the reference element
  std::mt19937 rng;
  rng.seed(0); // fixed seed
  std::uniform_real_distribution<double> coeff(-1.0, 1.0);

  double Q[4] = { 1.0+std::abs(coeff(rng)), coeff(rng), coeff(rng), 1.0+std::abs(coeff(rng)) };

  double Qt[4];
  transpose22(Q, Qt);

  double metric_exact[4];
  mxm22(Q,Qt,metric_exact);


  const auto& coordField = *static_cast<const VectorFieldType*>(meta.coordinate_field());
  std::vector<double> ws_coords(topo.num_nodes() * dim);
  const auto* nodes = bulk.begin_nodes(elem);
  for (unsigned j = 0; j < topo.num_nodes(); ++j) {
    const double* coords = stk::mesh::field_data(coordField, nodes[j]);
    matvec22(Q, coords, &ws_coords[j*dim]);
  }

  std::vector<double> contravariant_metric; std::vector<double> covariant_metric;
  std::tie(contravariant_metric, covariant_metric) = calculate_metric_tensor(*mescs, ws_coords);

  for (int ip = 0; ip < mescs->numIntPoints_; ++ip) {
    double identity[4] = {1.0,0.0,0.0,1.0};
    double shouldBeIdentity[4];
    mxm22(&contravariant_metric[4*ip], &covariant_metric[4*ip], shouldBeIdentity);
    for (unsigned k = 0; k < 4; ++k) {
      EXPECT_NEAR(contravariant_metric[4*ip+k], metric_exact[k], tol);
      EXPECT_NEAR(shouldBeIdentity[k], identity[k], tol);
    }
  }
}

void test_metric_for_topo_3D(stk::topology topo, double tol) {
  int dim = topo.dimension();
  ASSERT_EQ(dim,3);

  stk::mesh::MetaData meta(dim);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);
  stk::mesh::Entity elem = unit_test_utils::create_one_reference_element(bulk, topo);

  auto* mescs = unit_test_utils::get_surface_master_element(topo);

  // apply some arbitrary linear map the reference element
  std::mt19937 rng;
  rng.seed(0); // fixed seed
  std::uniform_real_distribution<double> coeff(-1.0, 1.0);

  double Q[9] = {
      1.0+std::abs(coeff(rng)), coeff(rng), coeff(rng),
      coeff(rng), 1.0+std::abs(coeff(rng)), coeff(rng),
      coeff(rng), coeff(rng), 1.0+std::abs(coeff(rng))
  };

  double Qt[9];
  transpose33(Q,Qt);

  double metric_exact[9];
  mxm33(Q,Qt,metric_exact);


  const auto& coordField = *static_cast<const VectorFieldType*>(meta.coordinate_field());

  std::vector<double> ws_coords(topo.num_nodes() * dim);
  const auto* nodes = bulk.begin_nodes(elem);
  for (unsigned j = 0; j < topo.num_nodes(); ++j) {
    const double* coords = stk::mesh::field_data(coordField, nodes[j]);
    matvec33(Q, coords, &ws_coords[j*dim]);
  }

  std::vector<double> contravariant_metric; std::vector<double> covariant_metric;
  std::tie(contravariant_metric, covariant_metric) = calculate_metric_tensor(*mescs, ws_coords);

  for (int ip = 0; ip < mescs->numIntPoints_; ++ip) {
    double identity[9] = {
        1.0,0.0,0.0,
        0.0,1.0,0.0,
        0.0,0.0,1.0
    };

    double shouldBeIdentity[9];
    mxm33(&contravariant_metric[9*ip], &covariant_metric[9*ip], shouldBeIdentity);
    for (unsigned k = 0; k < 9; ++k) {
      EXPECT_NEAR(contravariant_metric[9*ip+k], metric_exact[k], tol);
      EXPECT_NEAR(shouldBeIdentity[k], identity[k], tol);
    }
  }
}

}

TEST(MetricTensor, tri3)
{
  test_metric_for_topo_2D(stk::topology::TRIANGLE_3_2D, 1.0e-10);
}

TEST(MetricTensor, quad4)
{
  test_metric_for_topo_2D(stk::topology::QUADRILATERAL_4_2D,1.0e-10);
}

TEST(MetricTensor, quad9)
{
  test_metric_for_topo_2D(stk::topology::QUADRILATERAL_9_2D,1.0e-10);
}

TEST(MetricTensor, tet4)
{
  test_metric_for_topo_3D(stk::topology::TET_4,1.0e-10);
}

TEST(MetricTensor, wedge6)
{
  test_metric_for_topo_3D(stk::topology::WEDGE_6,1.0e-10);
}

TEST(MetricTensor, hex8)
{
  test_metric_for_topo_3D(stk::topology::HEX_8,1.0e-10);
}
TEST(MetricTensor, hex27)
{
  test_metric_for_topo_3D(stk::topology::HEX_27,1.0e-10);
}
