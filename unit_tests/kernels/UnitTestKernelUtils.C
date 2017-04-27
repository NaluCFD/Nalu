/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernels/UnitTestKernelUtils.h"
#include "UnitTestKokkosUtils.h"

#include <cmath>

namespace {

/** Trigonometric field functions for unit testing
 *
 * This convenience class generates trigonometric test functions for primitive
 * variables based on Steady Taylor Vortex and Steady Thermal MMS.
 */
struct TrigFieldFunction
{
  TrigFieldFunction()
    : pi(std::acos(-1.0))
  {}

  void velocity(const double* coords, double* qField) const
  {
    const double x = coords[0];
    const double y = coords[1];

    qField[0] = -unot * std::cos(a * pi * x) * std::sin(a * pi * y);
    qField[1] = +vnot * std::sin(a * pi * x) * std::cos(a * pi * y);
  }

  void pressure(const double* coords, double* qField) const
  {
    const double x = coords[0];
    const double y = coords[1];

    qField[0] = -pnot/4.0 * (
      std::cos(2.0 * a * pi * x) + std::cos(2.0 * a * pi * y));
  }

  void dpdx(const double* coords, double* qField) const
  {
    const double x = coords[0];
    const double y = coords[1];

    qField[0] = 0.5 * a * pi * std::sin(2.0 * a * pi * x);
    qField[1] = 0.5 * a * pi * std::sin(2.0 * a * pi * y);
  }

  void temperature(const double* coords, double* qField) const
  {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];

    qField[0] = (k/4.0) * (
      std::cos(2.0 * aT * pi * x) +
      std::cos(2.0 * aT * pi * y) +
      std::cos(2.0 * aT * pi * z));
  }

private:
  /// Factor for u-component of velocity
  static constexpr double unot{1.0};

  /// Factor for v-component of velocity
  static constexpr double vnot{1.0};

  /// Factor for  pressure field
  static constexpr double pnot{1.0};

  /// Factor for temperature field
  static constexpr double k{1.0};

  /// Frequency multiplier for velocity and pressure fields
  static constexpr double a{20.0};

  /// Frequency for temperature fields
  static constexpr double aT{1.0};

  const double pi;
};

/** Initialize the field array with the trigonometric test function
 *
 * @param bulk Reference to STK BulkData
 * @param qField Field to be initialized. Template specialized on this field
 */
template<typename T>
void init_trigonometric_field(
  const stk::mesh::BulkData& bulk,
  const VectorFieldType& coordinates,
  T& qField)
{
  using FieldInitFunction = void (TrigFieldFunction::*)(const double*, double*) const;
  const TrigFieldFunction stv;
  const auto fieldName = qField.name();
  FieldInitFunction funcPtr = nullptr;

  if (fieldName == "velocity")
    funcPtr = &TrigFieldFunction::velocity;
  else if (fieldName == "pressure")
    funcPtr = &TrigFieldFunction::pressure;
  else if (fieldName == "dpdx")
    funcPtr = &TrigFieldFunction::dpdx;
  else if (fieldName == "temperature")
    funcPtr = &TrigFieldFunction::temperature;
  else
    funcPtr = nullptr;

  EXPECT_TRUE(funcPtr != nullptr);

  const auto& meta = bulk.mesh_meta_data();
  EXPECT_EQ(meta.spatial_dimension(), 3u);

  const stk::mesh::Selector selector =
    meta.locally_owned_part() | meta.globally_shared_part();
  const auto& buckets = bulk.get_buckets(stk::topology::NODE_RANK, selector);

  kokkos_thread_team_bucket_loop(buckets, [&](stk::mesh::Entity node) {
      const double* coords = stk::mesh::field_data(coordinates, node);
      double* qNode = stk::mesh::field_data(qField, node);

      ((stv).*(funcPtr))(coords, qNode);
    });
}

} // anonymous namespace

namespace unit_test_kernel_utils {

void velocity_test_function(
  const stk::mesh::BulkData& bulk,
  const VectorFieldType& coordinates,
  VectorFieldType& velocity)
{
  // Add additional test functions in future?
  init_trigonometric_field(bulk, coordinates, velocity);
}

void pressure_test_function(
  const stk::mesh::BulkData& bulk,
  const VectorFieldType& coordinates,
  ScalarFieldType& pressure)
{
  init_trigonometric_field(bulk, coordinates, pressure);
}

void dpdx_test_function(
  const stk::mesh::BulkData& bulk,
  const VectorFieldType& coordinates,
  VectorFieldType& dpdx)
{
  init_trigonometric_field(bulk, coordinates, dpdx);
}

void temperature_test_function(
  const stk::mesh::BulkData& bulk,
  const VectorFieldType& coordinates,
  ScalarFieldType& temperature)
{
  init_trigonometric_field(bulk, coordinates, temperature);
}


void calc_mass_flow_rate_scs(
  const stk::mesh::BulkData& bulk,
  const stk::topology& topo,
  const VectorFieldType& coordinates,
  const ScalarFieldType& density,
  const VectorFieldType& velocity,
  const GenericFieldType& massFlowRate)
{
  const auto& meta = bulk.mesh_meta_data();
  const int ndim = meta.spatial_dimension();
  EXPECT_EQ(ndim, 3);

  sierra::nalu::ElemDataRequests dataNeeded;

  const ScalarFieldType& densityNp1 = density.field_of_state(stk::mesh::StateNP1);
  const VectorFieldType& velocityNp1 = velocity.field_of_state(stk::mesh::StateNP1);
  auto meSCS = sierra::nalu::get_surface_master_element(topo);

  dataNeeded.add_cvfem_surface_me(meSCS);
  dataNeeded.add_gathered_nodal_field(coordinates, ndim);
  dataNeeded.add_gathered_nodal_field(densityNp1, 1);
  dataNeeded.add_gathered_nodal_field(velocityNp1, ndim);
  dataNeeded.add_master_element_call(sierra::nalu::SCS_AREAV);

  const stk::mesh::Selector selector =
    meta.locally_owned_part() | meta.globally_shared_part();
  const auto& buckets = bulk.get_buckets(stk::topology::ELEM_RANK, selector);

  const int bytes_per_team = 0;
  const int bytes_per_thread = sierra::nalu::get_num_bytes_pre_req_data(
    dataNeeded, meta.spatial_dimension()) ;

  auto team_exec = sierra::nalu::get_team_policy(
    buckets.size(), bytes_per_team, bytes_per_thread);

  Kokkos::parallel_for(team_exec, [&](const sierra::nalu::TeamHandleType& team) {
      auto& b = *buckets[team.league_rank()];
      const auto length = b.size();

      EXPECT_EQ(b.topology(), topo);

      sierra::nalu::ScratchViews preReqData(
        team, bulk, topo, dataNeeded);

      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, length), [&](const size_t& k) {
          stk::mesh::Entity element = b[k];
          sierra::nalu::fill_pre_req_data(
            dataNeeded, bulk, topo, element, &coordinates, preReqData);

          std::vector<double> rhoU(ndim);
          std::vector<double> v_shape_function(
            meSCS->numIntPoints_ * meSCS->nodesPerElement_);
          meSCS->shape_fcn(v_shape_function.data());
          double *mdot = stk::mesh::field_data(massFlowRate, element);
          auto& v_densityNp1 = preReqData.get_scratch_view_1D(densityNp1);
          auto& v_velocityNp1 = preReqData.get_scratch_view_2D(velocityNp1);
          auto& v_scs_areav = preReqData.scs_areav;

          for (int ip=0; ip < meSCS->numIntPoints_; ++ip) {
            for (int j=0; j < ndim; ++j) {
              rhoU[j] = 0.0;
            }
            const int offset = ip * meSCS->nodesPerElement_;

            for (int ic=0; ic < meSCS->nodesPerElement_; ++ic) {
              const double r = v_shape_function[offset+ic];
              for (int j=0; j < ndim; ++j) {
                rhoU[j] += r * v_densityNp1(ic) * v_velocityNp1(ic, j);
              }
            }

            double tmdot = 0.0;
            for (int j=0; j < ndim; j++) {
              tmdot += rhoU[j] * v_scs_areav(ip, j);
            }
            mdot[ip] = tmdot;
          }
        });
    });
}

void expect_all_near(
  const sierra::nalu::SharedMemView<double*>& calcValue,
  const double* exactValue,
  const double tol)
{
  const int length = calcValue.dimension(0);

  for (int i=0; i < length; ++i) {
    EXPECT_NEAR(calcValue[i], exactValue[i], tol);
  }
}

void expect_all_near(
  const sierra::nalu::SharedMemView<double*>& calcValue,
  const double exactValue,
  const double tol)
{
  const int length = calcValue.dimension(0);

  for (int i=0; i < length; ++i) {
    EXPECT_NEAR(calcValue[i], exactValue, tol);
  }
}

void expect_all_near(
  const sierra::nalu::SharedMemView<double**>& calcValue,
  const double* exactValue,
  const double tol)
{
  const int dim1 = calcValue.dimension(0);
  const int dim2 = calcValue.dimension(1);

  for (int i=0; i < dim1; i++)
    for (int j=0; j < dim2; j++)
      EXPECT_NEAR(calcValue(i,j),exactValue[i*dim2+j], tol);
}

} // unit_test_kernel_utils
