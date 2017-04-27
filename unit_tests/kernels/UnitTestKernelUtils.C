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
} // unit_test_kernel_utils
