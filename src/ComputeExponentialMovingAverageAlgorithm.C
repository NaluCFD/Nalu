/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <ComputeExponentialMovingAverageAlgorithm.h>
#include <Algorithm.h>

#include <FieldTypeDef.h>
#include <Realm.h>
#include <TimeIntegrator.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

#include <limits>

namespace sierra{
namespace nalu{

//--------------------------------------------------------------------------
ExponentialMovingAverager::ExponentialMovingAverager(double timeScale)
: timeScale_(timeScale),
  isInit_(true),
  alpha_(-1)
{
  ThrowRequireMsg(timeScale_ > 1.0e6 * std::numeric_limits<double>::min(),
    "Time scale must be strictly positive");
}
//--------------------------------------------------------------------------
void
ExponentialMovingAverager::compute_and_set_alpha(double delta_t)
{
  ThrowAssert(delta_t >1.0e6 * std::numeric_limits<double>::min());
  alpha_ = (isInit_) ? 1.0 : std::min(1.0, delta_t / timeScale_);
}
//--------------------------------------------------------------------------
double
ExponentialMovingAverager::compute_updated_average(double oldAvg, double newVal)
{
  return (alpha_ * newVal  + (1 -  alpha_) * oldAvg);
}
//--------------------------------------------------------------------------
void
ExponentialMovingAverager::init_state(bool init)
{
  isInit_ = init;
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
ComputeExponentialMovingAverageAlgorithm::ComputeExponentialMovingAverageAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  stk::mesh::FieldBase* field,
  double timeScale)
  : Algorithm(realm, part),
    field_(field),
    averager_(ExponentialMovingAverager(timeScale))
{
  stk::mesh::MetaData& meta_data = realm_.meta_data();
  ThrowRequireMsg(field_ != nullptr, "Field argument invalid");

  ThrowRequireMsg(field_->type_is<double>(), "Only double precision-typed fields allowed");

  avgField_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "time_filtered_" + field->name());
  ThrowRequireMsg(avgField_ != nullptr, "time_filtered_" + field->name() + " field not registered");
}
//--------------------------------------------------------------------------
void
ComputeExponentialMovingAverageAlgorithm::execute()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  stk::mesh::Selector s_all_nodes = stk::mesh::selectUnion(partVec_);

  double dt  = realm_.timeIntegrator_->get_time_step();
  averager_.compute_and_set_alpha(dt);

  unsigned fieldLength = field_->max_size(stk::topology::NODE_RANK);
  ThrowRequire(fieldLength = avgField_->max_size(stk::topology::NODE_RANK));

  // initialize to something largely negative
  const auto& node_buckets = bulk_data.get_buckets(stk::topology::NODE_RANK, s_all_nodes);
  for (const auto* ib : node_buckets) {
    const auto& b = *ib;
    for (size_t k = 0u; k < b.size(); ++k) {
      const double* fieldVals = static_cast<double*>(stk::mesh::field_data(*field_, b[k]));
      double* avgFieldVals = static_cast<double*>(stk::mesh::field_data(*avgField_, b[k]));
      for (unsigned d = 0; d < fieldLength; ++d) {
        avgFieldVals[d] = averager_.compute_updated_average(avgFieldVals[d], fieldVals[d]);
      }
    }
  }
  averager_.init_state(false);
}

} // namespace nalu
} // namespace Sierra
