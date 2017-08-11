/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <MovingAveragePostProcessor.h>
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
ExponentialMovingAverager::ExponentialMovingAverager(double timeScale, bool isInit, double alpha)
: timeScale_(timeScale),
  isInit_(isInit),
  alpha_(alpha)
{ }
//--------------------------------------------------------------------------
void
ExponentialMovingAverager::compute_and_set_alpha(double delta_t)
{
  ThrowAssert(delta_t > 1.0e6 * std::numeric_limits<double>::min());
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
MovingAveragePostProcessor::MovingAveragePostProcessor(
  stk::mesh::BulkData& bulk,
  TimeIntegrator& timeIntegrator,
  bool isRestarted)
  :  bulk_(bulk),
     timeIntegrator_(timeIntegrator),
     isRestarted_(isRestarted)
{}
//--------------------------------------------------------------------------
void MovingAveragePostProcessor::add_fields(std::vector<std::string> fieldNames)
{
  for (const auto& fieldName : fieldNames) {
    auto& meta = bulk_.mesh_meta_data();
    auto* field = meta.get_field(stk::topology::NODE_RANK, fieldName);
    ThrowRequireMsg(field != nullptr, "Requested field `" + fieldName + "' not available for averaging");
    ThrowRequireMsg(field->type_is<double>(), "Only double precision-typed fields allowed");

    stk::mesh::FieldBase* avgField = meta.get_field(stk::topology::NODE_RANK, filtered_field_name(field->name()));
    ThrowRequireMsg(avgField != nullptr, filtered_field_name(field->name()) + " field not registered" );
    fieldMap_.insert({field, avgField});
  }
}
//--------------------------------------------------------------------------
void append_unique_parts(const stk::mesh::PartVector& partsToAdd, stk::mesh::PartVector& partVec)
{
  for (auto* part : partsToAdd) {
    ThrowRequire(part != nullptr);
    if (std::find(partVec.begin(), partVec.end(), part) == partVec.end()) {
      partVec.push_back(part);
    }
  }
}
//--------------------------------------------------------------------------
void
MovingAveragePostProcessor::add_parts_for_all_fields(stk::mesh::PartVector parts)
{
  for (const auto& fieldPair : fieldMap_) {
    append_unique_parts(parts,  partVecs_[fieldPair.first->name()]);
  }
}
//--------------------------------------------------------------------------
void
MovingAveragePostProcessor::set_time_scale(std::string fieldName, double timeScale)
{
  averagers_[fieldName] = ExponentialMovingAverager(timeScale, !isRestarted_);
}
//--------------------------------------------------------------------------
void
MovingAveragePostProcessor::set_time_scale(double timeScale)
{
  for (const auto& fieldPair : fieldMap_) {
    averagers_[fieldPair.first->name()]= ExponentialMovingAverager(timeScale, !isRestarted_);
  }
}
//--------------------------------------------------------------------------
void
MovingAveragePostProcessor::add_parts_for_field(std::string name, stk::mesh::PartVector parts)
{
  bool fieldIsInList = false;
  for (auto& fieldPair : fieldMap_) {
    if (fieldPair.first->name() == name) {
      fieldIsInList = true;
      auto& partVec = partVecs_[name];
      append_unique_parts(parts, partVec);
    }
  }
  ThrowRequireMsg(fieldIsInList, "Field `" +  name + "' not in list");
}
//--------------------------------------------------------------------------
void
MovingAveragePostProcessor::execute()
{
  double dt = timeIntegrator_.get_time_step();


  for (auto& fieldPair : fieldMap_) {
    const stk::mesh::FieldBase* field = fieldPair.first;
    const stk::mesh::FieldBase* avgField = fieldPair.second;
    auto& averager = averagers_.at(field->name());

    averager.compute_and_set_alpha(dt);

    // Average is performed on all point where the original field is defined,
    // unless a partvector is specified for the field, in which case that is used instead
    const auto& partVecForField = partVecs_[field->name()];
    stk::mesh::Selector selector = (partVecForField.empty()) ?
        stk::mesh::selectField(*field) : stk::mesh::selectUnion(partVecForField);

    const auto& node_buckets = bulk_.get_buckets(stk::topology::NODE_RANK, stk::mesh::selectField(*field));
    for (const auto* ib : node_buckets) {
      const auto& b = *ib;
      for (size_t k = 0u; k < b.size(); ++k) {
        const double* fieldVals = static_cast<double*>(stk::mesh::field_data(*field, b[k]));
        double* avgFieldVals = static_cast<double*>(stk::mesh::field_data(*avgField, b[k]));

        const unsigned fieldLength = field->max_size(stk::topology::NODE_RANK);
        ThrowAssert(fieldLength == avgField->max_size(stk::topology::NODE_RANK));

        for (unsigned d = 0; d < fieldLength; ++d) {
          avgFieldVals[d] = averager.compute_updated_average(avgFieldVals[d], fieldVals[d]);
        }
      }
    }
  }

  for (auto& averagerPair : averagers_) {
    averagerPair.second.init_state(false);
  }
}

} // namespace nalu
} // namespace Sierra
