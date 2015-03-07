/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef FieldFunctions_h
#define FieldFunctions_h

#include <stk_mesh/base/MetaData.hpp>

namespace stk{
namespace mesh{
class FieldBase;
class BulkData;
}
}

namespace sierra{
namespace nalu{

// y = alpha*x + beta*y
void field_axpby(
  const stk::mesh::MetaData & metaData,
  const stk::mesh::BulkData & bulkData,
  const double alpha,
  const stk::mesh::FieldBase & xField,
  const double beta,
  const stk::mesh::FieldBase & yField,
  const stk::topology::rank_t entityRankValue=stk::topology::NODE_RANK);

// x = alpha
void field_fill(
  const stk::mesh::MetaData & metaData,
  const stk::mesh::BulkData & bulkData,
  const double alpha,
  const stk::mesh::FieldBase & xField,
  const stk::topology::rank_t entityRankValue=stk::topology::NODE_RANK);

// x = alpha * x
void field_scale(
  const stk::mesh::MetaData & metaData,
  const stk::mesh::BulkData & bulkData,
  const double alpha,
  const stk::mesh::FieldBase & xField,
  const stk::topology::rank_t entityRankValue=stk::topology::NODE_RANK);

// y = x
void field_copy(
  const stk::mesh::MetaData & metaData,
  const stk::mesh::BulkData & bulkData,
  const stk::mesh::FieldBase & xField,
  const stk::mesh::FieldBase & yField,
  const stk::topology::rank_t entityRankValue=stk::topology::NODE_RANK);

// y[compJ] = x[compK]
void field_index_copy(
  const stk::mesh::MetaData & metaData,
  const stk::mesh::BulkData & bulkData,
  const stk::mesh::FieldBase & xField,
  const int xFieldIndex,
  const stk::mesh::FieldBase & yField,
  const int yFieldIndex,
  const stk::topology::rank_t entityRankValue=stk::topology::NODE_RANK);

} // namespace nalu
} // namespace Sierra

#endif
