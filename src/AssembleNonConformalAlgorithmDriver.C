/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <AssembleNonConformalAlgorithmDriver.h>
#include <Algorithm.h>
#include <AlgorithmDriver.h>
#include <FieldTypeDef.h>
#include <Realm.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

//==========================================================================
// Class Definition
//==========================================================================
// AssembleNonConformalAlgorithmDriver - Drives non-conformal algs
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleNonConformalAlgorithmDriver::AssembleNonConformalAlgorithmDriver(
  Realm &realm,
  stk::mesh::FieldBase *ncNormalFlux,
  stk::mesh::FieldBase *ncPenalty,
  ScalarFieldType *ncArea,
  const int fluxFieldSize)
  : AlgorithmDriver(realm),
    ncNormalFlux_(ncNormalFlux),
    ncPenalty_(ncPenalty),
    ncArea_(ncArea),
    fluxFieldSize_(fluxFieldSize)
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleNonConformalAlgorithmDriver::~AssembleNonConformalAlgorithmDriver()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- pre_work --------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNonConformalAlgorithmDriver::pre_work()
{
  // need to zero out the fields
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // define some common selectors; select all nodes (locally and shared)
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*ncPenalty_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * ncPenalty = (double*)stk::mesh::field_data(*ncPenalty_, b);
    double * ncNormalFlux = (double*)stk::mesh::field_data(*ncNormalFlux_, b);
    double * ncArea = (double*)stk::mesh::field_data(*ncArea_, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      ncArea[k] = 0.0;      
      ncPenalty[k] = 0.0;
      const int offSet = k*fluxFieldSize_;
      for ( unsigned i = 0; i < fluxFieldSize_; ++i ) {
        ncNormalFlux[offSet+i] = 0.0;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- post_work -------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNonConformalAlgorithmDriver::post_work()
{
  // assemble and normalize
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  std::vector<stk::mesh::FieldBase*> fields;
  fields.push_back(ncPenalty_);
  fields.push_back(ncNormalFlux_);
  fields.push_back(ncArea_);
  stk::mesh::parallel_sum(bulk_data, fields);

  if ( realm_.hasPeriodic_) {
    const int sizeOfField = 1;
    realm_.periodic_field_update(ncPenalty_, fluxFieldSize_);
    realm_.periodic_field_update(ncNormalFlux_, fluxFieldSize_);
    realm_.periodic_field_update(ncArea_, sizeOfField);
  }

  // normalize
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*ncPenalty_);
  
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * ncNormalFlux = (double*)stk::mesh::field_data(*ncNormalFlux_, b);
    double * ncPenalty = (double*)stk::mesh::field_data(*ncPenalty_, b);
    const double * ncArea = (double*)stk::mesh::field_data(*ncArea_, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const double theArea = ncArea[k];
      ncPenalty[k] /= theArea;
      const int offSet = k*fluxFieldSize_;
      for (unsigned i = 0; i < fluxFieldSize_; ++i) {
        ncNormalFlux[offSet+i] /= theArea;
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
