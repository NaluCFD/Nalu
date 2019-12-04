/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <WallFunctionParamsAlgorithmDriver.h>
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
// WallFunctionParamsAlgorithmDriver - Drives nodal grad algorithms
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
WallFunctionParamsAlgorithmDriver::WallFunctionParamsAlgorithmDriver(
  Realm &realm)
  : AlgorithmDriver(realm),
    assembledWallArea_(NULL),
    assembledWallNormalDistance_(NULL)
{
  // register the fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  assembledWallArea_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_wall_area_wf");
  assembledWallNormalDistance_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_wall_normal_distance");
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
WallFunctionParamsAlgorithmDriver::~WallFunctionParamsAlgorithmDriver()
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- pre_work --------------------------------------------------------
//--------------------------------------------------------------------------
void
WallFunctionParamsAlgorithmDriver::pre_work()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // define some common selectors; select all nodes (locally and shared)
  // where assembledWallArea is defined
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*assembledWallArea_);

  //===========================================================
  // zero out nodal fields
  //===========================================================

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    const stk::mesh::Bucket::size_type length   = b.size();
    double * assembledWallArea = stk::mesh::field_data(*assembledWallArea_, b);
    double * assembledWallNormalDistance = stk::mesh::field_data(*assembledWallNormalDistance_, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      assembledWallArea[k] = 0.0;
      assembledWallNormalDistance[k] = 0.0;
    }
  }
}

//--------------------------------------------------------------------------
//-------- post_work -------------------------------------------------------
//--------------------------------------------------------------------------
void
WallFunctionParamsAlgorithmDriver::post_work()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // parallel assemble
  stk::mesh::parallel_sum(bulk_data, {assembledWallArea_, assembledWallNormalDistance_});
  
  // periodic assemble
  if ( realm_.hasPeriodic_) {
    const unsigned fieldSize = 1;
    const bool bypassFieldCheck = false; // fields are not defined at all slave/master node pairs
    realm_.periodic_field_update(assembledWallArea_, fieldSize, bypassFieldCheck);
    realm_.periodic_field_update(assembledWallNormalDistance_, fieldSize, bypassFieldCheck);
  }

  // normalize
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*assembledWallArea_);
  
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    const stk::mesh::Bucket::size_type length   = b.size();
    const double * assembledWallArea = stk::mesh::field_data(*assembledWallArea_, b);
    double * assembledWallNormalDistance = stk::mesh::field_data(*assembledWallNormalDistance_, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      assembledWallNormalDistance[k] /= assembledWallArea[k];
    }
  }
}

} // namespace nalu
} // namespace Sierra
