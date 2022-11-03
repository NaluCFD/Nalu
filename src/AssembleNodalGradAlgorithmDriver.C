/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <AssembleNodalGradAlgorithmDriver.h>
#include <Algorithm.h>
#include <AlgorithmDriver.h>
#include <FieldTypeDef.h>
#include <FieldFunctions.h>
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
// AssembleNodalGradAlgorithmDriver - Drives nodal grad algorithms
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleNodalGradAlgorithmDriver::AssembleNodalGradAlgorithmDriver(
  Realm &realm,
  const std::string & scalarQName,
  const std::string & dqdxName,
  const std::string & areaWeightName,
  const bool areaWeight)
  : AlgorithmDriver(realm),
    scalarQName_(scalarQName),
    dqdxName_(dqdxName),
    areaWeightName_(areaWeightName),
    areaWeight_(areaWeight)
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
AssembleNodalGradAlgorithmDriver::~AssembleNodalGradAlgorithmDriver()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- pre_work --------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNodalGradAlgorithmDriver::pre_work()
{
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  const int nDim = metaData.spatial_dimension();

  // extract fields
  VectorFieldType *dqdx = metaData.get_field<double>(stk::topology::NODE_RANK, dqdxName_);

  // define some common selectors; select all nodes (locally and shared)
  // where dqdx is defined
  stk::mesh::Selector s_all_nodes
    = (metaData.locally_owned_part() | metaData.globally_shared_part())
    &stk::mesh::selectField(*dqdx);

  //===========================================================
  // zero out nodal gradient
  //===========================================================

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    const stk::mesh::Bucket::size_type length   = b.size();
    double * gq = stk::mesh::field_data(*dqdx, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int offSet = k*nDim;

      for ( int j = 0; j < nDim; ++j ) {
        gq[offSet+j] = 0.0;
      }
    }
  }

  if ( areaWeight_ ) {
    VectorFieldType *areaWeight = metaData.get_field<double>(stk::topology::NODE_RANK, areaWeightName_);
    field_fill( metaData, bulkData, 0.0, *areaWeight, realm_.get_activate_aura());
  }
}

//--------------------------------------------------------------------------
//-------- post_work -------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNodalGradAlgorithmDriver::post_work()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const unsigned nDim = meta_data.spatial_dimension();

  // extract fields
  VectorFieldType *dqdx = meta_data.get_field<double>(stk::topology::NODE_RANK, dqdxName_);
  std::vector<const stk::mesh::FieldBase*> sum_fields(1, dqdx);
  stk::mesh::parallel_sum(bulk_data, sum_fields);

  if ( realm_.hasPeriodic_) {
    realm_.periodic_field_update(dqdx, nDim);
  }

  // allow for area weighting
  if ( areaWeight_ ) {
    VectorFieldType *areaWeight = meta_data.get_field<double>(stk::topology::NODE_RANK, areaWeightName_);
    std::vector<const stk::mesh::FieldBase*> sum_area(1, areaWeight);
    stk::mesh::parallel_sum(bulk_data, sum_area);
    if ( realm_.hasPeriodic_) {
      realm_.periodic_field_update(areaWeight, nDim);
    }
    normalize_by_area();
  }

  // assemble the projected nodal gradient
  if ( realm_.hasOverset_ ) {
    realm_.overset_constraint_node_field_update(dqdx, 1, nDim);
  }
}

//--------------------------------------------------------------------------
//-------- normalize_by_area -----------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNodalGradAlgorithmDriver::normalize_by_area()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // extract fields
  VectorFieldType *dqdx = meta_data.get_field<double>(stk::topology::NODE_RANK, dqdxName_);
  VectorFieldType *areaWeight = meta_data.get_field<double>(stk::topology::NODE_RANK, areaWeightName_);

  // define some common selectors; select all nodes (locally and shared)
  // where dqdx is defined
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*dqdx);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    const stk::mesh::Bucket::size_type length   = b.size();
    double * gq = stk::mesh::field_data(*dqdx, b);
    const double * aw = stk::mesh::field_data(*areaWeight, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int offSet = k*nDim;

      for ( int j = 0; j < nDim; ++j ) {
        gq[offSet+j] /= aw[offSet+j];
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
