/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ComputeGeometryAlgorithmDriver.h>

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
// ComputeGeometryAlgorithmDriver - Drives nodal grad algorithms
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ComputeGeometryAlgorithmDriver::ComputeGeometryAlgorithmDriver(
  Realm &realm)
  : AlgorithmDriver(realm)
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- pre_work --------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeGeometryAlgorithmDriver::pre_work()
{
  // meta and bulk data
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // extract field that is always germane
  ScalarFieldType *dualNodalVolume = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  // define some common selectors
  stk::mesh::Selector s_locally_owned = meta_data.locally_owned_part();

  //====================================================
  // Initialize nodal volume and area vector to zero
  //====================================================

  // nodal fields first
  stk::mesh::Selector s_all_vol
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*dualNodalVolume);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_vol );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * nv = stk::mesh::field_data( *dualNodalVolume, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      nv[k] = 0.0;
    }
  }

  // handle case for realm using edge-based
  if ( realm_.realmUsesEdges_ ) {

    VectorFieldType *edgeAreaVec = meta_data.get_field<VectorFieldType>(stk::topology::EDGE_RANK, "edge_area_vector");

    // edge fields second
    stk::mesh::Selector s_all_area
      = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
      &stk::mesh::selectField(*edgeAreaVec);
    stk::mesh::BucketVector const& edge_buckets =
      realm_.get_buckets( stk::topology::EDGE_RANK, s_all_area );
    for ( stk::mesh::BucketVector::const_iterator ib = edge_buckets.begin();
          ib != edge_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      const stk::mesh::Bucket::size_type length   = b.size();
      double * av = stk::mesh::field_data(*edgeAreaVec, b);
      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        const int offSet = k*nDim;
        for ( int j = 0; j < nDim; ++j )
          av[offSet+j] = 0.0;
      }
    }
  }

}

//--------------------------------------------------------------------------
//-------- post_work -------------------------------------------------------
//--------------------------------------------------------------------------
void
ComputeGeometryAlgorithmDriver::post_work()
{

  // meta and bulk data
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  std::vector<stk::mesh::FieldBase*> sum_fields;

  // extract field always germane
  ScalarFieldType *dualNodalVolume = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  sum_fields.push_back(dualNodalVolume);

  // handle case for realm using edge-based
  if ( realm_.realmUsesEdges_ ) {
    VectorFieldType *edgeAreaVec = meta_data.get_field<VectorFieldType>(stk::topology::EDGE_RANK, "edge_area_vector");
    sum_fields.push_back(edgeAreaVec);
  }

  // deal with parallel
  stk::mesh::parallel_sum(bulk_data, sum_fields);

  if ( realm_.hasPeriodic_) {
    const unsigned fieldSize = 1;
    realm_.periodic_field_update(dualNodalVolume, fieldSize);
  }

}

} // namespace nalu
} // namespace Sierra
