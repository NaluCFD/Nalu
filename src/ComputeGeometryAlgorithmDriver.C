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
#include <SolutionOptions.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

void compute_open_mean_normal(const stk::mesh::BulkData& bulk)
{
  const auto& meta = bulk.mesh_meta_data();
  GenericFieldType& areavField = *meta.get_field<GenericFieldType>(meta.side_rank(), "exposed_area_vector");

  ThrowRequireMsg(meta.get_field<VectorFieldType>(stk::topology::NODE_RANK, "average_open_normal") != nullptr,
    "average_open_normal field required");

  VectorFieldType& meanNormalField = *meta.get_field<VectorFieldType>(stk::topology::NODE_RANK, "average_open_normal");

  ThrowRequireMsg(meta.get_field<ScalarIntFieldType>(stk::topology::NODE_RANK, "open_face_connection_count") != nullptr,
    "open_face_connection_count field required");
  ScalarIntFieldType& faceConnectionCount = *meta.get_field<ScalarIntFieldType>(stk::topology::NODE_RANK, "open_face_connection_count");

  const stk::mesh::Selector& locally_owned_open = meta.locally_owned_part() & stk::mesh::selectField(meanNormalField);
  const stk::mesh::BucketVector& face_buckets = bulk.get_buckets(meta.side_rank(), locally_owned_open);

  const int nDim = meta.spatial_dimension();

  for (const auto* ib : face_buckets) {
    const auto& b = *ib;
    const size_t length = b.size();

    MasterElement& meFC = *MasterElementRepo::get_surface_master_element(b.topology());
    const int numIp = meFC.numIntPoints_;
    const int* ipNodeMap = meFC.ipNodeMap();

    for (size_t k = 0u; k < length; ++k) {
      const double* areaVecs = stk::mesh::field_data(areavField, b, k);
      const auto* face_node_rels = bulk.begin_nodes(b[k]);

      for (int ip = 0; ip < numIp; ++ip) {
        const double* const areav = &areaVecs[ip*nDim];

        double aMag = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          aMag += areav[j] * areav[j];
        }
        ThrowAssert(aMag > std::numeric_limits<double>::min());
        const double inv_aMag = 1.0 / std::sqrt(aMag);

        const stk::mesh::Entity nearestNode = face_node_rels[ipNodeMap[ip]];
        double* meanNormal = stk::mesh::field_data(meanNormalField, nearestNode);
        ThrowRequire(meanNormal != nullptr);

        for (int j = 0; j < nDim; ++j) {
          meanNormal[j] += areav[j] * inv_aMag;
        }
        ThrowRequire(stk::mesh::field_data(faceConnectionCount, nearestNode) != nullptr);
        *stk::mesh::field_data(faceConnectionCount, nearestNode) += 1;
      }
    }
  }

  stk::mesh::parallel_sum(bulk, {&meanNormalField});
  stk::mesh::parallel_sum(bulk, {&faceConnectionCount});

  const stk::mesh::Selector& owned_or_shared_open = (meta.locally_owned_part() | meta.globally_shared_part()) & stk::mesh::selectField(meanNormalField);
  const stk::mesh::BucketVector& node_buckets = bulk.get_buckets(stk::topology::NODE_RANK, owned_or_shared_open);
  for (const auto* ib : node_buckets) {
    const auto& b = *ib;
    const size_t length = b.size();
    for (size_t k = 0u; k < length; ++k) {
      double* meanNormal = stk::mesh::field_data(meanNormalField, b, k);
      const double avgFactor = *stk::mesh::field_data(faceConnectionCount, b, k);
      ThrowAssert(avgFactor >= 1);
      for (int j = 0; j < nDim; ++j) {
        meanNormal[j] /= avgFactor;
      }
    }
  }
}

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

  if (realm_.solutionOptions_->explicitlyZeroOpenPressureGradient_) {
    ThrowRequireMsg(meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "average_open_normal") != nullptr,
      "average_open_normal field required");

    VectorFieldType& meanNormalField = *meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "average_open_normal");
    stk::mesh::field_fill(0.0, meanNormalField );

    ThrowRequireMsg(meta_data.get_field<ScalarIntFieldType>(stk::topology::NODE_RANK, "open_face_connection_count") != nullptr,
      "open_face_connection_count field required");
    ScalarIntFieldType& faceConnectionCount = *meta_data.get_field<ScalarIntFieldType>(stk::topology::NODE_RANK, "open_face_connection_count");
    stk::mesh::field_fill(0, faceConnectionCount );
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

  if (realm_.solutionOptions_->explicitlyZeroOpenPressureGradient_) {
    compute_open_mean_normal(realm_.bulk_data());
  }
}

} // namespace nalu
} // namespace Sierra
