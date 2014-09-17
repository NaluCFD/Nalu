/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <ExtrusionMeshDistanceBoundaryAlgorithm.h>

#include <Realm.h>
#include <FieldTypeDef.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_io
#include <stk_io/StkMeshIoBroker.hpp>

// basic c++
#include <iostream>
#include <math.h>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// ExtrusionMeshDistanceBoundaryAlgorithm - compute normal direction and
// distance for extrusion alg; uses equal area weight for extrusion direction
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ExtrusionMeshDistanceBoundaryAlgorithm::ExtrusionMeshDistanceBoundaryAlgorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : Algorithm(realm, part)
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
ExtrusionMeshDistanceBoundaryAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.fixture_->meta_data();
  stk::mesh::BulkData & bulk_data = realm_.fixture_->bulk_data();

  const int nDim = meta_data.spatial_dimension();

  //============================
  // zero out haloDxj
  //============================
  VectorFieldType *haloDxj = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "halo_dxj");

  stk::mesh::Selector s_all = stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& all_node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all );
  for ( stk::mesh::BucketVector::const_iterator ib = all_node_buckets.begin();
        ib != all_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // pointer to data
    double * hDxj = stk::mesh::field_data(*haloDxj, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      // get node
      const size_t offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j )
        hDxj[offSet+j] = 0.0;
    }
  }

  //==============================
  // compute exposed area vector
  //==============================
  GenericFieldType *exposedAreaVec = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
  VectorFieldType *coordinates = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);
  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );

  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // extract master element
    MasterElement *meFC = realm_.get_surface_master_element(b.topology());

    // extract master element specifics
    const int nodesPerElement = meFC->nodesPerElement_;
    const int numScsIp = meFC->numIntPoints_;

    // define scratch field
    std::vector<double > ws_coordinates(nodesPerElement*nDim);
    std::vector<double > ws_scs_areav(numScsIp*nDim);

    const stk::mesh::Bucket::size_type length   = b.size();
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // face data
      double * areaVec = stk::mesh::field_data(*exposedAreaVec, b, k);

      // face node relations for nodal gather
      stk::mesh::Entity const * face_node_rels = b.begin_nodes(k);

      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      int num_nodes = b.num_nodes(k);
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        double * coords = stk::mesh::field_data(*coordinates, node);
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          ws_coordinates[offSet+j] = coords[j];
        }
      }

      // compute scs integration point areavec
      double scs_error = 0.0;
      meFC->determinant(1, &ws_coordinates[0], &ws_scs_areav[0], &scs_error);

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          areaVec[offSet+j] = ws_scs_areav[offSet+j];
        }
      }
    }
  }

  //===========================
  // assemble normal to nodes
  //===========================
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // size some things that are useful
    const int num_face_nodes = b.topology().num_nodes();

    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get face
      stk::mesh::Entity face = b[k];

      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec, face);

      // face node relations for nodal gather
      stk::mesh::Entity const* face_node_rels = bulk_data.begin_nodes(face);

      // one to one mapping between ips and nodes
      for ( int ip = 0; ip < num_face_nodes; ++ip ) {

        stk::mesh::Entity node = face_node_rels[ip];
        double * hDxj = stk::mesh::field_data(*haloDxj, node);

        // offset for bip area vector
        const int faceOffSet = ip*nDim;

        double aMag = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double axj = areaVec[faceOffSet+j];
          aMag += axj*axj;
        }
        aMag = std::sqrt(aMag);

        for ( int j = 0; j < nDim; ++j ) {
          hDxj[j] += areaVec[faceOffSet+j]/aMag;
        }
      }
    }
  }

  // parallel assemble
  std::vector<stk::mesh::FieldBase*> sum_fields(1, haloDxj);
  stk::mesh::parallel_sum(bulk_data, sum_fields);

  //==============================================
  // now normalize and scale by extrusion distance
  //==============================================
  ScalarFieldType *extrusionDistance = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "extrusion_distance");
  for ( stk::mesh::BucketVector::const_iterator ib = all_node_buckets.begin();
        ib != all_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // pointer to data
    double * hDxj = stk::mesh::field_data(*haloDxj, b);
    const double * extDist = stk::mesh::field_data(*extrusionDistance, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      const size_t offSet = k*nDim;

      double dxMag = 0.0;
      for ( int j = 0; j < nDim; ++j ) {
        const double hdxj = hDxj[offSet+j];
        dxMag += hdxj*hdxj;
      }
      dxMag = std::sqrt(dxMag);
      const double fac = extDist[k]/dxMag;
      for ( int j = 0; j < nDim; ++j ) {
        hDxj[offSet+j] *= fac;
      }
    }
  }
}


} // namespace nalu
} // namespace Sierra
