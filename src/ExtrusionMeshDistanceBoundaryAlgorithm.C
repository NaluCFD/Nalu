/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <ExtrusionMeshDistanceBoundaryAlgorithm.h>

#include <Realm.h>
#include <FieldTypeDef.h>
#include <SolutionOptions.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// basic c++
#include <cmath>
#include <vector>

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

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();

  // extract extrusion factor
  const double extrusionCorrectionFac = realm_.solutionOptions_->extrusionCorrectionFac_;
  const double om_extrusioneCorrectionFac = 1.0 - extrusionCorrectionFac;

  //============================================
  // zero out haloNormal and correction factors
  //=============================================
  VectorFieldType *haloNormal 
    = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "halo_normal");
  ScalarFieldType *extDistCorrFac 
    = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "extrusion_distance_correct_fac");
  ScalarFieldType *extDistCorrCount 
    = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "extrusion_distance_correct_count");

  stk::mesh::Selector s_all = stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& all_node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all );
  for ( stk::mesh::BucketVector::const_iterator ib = all_node_buckets.begin();
        ib != all_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // pointer to data
    double * hNormal = stk::mesh::field_data(*haloNormal, b);
    double * extDCF = stk::mesh::field_data(*extDistCorrFac,b);
    double * extDCC = stk::mesh::field_data(*extDistCorrCount,b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      extDCF[k] = 0.0;
      extDCC[k] = 0.0;
      const size_t offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j )
        hNormal[offSet+j] = 0.0;
    }
  }

  //=====================================================================
  // compute exposed area vector; increment faces touching contact nodes
  //=====================================================================
  GenericFieldType *exposedAreaVec 
    = meta_data.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
  VectorFieldType *coordinates 
    = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

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

      for ( int ip = 0; ip < num_nodes; ++ip ) {
   
        // offset for bip area vector
        const int faceOffSet = ip*nDim;

        double aMag = 0.0;
        for ( int j=0; j < nDim; ++j ) {
          const double axj = ws_scs_areav[faceOffSet+j];
          areaVec[faceOffSet+j] = axj;
          aMag += axj*axj;
        }
        aMag = std::sqrt(aMag);

        // extract node; increment count; assemble normal
        stk::mesh::Entity node = face_node_rels[ip];

        // get nodal fields tocuhing this ip (nearest node); zeroed above
        double * extDCC = stk::mesh::field_data(*extDistCorrCount,node);
        double * hNormal = stk::mesh::field_data(*haloNormal, node);

        // accumulate number of faces that touch this node
        *extDCC += 1.0;

        // accumulate normal
        for ( int j = 0; j < nDim; ++j ) {
          hNormal[j] += areaVec[faceOffSet+j]/aMag;
        }
      }
    }
  }

  // parallel assemble correction count and nodal normals
  std::vector<stk::mesh::FieldBase*> sum_halo_vec;
  sum_halo_vec.push_back(extDistCorrCount);
  sum_halo_vec.push_back(haloNormal);
  stk::mesh::parallel_sum(bulk_data, sum_halo_vec);

  // normalize hNormal
  for ( stk::mesh::BucketVector::const_iterator ib = all_node_buckets.begin();
        ib != all_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // pointer to data; to be modified and constant
    double * hNormal  = stk::mesh::field_data(*haloNormal, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      const size_t offSet = k*nDim;
      double nMag = 0.0;
      for ( int j = 0; j < nDim; ++j ) {
        nMag += hNormal[offSet+j]*hNormal[offSet+j];
      }
      nMag = std::sqrt(nMag);
      for ( int j = 0; j < nDim; ++j ) {
        hNormal[offSet+j] /= nMag;
      }
    }
  }

  //==============================================
  // assemble n^nodal_k (dot) n^bip_k
  //==============================================
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

        // offset for bip area vector
        const int faceOffSet = ip*nDim;

        double aMag = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double axj = areaVec[faceOffSet+j];
          aMag += axj*axj;
        }
        aMag = std::sqrt(aMag);

        stk::mesh::Entity node = face_node_rels[ip];
        double *hNormal = stk::mesh::field_data(*haloNormal, node);
        double * eDCF  = stk::mesh::field_data(*extDistCorrFac, node);
     
        double sum = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double nj = areaVec[faceOffSet+j]/aMag;
          const double nodal_nj = hNormal[j];
          sum += nodal_nj*nj;
        }
        *eDCF += sum;
      }
    }
  }

  // parallel assemble correction factor; not yet normalized and in inverse state
  std::vector<stk::mesh::FieldBase*> sum_fields(1, extDistCorrFac);
  stk::mesh::parallel_sum(bulk_data, sum_fields);

  //==============================================
  // compute final nodal extrusion distance
  //==============================================
  VectorFieldType *haloDxj = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "halo_dxj");
  ScalarFieldType *extrusionDistance = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "extrusion_distance");
  for ( stk::mesh::BucketVector::const_iterator ib = all_node_buckets.begin();
        ib != all_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // pointer to data; to be modified and constant
    double * hDxj = stk::mesh::field_data(*haloDxj, b);
    double * eDCF  = stk::mesh::field_data(*extDistCorrFac, b);
    const double * extDist = stk::mesh::field_data(*extrusionDistance, b);
    const double * eDCC  = stk::mesh::field_data(*extDistCorrCount, b);
    const double * hNormal  = stk::mesh::field_data(*haloNormal, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // choose non-planar correction; user controled to shut it off
      const double nonPlanarCorrection 
        = eDCC[k]/eDCF[k]*extrusionCorrectionFac + om_extrusioneCorrectionFac;

      eDCF[k] = nonPlanarCorrection;
      const double fac = extDist[k]*nonPlanarCorrection;
      const size_t offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j ) {
        hDxj[offSet+j] = fac*hNormal[offSet+j];
      }
    }
  }
}


} // namespace nalu
} // namespace Sierra
