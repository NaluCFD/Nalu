/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleNodalGradUEdgeAlgorithm.h>
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

//==========================================================================
// Class Definition
//==========================================================================
// AssembleNodalGradUEdgeAlgorithm - base class equation system
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleNodalGradUEdgeAlgorithm::AssembleNodalGradUEdgeAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  VectorFieldType *velocity,
  GenericFieldType *dudx)
  : Algorithm(realm, part),
    velocity_(velocity),
    dudx_(dudx)
{
  // extract fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  edgeAreaVec_ = meta_data.get_field<VectorFieldType>(stk::topology::EDGE_RANK, "edge_area_vector");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNodalGradUEdgeAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    & stk::mesh::selectUnion(partVec_)  
    & !(realm_.get_inactive_selector());

  //===========================================================
  // assemble edge-based gradient operator to the node
  //===========================================================

  stk::mesh::BucketVector const& edge_buckets =
    realm_.get_buckets( stk::topology::EDGE_RANK, s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = edge_buckets.begin();
        ib != edge_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // pointer to edge area vector
    double * av = stk::mesh::field_data(*edgeAreaVec_, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      stk::mesh::Entity const * edge_node_rels = b.begin_nodes(k);

      // sanity check on number or nodes
      ThrowAssert( b.num_nodes(k) == 2 );

      // left and right nodes
      stk::mesh::Entity nodeL = edge_node_rels[0];
      stk::mesh::Entity nodeR = edge_node_rels[1];

      // dudx phi at nodes
      double * dudxL = stk::mesh::field_data( *dudx_, nodeL );
      double * dudxR = stk::mesh::field_data( *dudx_, nodeR );

      // dual volume at nodes
      const double volL = *stk::mesh::field_data( *dualNodalVolume_, nodeL );
      const double volR = *stk::mesh::field_data( *dualNodalVolume_, nodeR );

      // velocity at nodes
      double *uL = stk::mesh::field_data( *velocity_, nodeL );
      double *uR = stk::mesh::field_data( *velocity_, nodeR );

      // start the work...
      const size_t offSet = k*nDim;

      const double invVolL = 1.0/volL;
      const double invVolR = 1.0/volR;

      int counter = 0;
      for ( int i = 0; i < nDim; ++i ) {
        const double uip = 0.5*(uL[i] + uR[i]);
        for ( int j = 0; j < nDim; ++j ) {
          const double aj = av[offSet+j];
          const double ajUip = aj*uip;
          const int cj = counter++;
          dudxL[cj] += ajUip*invVolL;
          dudxR[cj] -= ajUip*invVolR;
        }
      }

    }
  }
}

} // namespace nalu
} // namespace Sierra
