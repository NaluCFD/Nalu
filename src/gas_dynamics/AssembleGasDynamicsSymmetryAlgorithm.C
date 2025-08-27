/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <gas_dynamics/AssembleGasDynamicsSymmetryAlgorithm.h>
#include <Realm.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleGasDynamicsSymmetryAlgorithm - assembles RHS for gas dynamics; sym
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleGasDynamicsSymmetryAlgorithm::AssembleGasDynamicsSymmetryAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *pressure,
  GenericFieldType *rhsGasDyn)
  : Algorithm(realm, part),
    pressure_(pressure),
    rhsGasDyn_(rhsGasDyn),
    exposedAreaVec_(NULL)
{
  // save off mising fields
  stk::mesh::MetaData & metaData = realm_.meta_data();
  exposedAreaVec_ = metaData.get_field<double>(metaData.side_rank(), "exposed_area_vector");
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleGasDynamicsSymmetryAlgorithm::execute()
{
  stk::mesh::MetaData & metaData = realm_.meta_data();

  // sizes
  const int nDim = metaData.spatial_dimension();

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = metaData.locally_owned_part()
    & stk::mesh::selectUnion(partVec_) 
    & !(realm_.get_inactive_selector());

  //===========================================================
  // assemble edge-based flux operator to the node
  //===========================================================

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( metaData.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());

    // extract master element specifics
    const int numScsIp = meFC->numIntPoints_;
    const int *ipNodeMap = meFC->ipNodeMap();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, b, k);

      stk::mesh::Entity const * face_node_rels = b.begin_nodes(k);

      // start assembly
      for ( int ip = 0; ip < numScsIp; ++ip ) {
        
        // nearest node
        const int nn = ipNodeMap[ip];
        
        stk::mesh::Entity nodeNN = face_node_rels[nn];

        // pointer to fields to assemble
        double *rhsGasDynNN = stk::mesh::field_data(*rhsGasDyn_, nodeNN);

        // suplemental
        double pressureNN = *stk::mesh::field_data(*pressure_, nodeNN);
        
        for ( int i = 0; i < nDim; ++i ) {
          rhsGasDynNN[i] -= pressureNN*areaVec[ip*nDim+i];
        }
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
