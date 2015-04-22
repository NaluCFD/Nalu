/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleNodalGradUEdgeContactAlgorithm.h>
#include <SolverAlgorithm.h>

#include <ContactInfo.h>
#include <ContactManager.h>
#include <FieldTypeDef.h>
#include <HaloInfo.h>
#include <LinearSystem.h>
#include <Realm.h>
#include <TimeIntegrator.h>

#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleNodalGradUEdgeContactAlgorithm - contact duidxj
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleNodalGradUEdgeContactAlgorithm::AssembleNodalGradUEdgeContactAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  VectorFieldType *velocity,
  GenericFieldType *dudx)
  : Algorithm(realm, part),
    velocity_(velocity),
    dudx_(dudx)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  // populate fieldVec; no state
  ghostFieldVec_.push_back(velocity_);
  ghostFieldVec_.push_back(dudx_);

}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNodalGradUEdgeContactAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();

  // space for interpolated right state (halo)
  std::vector<double> uNp1R(nDim);
  
  // parallel communicate ghosted entities
  if ( NULL != realm_.contactManager_->contactGhosting_ )
    stk::mesh::communicate_field_data(*(realm_.contactManager_->contactGhosting_), ghostFieldVec_);
  
  // iterate contactInfoVec_
  std::vector<ContactInfo *>::iterator ii;
  for( ii=realm_.contactManager_->contactInfoVec_.begin();
       ii!=realm_.contactManager_->contactInfoVec_.end(); ++ii ) {
    
    // get master element type for this contactInfo
    MasterElement *meSCS  = (*ii)->meSCS_;
    const int nodesPerElement = meSCS->nodesPerElement_;
    std::vector <double > elemNodalUnp1(nDim*nodesPerElement);
    std::vector <double > shpfc(nodesPerElement);

    // iterate halo face nodes
    std::map<uint64_t, HaloInfo *>::iterator iterHalo;
    for (iterHalo  = (*ii)->haloInfoMap_.begin();
         iterHalo != (*ii)->haloInfoMap_.end();
         ++iterHalo) {

      // halo info object of interest
      HaloInfo * infoObject = (*iterHalo).second;

      // extract element mesh object and global id for face node
      stk::mesh::Entity elem  = infoObject->owningElement_;

      // extract halo info stuff; can remove edge area vec...
      const double *p_areaVec = &infoObject->haloEdgeAreaVec_[0];

      stk::mesh::Entity const* elem_node_rels = bulk_data.begin_nodes(elem);
      const int num_nodes = bulk_data.num_nodes(elem);

      // now load the elemental values for future interpolation
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];
        const double *uNp1 = stk::mesh::field_data(*velocity_, node );
        for ( int j = 0; j < nDim; ++j ) {
          const int offSet = j*nodesPerElement +ni;
          elemNodalUnp1[offSet] = uNp1[j];
        }
      }

      // extract nodal fields; right state is Halo and requires interpolation
      const double volL = *stk::mesh::field_data(*dualNodalVolume_, infoObject->faceNode_);
      double *dudxL = stk::mesh::field_data(*dudx_, infoObject->faceNode_);

      const double *uNp1L = stk::mesh::field_data(*velocity_, infoObject->faceNode_);
      const int sizeOfVectorField = nDim;
      meSCS->interpolatePoint(
        sizeOfVectorField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalUnp1[0],
        &(uNp1R[0]));

      // recall that halo node is "right" state; face node is "left"
      for ( int i = 0; i < nDim; ++ i ) {
        const double uIp = 0.5*(uNp1L[i] + uNp1R[i]);
        const int row_gradQ = i*nDim;
        for (int j = 0; j < nDim; ++j ) {
          const double ajUip = uIp*p_areaVec[j];
          dudxL[row_gradQ+j] += ajUip/volL;
        }
      }
    }
  }
}


} // namespace nalu
} // namespace Sierra
