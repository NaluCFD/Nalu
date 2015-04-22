/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleNodalGradEdgeContactAlgorithm.h>
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
// AssembleNodalGradEdgeContactAlgorithm - contact scalar dqdx
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleNodalGradEdgeContactAlgorithm::AssembleNodalGradEdgeContactAlgorithm(
    Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *scalarQ,
  VectorFieldType *dqdx)
  : Algorithm(realm, part),
    scalarQ_(scalarQ),
    dqdx_(dqdx)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  // populate fieldVec
  ghostFieldVec_.push_back(scalarQ_);
  ghostFieldVec_.push_back(dqdx_);

}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleNodalGradEdgeContactAlgorithm::execute()
{

  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  const int nDim = meta_data.spatial_dimension();

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
    std::vector <double > elemNodalScalarQ(nodesPerElement);
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

      stk::mesh::Entity const *elem_node_rels = bulk_data.begin_nodes(elem);
      const unsigned num_nodes = bulk_data.num_nodes(elem);

      // now load the elemental values for future interpolation
      for ( unsigned ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = elem_node_rels[ni];
        elemNodalScalarQ[ni] = *stk::mesh::field_data(*scalarQ_, node);
      }

      // extract nodal fields; right state is Halo and requires inperpolation
      const double scalarQL  = *stk::mesh::field_data(*scalarQ_, infoObject->faceNode_);
      const int sizeOfScalarField = 1;
      double scalarQR = 0.0;
      meSCS->interpolatePoint(
        sizeOfScalarField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalScalarQ[0],
        &scalarQR);

      // extract face values
      const double volL = *stk::mesh::field_data(*dualNodalVolume_, infoObject->faceNode_);
      double *dqdxL = stk::mesh::field_data(*dqdx_, infoObject->faceNode_);

      // recall that halo node is "right" state; face node is "left"
      const double scalarQIp = 0.5*(scalarQL + scalarQR);
      for (int j = 0; j < nDim; ++j ) {
        const double scalarQ_ax = scalarQIp*p_areaVec[j];
        dqdxL[j] += scalarQ_ax/volL;
      }
    }
  }
}


} // namespace nalu
} // namespace Sierra
