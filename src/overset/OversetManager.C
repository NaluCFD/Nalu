/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "overset/OversetManager.h"

#include "Realm.h"
#include "master_element/MasterElement.h"
#include "overset/OversetInfo.h"

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>

namespace sierra {
namespace nalu {

OversetManager::OversetManager(Realm& realm)
  : realm_(realm),
    metaData_(&realm.meta_data()),
    bulkData_(&realm_.bulk_data())
{}

OversetManager::~OversetManager()
{
  delete_info_vec();
}

void
OversetManager::delete_info_vec()
{
  for (auto ii=oversetInfoVec_.begin();
       ii != oversetInfoVec_.end(); ++ii)
    delete (*ii);

  for (auto ii=fringeInfoVec_.begin();
       ii != fringeInfoVec_.end(); ++ii)
    delete (*ii);
}

void
OversetManager::setup()
{}

stk::mesh::Selector
OversetManager::get_inactive_selector()
{
  return (stk::mesh::Selector(*inActivePart_) &
          !(stk::mesh::selectUnion(constraintPointSurfaceVecBackground_)));
}

void
OversetManager::overset_constraint_node_field_update(
  stk::mesh::FieldBase *theField,
  const int sizeRow,
  const int sizeCol)
{
  // operates on known oversetInfoVec_
  overset_constraint_node_field_update_gen(
    theField,
    sizeRow,
    sizeCol,
    oversetInfoVec_);
}

void
OversetManager::overset_constraint_node_field_update_post(
  stk::mesh::FieldBase *theField,
  const int sizeRow,
  const int sizeCol)
{
  // operates on known fringeInfoVec_
  overset_constraint_node_field_update_gen(
    theField,
    sizeRow,
    sizeCol,
    fringeInfoVec_);
}

void
OversetManager::overset_constraint_node_field_update_gen(
  stk::mesh::FieldBase *theField,
  const int sizeRow,
  const int sizeCol,
  std::vector<OversetInfo *> &infoVec)
{
  const unsigned sizeOfField = sizeRow*sizeCol;
  std::vector <double > constraintNodalQ(sizeOfField);

  // parallel communicate ghosted entities
  if ( NULL != oversetGhosting_ ) {
    std::vector< const stk::mesh::FieldBase *> fieldVec(1, theField);
    stk::mesh::communicate_field_data(*oversetGhosting_, fieldVec);
  }

  // iterate oversetInfoVec_
  std::vector<OversetInfo *>::iterator ii;
  for( ii=infoVec.begin();
       ii!=infoVec.end(); ++ii ) {

    // overset info object of interest
    OversetInfo * infoObject = (*ii);

    // extract element and node mesh object
    stk::mesh::Entity owningElement = infoObject->owningElement_;
    stk::mesh::Entity constraintNode = infoObject->constraintNode_;

    // get master element type for this contactInfo
    MasterElement *meSCS  = infoObject->meSCS_;
    const int nodesPerElement = meSCS->nodesPerElement_;
    std::vector <double > elemNodalQ(nodesPerElement*sizeRow*sizeCol);
    std::vector <double > shpfc(nodesPerElement);

    // start on gather
    stk::mesh::Entity const* elem_node_rels = bulkData_->begin_nodes(owningElement);
    const int num_nodes = bulkData_->num_nodes(owningElement);

    // sanity check on num nodes
    STK_ThrowAssert( num_nodes == nodesPerElement );

    for ( int ni = 0; ni < num_nodes; ++ni ) {
      stk::mesh::Entity node = elem_node_rels[ni];

      // load up scalar/vectors/tensor
      const double *fieldQ = (double *) stk::mesh::field_data(*theField, node );

      for ( int i = 0; i < sizeRow; ++i ) {
        const int rowI = i*sizeCol;
        const int offSetT = i*nodesPerElement*sizeRow;
        for ( int j = 0; j < sizeCol; ++j ) {
          elemNodalQ[offSetT+j*nodesPerElement+ni] = fieldQ[rowI+j];
        }
      }
    }

    // interpolate to node
    meSCS->interpolatePoint(
        sizeOfField,
        &(infoObject->isoParCoords_[0]),
        &elemNodalQ[0],
        &(constraintNodalQ[0]));

    // populate constraint node
    double *constraintQ = (double *) stk::mesh::field_data(*theField, constraintNode);
    for ( int i = 0; i < sizeRow; ++i ) {
      const int rowI = i*sizeCol;
      for ( int j = 0; j < sizeCol; ++j ) {
        constraintQ[rowI+j] = constraintNodalQ[rowI+j];
      }
    }
  }
}


}  // nalu
}  // sierra
