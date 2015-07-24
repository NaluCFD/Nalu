/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <overset/AssembleScalarOversetSolverAlgorithm.h>
#include <EquationSystem.h>
#include <SolverAlgorithm.h>

#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <Realm.h>
#include <TimeIntegrator.h>

// master element
#include <master_element/MasterElement.h>

// overset
#include <overset/OversetManager.h>
#include <overset/OversetInfo.h>

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
// AssembleScalarOversetSolverAlgorithm - add LHS/RHS for scalar
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleScalarOversetSolverAlgorithm::AssembleScalarOversetSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  ScalarFieldType *scalarQ)
  : SolverAlgorithm(realm, part, eqSystem),
    scalarQ_(scalarQ)
{
  // populate fieldVec
  ghostFieldVec_.push_back(scalarQ_);
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleScalarOversetSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildOversetNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleScalarOversetSolverAlgorithm::execute()
{
  // first thing to do is to zero out the row (lhs and rhs)
  prepare_constraints();
  
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  // space for LHS/RHS (nodesPerElem+1)*(nodesPerElem+1); nodesPerElem+1
  std::vector<double> lhs;
  std::vector<double> rhs;
  std::vector<stk::mesh::Entity> connected_nodes;

  // master element data
  std::vector <double > ws_general_shape_function;
 
  // deal with state
  ScalarFieldType &scalarQNp1  = scalarQ_->field_of_state(stk::mesh::StateNP1);

  // interpolate nodal values to point-in-elem
  const int sizeOfScalarField = 1;
 
  // parallel communicate ghosted entities
  if ( NULL != realm_.oversetManager_->oversetGhosting_ )
    stk::mesh::communicate_field_data(*(realm_.oversetManager_->oversetGhosting_), ghostFieldVec_);  

  // iterate oversetInfoVec_
  std::vector<OversetInfo *>::iterator ii;
  for( ii=realm_.oversetManager_->oversetInfoVec_.begin();
       ii!=realm_.oversetManager_->oversetInfoVec_.end(); ++ii ) {

    // overset info object of interest
    OversetInfo * infoObject = (*ii);

    // extract element and node mesh object
    stk::mesh::Entity owningElement = infoObject->owningElement_;
    stk::mesh::Entity orphanNode = infoObject->orphanNode_;

    // get master element type for this contactInfo
    MasterElement *meSCS  = infoObject->meSCS_;
    const int nodesPerElement = meSCS->nodesPerElement_;
    std::vector <double > elemNodalQ(nodesPerElement);
    std::vector <double > shpfc(nodesPerElement);

    // resize some things; matrix related
    const int npePlusOne = nodesPerElement+1;
    const int lhsSize = npePlusOne*npePlusOne;
    const int rhsSize = npePlusOne;
    lhs.resize(lhsSize);
    rhs.resize(rhsSize);
    connected_nodes.resize(npePlusOne);

    // algorithm related; element 
    ws_general_shape_function.resize(nodesPerElement);
   
    // pointer to lhs/rhs
    double *p_lhs = &lhs[0];
    double *p_rhs = &rhs[0];
    
    // zeroing of lhs/rhs
    for ( int k = 0; k < lhsSize; ++k ) {
      p_lhs[k] = 0.0;
    }
    for ( int k = 0; k < rhsSize; ++k ) {
      p_rhs[k] = 0.0;
    }
    
    // extract nodal value for scalarQ
    const double qNp1Nodal = *stk::mesh::field_data(scalarQNp1, orphanNode);
    
    stk::mesh::Entity const* elem_node_rels = bulk_data.begin_nodes(owningElement);
    const int num_nodes = bulk_data.num_nodes(owningElement);

    // now load the elemental values for future interpolation; fill in connected nodes; first connected node is orhpan
    connected_nodes[0] = orphanNode;
    for ( int ni = 0; ni < num_nodes; ++ni ) {
      stk::mesh::Entity node = elem_node_rels[ni];
      connected_nodes[ni+1] = node;
      elemNodalQ[ni] = *stk::mesh::field_data(scalarQNp1, node);
    }

    // interpolate dof to elemental ips
    double qNp1Orphan = 0.0;
    meSCS->interpolatePoint(
      sizeOfScalarField,
      &(infoObject->isoParCoords_[0]),
      &elemNodalQ[0],
      &qNp1Orphan);
    
    // rhs...
    const double residual = qNp1Nodal - qNp1Orphan;
    p_rhs[0] = -residual;
    
    // lhs; extract general shape function
    meSCS->general_shape_fcn(1, &(infoObject->isoParCoords_[0]), &ws_general_shape_function[0]);
    
    // row is zero by design (first connected node is the orphan node)
    int rowR = 0;
    p_lhs[0] += 1.0;
    for ( int ic = 0; ic < nodesPerElement; ++ic ) {
      p_lhs[ic+1] -= ws_general_shape_function[ic];
    }
    
    // apply to linear system
    apply_coeff(connected_nodes, rhs, lhs, __FILE__);
  }
}

//--------------------------------------------------------------------------
//-------- prepare_constraints ---------------------------------------------
//--------------------------------------------------------------------------
void
AssembleScalarOversetSolverAlgorithm::prepare_constraints()
{
  eqSystem_->linsys_->prepareConstraints(0,1);
}

} // namespace nalu
} // namespace Sierra
