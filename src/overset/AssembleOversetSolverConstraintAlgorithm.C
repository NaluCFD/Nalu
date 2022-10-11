/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <overset/AssembleOversetSolverConstraintAlgorithm.h>
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
// AssembleOversetSolverConstraintAlgorithm - add LHS/RHS for scalar
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleOversetSolverConstraintAlgorithm::AssembleOversetSolverConstraintAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  stk::mesh::FieldBase *fieldQ,
  const bool densityScaling)
  : SolverAlgorithm(realm, part, eqSystem),
    fieldQ_(fieldQ),
    density_(nullptr),
    dualNodalVolume_(nullptr),
    scaleFac_(densityScaling ? 1.0 : 0.0 )
{
  // populate fieldVec
  ghostFieldVec_.push_back(fieldQ_);

  // extract fields (only used for density scaling)
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  density_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "density");
  dualNodalVolume_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "dual_nodal_volume");
  if ( densityScaling ) {
    ghostFieldVec_.push_back(density_);
    ghostFieldVec_.push_back(dualNodalVolume_);    
    NaluEnv::self().naluOutputP0() << "Overset density scaling active " << scaleFac_ << std::endl;
  }
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleOversetSolverConstraintAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildOversetNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleOversetSolverConstraintAlgorithm::execute()
{
  // first thing to do is to zero out the row (lhs and rhs)
  prepare_constraints();

  // parallel communicate ghosted entities
  if ( NULL != realm_.oversetManager_->oversetGhosting_ )
    stk::mesh::communicate_field_data(*(realm_.oversetManager_->oversetGhosting_), ghostFieldVec_);

  // extract the rank
  const int theRank = NaluEnv::self().parallel_rank();

  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  const int nDim = realm_.meta_data().spatial_dimension();

  // space for LHS/RHS (nodesPerElem+1)*numDof*(nodesPerElem+1)*numDof; (nodesPerElem+1)*numDof
  std::vector<double> lhs;
  std::vector<double> rhs;
  std::vector<int> scratchIds;
  std::vector<double> scratchVals;
  std::vector<stk::mesh::Entity> connected_nodes;

  // master element data
  std::vector <double > ws_general_shape_function;
 
  // interpolate nodal values to point-in-elem
  const int sizeOfDof = eqSystem_->linsys_->numDof();
 
  // size interpolated value
  std::vector<double> qNp1Constraint(sizeOfDof, 0.0);

  // Parallel communication of ghosted entities has been already handled in
  // EquationSystems::pre_iter_work

  // iterate oversetInfoVec_
  std::vector<OversetInfo *>::iterator ii;
  for( ii=realm_.oversetManager_->oversetInfoVec_.begin();
       ii!=realm_.oversetManager_->oversetInfoVec_.end(); ++ii ) {

    // overset info object of interest
    OversetInfo * infoObject = (*ii);

    // extract element and node mesh object
    stk::mesh::Entity owningElement = infoObject->owningElement_;
    stk::mesh::Entity constraintNode = infoObject->constraintNode_;

    // extract the owning rank for this node
    const int nodeRank = bulkData.parallel_owner_rank(constraintNode);
    
    // check to see if this node is locally owned by this rank; we only want to process locally owned nodes, not shared
    if ( theRank != nodeRank )
      continue;

    // get master element type for this contactInfo
    MasterElement *meSCS  = infoObject->meSCS_;
    const int nodesPerElement = meSCS->nodesPerElement_;
    std::vector <double > elemNodalQ(nodesPerElement*sizeOfDof);
    std::vector <double > shpfc(nodesPerElement);

    // resize some things; matrix related
    const int npePlusOne = nodesPerElement+1;
    const int lhsSize = npePlusOne*sizeOfDof*npePlusOne*sizeOfDof;
    const int rhsSize = npePlusOne*sizeOfDof;
    lhs.resize(lhsSize);
    rhs.resize(rhsSize);
    scratchIds.resize(rhsSize);
    scratchVals.resize(rhsSize);
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
    const double *qNp1Nodal = (double *)stk::mesh::field_data(*fieldQ_, constraintNode);
    
    stk::mesh::Entity const* elem_node_rels = bulkData.begin_nodes(owningElement);
    const int num_nodes = bulkData.num_nodes(owningElement);

    // assemble element average
    double elemRho = 0.0;
    double elemDnv = 0.0;

    // now load the elemental values for future interpolation; fill in connected nodes; first connected node is orhpan
    connected_nodes[0] = constraintNode;
    for ( int ni = 0; ni < num_nodes; ++ni ) {
      stk::mesh::Entity node = elem_node_rels[ni];
      connected_nodes[ni+1] = node;

      const double *qNp1 = (double *)stk::mesh::field_data(*fieldQ_, node );
      for ( int i = 0; i < sizeOfDof; ++i ) {
        elemNodalQ[i*nodesPerElement + ni] = qNp1[i];
      }
      
      // density scaling fields
      const double rho = *stk::mesh::field_data(*density_, node );
      const double dnv = *stk::mesh::field_data(*dualNodalVolume_, node );
      elemRho += rho;
      elemDnv += dnv;
    }
    elemRho /= num_nodes;
    elemDnv /= num_nodes;

    // continuity equation may benifit from scaling when using VoF BF
    const double systemScaling = scaleFac_*std::pow(elemDnv, 1.0/nDim)/elemRho + (1.0 - scaleFac_);

    // interpolate dof to elemental ips (assigns qNp1)
    meSCS->interpolatePoint(
      sizeOfDof,
      &(infoObject->isoParCoords_[0]),
      &elemNodalQ[0],
      &qNp1Constraint[0]);
    
    // lhs; extract general shape function
    meSCS->general_shape_fcn(1, &(infoObject->isoParCoords_[0]), &ws_general_shape_function[0]);

    // rhs; constraint node is defined to be the zeroth connected node
    for ( int i = 0; i < sizeOfDof; ++i) {
      const int rowOi = i * npePlusOne * sizeOfDof;
      const double residual = qNp1Nodal[i] - qNp1Constraint[i];
      p_rhs[i] = -residual*systemScaling;

      // row is zero by design (first connected node is the constraint node); assign it fully
      p_lhs[rowOi+i] += 1.0*systemScaling;

      for ( int ic = 0; ic < nodesPerElement; ++ic ) {
        const int indexR = i + sizeOfDof*(ic+1);
        const int rOiR = rowOi+indexR;
        p_lhs[rOiR] -= ws_general_shape_function[ic]*systemScaling;
      }
    }

    // apply to linear system
    apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);
  }
}

//--------------------------------------------------------------------------
//-------- prepare_constraints ---------------------------------------------
//--------------------------------------------------------------------------
void
AssembleOversetSolverConstraintAlgorithm::prepare_constraints()
{
  const int sysStart = 0;
  const int sysEnd = eqSystem_->linsys_->numDof();
  eqSystem_->linsys_->prepareConstraints(sysStart,sysEnd);
}

} // namespace nalu
} // namespace Sierra
