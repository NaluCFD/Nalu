/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleElemSolverAlgorithm.h>
#include <EquationSystem.h>
#include <SolverAlgorithm.h>
#include <master_element/MasterElement.h>

#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <Realm.h>
#include <SupplementalAlgorithm.h>
#include <TimeIntegrator.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleElemSolverAlgorithm - add LHS/RHS for element-based contribution
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleElemSolverAlgorithm::AssembleElemSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  const stk::topology &theTopo)
  : SolverAlgorithm(realm, part, eqSystem)
{
  // size some things; matrix related
  const int nodesPerElement = realm.get_surface_master_element(theTopo)->nodesPerElement_;
  const int rhsSize = nodesPerElement*(eqSystem->linsys_->numDof());
  lhs_.resize(rhsSize*rhsSize);
  rhs_.resize(rhsSize);
  scratchIds_.resize(rhsSize);
  scratchVals_.resize(rhsSize);
  connectedNodes_.resize(nodesPerElement);
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleElemSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildElemToNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleElemSolverAlgorithm::execute()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // set any data
  const size_t supplementalAlgSize = supplementalAlg_.size();
  for ( size_t i = 0; i < supplementalAlgSize; ++i )
    supplementalAlg_[i]->setup();

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // get element
      stk::mesh::Entity element = b[k];

      // extract node relations and provide connected nodes
      stk::mesh::Entity const * node_rels = b.begin_nodes(k);
      int num_nodes = b.num_nodes(k);

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];
        // set connected nodes
        connectedNodes_[ni] = node;
      }

      for ( size_t i = 0; i < lhs_.size(); ++i )
        lhs_[i] = 0.0;
      for ( size_t i = 0; i < rhs_.size(); ++i )
        rhs_[i] = 0.0;

      // call supplemental; gathers happen inside the elem_execute method
      for ( size_t i = 0; i < supplementalAlgSize; ++i )
        supplementalAlg_[i]->element_execute( &lhs_[0], &rhs_[0], element );
      
      apply_coeff(connectedNodes_, scratchIds_, scratchVals_, rhs_, lhs_, __FILE__);
    }
  }
}

} // namespace nalu
} // namespace Sierra
