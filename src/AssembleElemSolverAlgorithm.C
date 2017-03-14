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

// stk topo
#include <stk_topology/topology.hpp>

#include <KokkosInterface.h>
#include <ScratchViews.h>

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
  : SolverAlgorithm(realm, part, eqSystem),
    topo_(theTopo)
{
  // size some things; matrix related; LHS/RHS not yet a View..
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
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  stk::mesh::FieldBase* coordField = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  // set any data
  const size_t supplementalAlgSize = supplementalAlg_.size();
  for ( size_t i = 0; i < supplementalAlgSize; ++i )
    supplementalAlg_[i]->setup();

  // fixed size for this homogeneous algorithm
  const int bytes_per_team = 0;
  const int bytes_per_thread = get_num_bytes_pre_req_data(dataNeededBySuppAlgs_, meta_data.spatial_dimension());

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );

  auto team_exec = sierra::nalu::get_team_policy(elem_buckets.size(), bytes_per_team, bytes_per_thread);
  Kokkos::parallel_for(team_exec, [&](const sierra::nalu::TeamHandleType& team)
  {
    stk::mesh::Bucket & b = *elem_buckets[team.league_rank()];
    
    ThrowAssert(b.topology() == topo_);

    sierra::nalu::ScratchViews prereqData(team, bulk_data, topo_, dataNeededBySuppAlgs_);

    const stk::mesh::Bucket::size_type length   = b.size();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, length), [&](const size_t& k)
    {
      // get element
      stk::mesh::Entity element = b[k];
      fill_pre_req_data(dataNeededBySuppAlgs_, bulk_data, topo_, element,
                        coordField, prereqData);

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
        supplementalAlg_[i]->element_execute( &lhs_[0], &rhs_[0], element, prereqData );
      
      apply_coeff(connectedNodes_, scratchIds_, scratchVals_, rhs_, lhs_, __FILE__);
    });
  });
}

} // namespace nalu
} // namespace Sierra
