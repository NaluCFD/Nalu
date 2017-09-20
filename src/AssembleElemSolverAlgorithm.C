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
#include <Kernel.h>
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
#include <CopyAndInterleave.h>

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
  const stk::topology &theTopo,
  bool interleaveMEViews)
  : SolverAlgorithm(realm, part, eqSystem),
    topo_(theTopo),
    rhsSize_(theTopo.num_nodes()*eqSystem->linsys_->numDof()),
    interleaveMEViews_(interleaveMEViews)
{
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

  // set any data
  const size_t activeKernelsSize = activeKernels_.size();
  for ( size_t i = 0; i < activeKernelsSize; ++i )
    activeKernels_[i]->setup(*realm_.timeIntegrator_);

  const int lhsSize = rhsSize_*rhsSize_;
  const int scratchIdsSize = rhsSize_;

  // fixed size for this homogeneous algorithm
  const int bytes_per_team = 0;
  // Bytes per thread =
  //    scratch views size + LHS + RHS + node IDs + padding for alignment
  int bytes_per_thread =
    (rhsSize_ + lhsSize)*sizeof(double) + (2*scratchIdsSize)*sizeof(int) +
    get_num_bytes_pre_req_data<double>(dataNeededByKernels_, meta_data.spatial_dimension());
  constexpr int simdLen = stk::simd::ndoubles;
  bytes_per_thread *= 3*simdLen;//3 is wrong, still need to debug this!!! should be 2 but seems to have memory problems.

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );

  auto team_exec = sierra::nalu::get_team_policy(elem_buckets.size(), bytes_per_team, bytes_per_thread);
  Kokkos::parallel_for(team_exec, [&](const sierra::nalu::TeamHandleType& team)
  {
    stk::mesh::Bucket & b = *elem_buckets[team.league_rank()];
    
    ThrowAssertMsg(b.topology() == topo_,"topo_ = "<<topo_<<", b.topology() = "<<b.topology());

    std::vector<sierra::nalu::ScratchViews<double>*> prereqData(simdLen, nullptr);

    for(int simdIndex=0; simdIndex<simdLen; ++simdIndex) {
      prereqData[simdIndex] = new sierra::nalu::ScratchViews<double>(team, bulk_data, topo_, dataNeededByKernels_);
    }

    sierra::nalu::ScratchViews<DoubleType> simdPrereqData(team, bulk_data, topo_, dataNeededByKernels_);

    SharedMemView<DoubleType*> simdrhs = get_shmem_view_1D<DoubleType>(team, rhsSize_);
    SharedMemView<DoubleType**> simdlhs = get_shmem_view_2D<DoubleType>(team, rhsSize_, rhsSize_);
    SharedMemView<double*> rhs = get_shmem_view_1D<double>(team, rhsSize_);
    SharedMemView<double**> lhs = get_shmem_view_2D<double>(team, rhsSize_, rhsSize_);

    SharedMemView<int*> scratchIds = get_int_shmem_view_1D(team, scratchIdsSize);
    SharedMemView<int*> sortPermutation = get_int_shmem_view_1D(team, scratchIdsSize);

    const stk::mesh::Entity* elemNodes[simdLen];
    unsigned num_nodes = b.topology().num_nodes();
    const stk::mesh::Bucket::size_type length   = b.size();
    stk::mesh::Bucket::size_type simdBucketLen = length/simdLen;
    const stk::mesh::Bucket::size_type remainder = length%simdLen;
    if (remainder > 0) {
      simdBucketLen += 1;
    }


    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, simdBucketLen), [&](const size_t& bktIndex)
    {
      int simdElems = simdLen;
      if (length - bktIndex*simdLen < simdLen) {
          simdElems = length - bktIndex*simdLen;
      }
      if (simdElems < 0 || simdElems > simdLen) {
         std::cout<<"ERROR, simdElems="<<simdElems<<" shouldn't happen!!"<<std::endl;
         simdElems = 0;
      }

      stk::mesh::Entity element;
      for(int simdElemIndex=0; simdElemIndex<simdElems; ++simdElemIndex) {
        // get element
        element = b[bktIndex*simdLen + simdElemIndex];
        elemNodes[simdElemIndex] = bulk_data.begin_nodes(element);
        fill_pre_req_data(dataNeededByKernels_, bulk_data, topo_, element,
                          *prereqData[simdElemIndex], interleaveMEViews_);
      }

      copy_and_interleave(prereqData, simdElems, simdPrereqData, interleaveMEViews_);

      if (!interleaveMEViews_) {
        fill_master_element_views(dataNeededByKernels_, bulk_data, topo_, element, simdPrereqData);
      }

      for ( int i = 0; i < rhsSize_; ++i ) {
        simdrhs(i) = 0.0;
      }

      for ( int i = 0; i < rhsSize_; ++i ) {
        for ( int j = 0; j < rhsSize_; ++j ) {
          simdlhs(i,j) = 0.0;
        }
      }

      // call supplemental; gathers happen inside the elem_execute method
      for ( size_t i = 0; i < activeKernelsSize; ++i )
        activeKernels_[i]->execute( simdlhs, simdrhs, simdPrereqData );

      for(int simdElemIndex=0; simdElemIndex<simdElems; ++simdElemIndex) {
        extract_vector_lane(simdrhs, simdElemIndex, rhs);
        extract_vector_lane(simdlhs, simdElemIndex, lhs);
        apply_coeff(num_nodes, elemNodes[simdElemIndex], scratchIds, sortPermutation, rhs, lhs, __FILE__);
      }
    });

    for(int simdIndex=0; simdIndex<simdLen; ++simdIndex) {
       delete prereqData[simdIndex];
    }
  });
}

} // namespace nalu
} // namespace Sierra
