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
  stk::mesh::EntityRank entityRank,
  unsigned nodesPerEntity,
  bool interleaveMEViews)
  : SolverAlgorithm(realm, part, eqSystem),
    entityRank_(entityRank),
    nodesPerEntity_(nodesPerEntity),
    rhsSize_(nodesPerEntity*eqSystem->linsys_->numDof()),
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

static constexpr int simdLen = stk::simd::ndoubles;

int
calculate_shared_mem_bytes_per_thread(int lhsSize, int rhsSize, int scratchIdsSize, int nDim,
                                      ElemDataRequests& dataNeededByKernels)
{
    int bytes_per_thread = (rhsSize + lhsSize)*sizeof(double) + (2*scratchIdsSize)*sizeof(int) +
                           get_num_bytes_pre_req_data<double>(dataNeededByKernels, nDim);
    bytes_per_thread *= 2*simdLen;
    return bytes_per_thread;
}

template<typename T>
void set_zero(T* values, unsigned length)
{
    for(unsigned i=0; i<length; ++i) {
        values[i] = 0;
    }
}

int get_next_num_elems_simd(int bktIndex, int bktLength)
{
  int numElems = simdLen;
  if (bktLength - bktIndex*simdLen < simdLen) {
    numElems = bktLength - bktIndex*simdLen;
  }
  if (numElems < 0 || numElems > simdLen) {
    std::cout<<"ERROR, simdElems="<<numElems<<" shouldn't happen!!"<<std::endl;
    numElems = 0;
  }
  return numElems;
}

size_t get_simd_bucket_length(size_t bktLength)
{
    size_t simdBucketLen = bktLength/simdLen;
    const size_t remainder = bktLength%simdLen;
    if (remainder > 0) {
      simdBucketLen += 1;
    }
    return simdBucketLen;
}

struct SharedMemData {
    SharedMemData(const sierra::nalu::TeamHandleType& team,
         const stk::mesh::BulkData& bulk,
         const ElemDataRequests& dataNeededByKernels,
         unsigned nodesPerEntity,
         unsigned rhsSize)
     : simdPrereqData(team, bulk, nodesPerEntity, dataNeededByKernels)
    {
        for(int simdIndex=0; simdIndex<simdLen; ++simdIndex) {
          prereqData[simdIndex] = std::unique_ptr<ScratchViews<double> >(new ScratchViews<double>(team, bulk, nodesPerEntity, dataNeededByKernels));
        }
        simdrhs = get_shmem_view_1D<DoubleType>(team, rhsSize);
        simdlhs = get_shmem_view_2D<DoubleType>(team, rhsSize, rhsSize);
        rhs = get_shmem_view_1D<double>(team, rhsSize);
        lhs = get_shmem_view_2D<double>(team, rhsSize, rhsSize);

        scratchIds = get_int_shmem_view_1D(team, rhsSize);
        sortPermutation = get_int_shmem_view_1D(team, rhsSize);
    }

    std::unique_ptr<ScratchViews<double>> prereqData[simdLen];
    ScratchViews<DoubleType> simdPrereqData;
    SharedMemView<DoubleType*> simdrhs;
    SharedMemView<DoubleType**> simdlhs;
    SharedMemView<double*> rhs;
    SharedMemView<double**> lhs;

    SharedMemView<int*> scratchIds;
    SharedMemView<int*> sortPermutation;
};

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
  const int bytes_per_thread = calculate_shared_mem_bytes_per_thread(lhsSize, rhsSize_, scratchIdsSize,
                                                                   meta_data.spatial_dimension(), dataNeededByKernels_);

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union =
          meta_data.locally_owned_part() & stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& elem_buckets =
          realm_.get_buckets(entityRank_, s_locally_owned_union );

  auto team_exec = sierra::nalu::get_team_policy(elem_buckets.size(), bytes_per_team, bytes_per_thread);
  Kokkos::parallel_for(team_exec, [&](const sierra::nalu::TeamHandleType& team)
  {
    stk::mesh::Bucket & b = *elem_buckets[team.league_rank()];
    
    ThrowAssertMsg(b.topology().num_nodes() == (unsigned)nodesPerEntity_,
                   "AssembleElemSolverAlgorithm expected nodesPerEntity_ = "
                   <<nodesPerEntity_<<", but b.topology().num_nodes() = "<<b.topology().num_nodes());

    SharedMemData smdata(team, bulk_data, dataNeededByKernels_, nodesPerEntity_, rhsSize_);

    const stk::mesh::Entity* entityNodes[simdLen];
    const size_t bucketLen   = b.size();
    const size_t simdBucketLen = get_simd_bucket_length(bucketLen);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, simdBucketLen), [&](const size_t& bktIndex)
    {
      int numSimdElems = get_next_num_elems_simd(bktIndex, bucketLen);

      for(int simdElemIndex=0; simdElemIndex<numSimdElems; ++simdElemIndex) {
        stk::mesh::Entity element = b[bktIndex*simdLen + simdElemIndex];
        entityNodes[simdElemIndex] = bulk_data.begin_nodes(element);
        fill_pre_req_data(dataNeededByKernels_, bulk_data, element,
                          *smdata.prereqData[simdElemIndex], interleaveMEViews_);
      }

      copy_and_interleave(smdata.prereqData, numSimdElems, smdata.simdPrereqData, interleaveMEViews_);

      if (!interleaveMEViews_) {
        fill_master_element_views(dataNeededByKernels_, bulk_data, smdata.simdPrereqData);
      }

      set_zero(smdata.simdrhs.data(), smdata.simdrhs.size());
      set_zero(smdata.simdlhs.data(), smdata.simdlhs.size());

      // call supplemental; gathers happen inside the elem_execute method
      for ( size_t i = 0; i < activeKernelsSize; ++i )
        activeKernels_[i]->execute( smdata.simdlhs, smdata.simdrhs, smdata.simdPrereqData );

      for(int simdElemIndex=0; simdElemIndex<numSimdElems; ++simdElemIndex) {
        extract_vector_lane(smdata.simdrhs, simdElemIndex, smdata.rhs);
        extract_vector_lane(smdata.simdlhs, simdElemIndex, smdata.lhs);
        apply_coeff(nodesPerEntity_, entityNodes[simdElemIndex],
                    smdata.scratchIds, smdata.sortPermutation, smdata.rhs, smdata.lhs, __FILE__);
      }
    });
  });
}

} // namespace nalu
} // namespace Sierra
