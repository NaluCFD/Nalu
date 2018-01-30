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

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleElemSolverAlgorithm::execute()
{
  stk::mesh::BulkData & bulk_data = realm_.bulk_data();

  // set any data
  const size_t activeKernelsSize = activeKernels_.size();
  for ( size_t i = 0; i < activeKernelsSize; ++i )
    activeKernels_[i]->setup(*realm_.timeIntegrator_);

  run_algorithm(bulk_data, [&](SharedMemData& smdata, unsigned nodesPerEntity, const stk::mesh::Entity* entityNodes[])
  {
      set_zero(smdata.simdrhs.data(), smdata.simdrhs.size());
      set_zero(smdata.simdlhs.data(), smdata.simdlhs.size());

      // call supplemental; gathers happen inside the elem_execute method
      for ( size_t i = 0; i < activeKernelsSize; ++i )
        activeKernels_[i]->execute( smdata.simdlhs, smdata.simdrhs, smdata.simdPrereqData );

      for(int simdElemIndex=0; simdElemIndex<smdata.numSimdElems; ++simdElemIndex) {
        extract_vector_lane(smdata.simdrhs, simdElemIndex, smdata.rhs);
        extract_vector_lane(smdata.simdlhs, simdElemIndex, smdata.lhs);
        apply_coeff(nodesPerEntity, entityNodes[simdElemIndex],
                    smdata.scratchIds, smdata.sortPermutation, smdata.rhs, smdata.lhs, __FILE__);
      }
  });
}

} // namespace nalu
} // namespace Sierra
