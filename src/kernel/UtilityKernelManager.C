/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <kernel/UtilityKernelManager.h>

#include <kernel/NodalUtilityKernelManager.h>

#include <EquationSystem.h>
#include <SolverAlgorithm.h>
#include <master_element/MasterElement.h>

#include <NodalScratchData.h>

#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <Realm.h>
#include <kernel/Kernel.h>
#include <TimeIntegrator.h>

#include <master_element/TensorOps.h>

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

UtilityKernelManager::UtilityKernelManager(
  const stk::mesh::BulkData& bulk,
  const SolutionOptions& solnOpts,
  const TimeIntegrator& timeIntegrator,
  ElemDataRequests& in,
  ElemDataRequests& out)
: bulk_(bulk),
  solnOpts_(solnOpts),
  timeIntegrator_(timeIntegrator),
  neededData_(in),
  outputData_(out)
{}

bool UtilityKernelManager::add_part(const stk::mesh::Part& part) {
  bool partWasAdded = false;
  if (std::find(parts_.begin(), parts_.end(), &part) == parts_.end()) {
    parts_.push_back(&part);
    partWasAdded = true;
  }
  return partWasAdded;
}

std::unique_ptr<UtilityKernelManager> UtilityKernelManager::create(
  stk::topology topo,
  const stk::mesh::BulkData& bulk,
  const SolutionOptions& solnOpts,
  const TimeIntegrator& timeIntegrator,
  ElemDataRequests& in,
  ElemDataRequests& out)
{
  switch (topo.value()) {
    case stk::topology::NODE: return make_unique<NodalUtilityKernelManager>(bulk, solnOpts, timeIntegrator, in, out);
    default: return nullptr;
  }
}


} // namespace nalu
} // namespace Sierra
