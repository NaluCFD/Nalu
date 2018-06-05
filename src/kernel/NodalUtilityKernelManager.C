/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
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

NodalUtilityKernelManager::NodalUtilityKernelManager(
  const stk::mesh::BulkData& bulk,
  const SolutionOptions& solnOpts,
  const TimeIntegrator& timeIntegrator,
  ElemDataRequests& in,
  ElemDataRequests& out)
: UtilityKernelManager(bulk, solnOpts, timeIntegrator, in, out)
{ }

void
NodalUtilityKernelManager::execute_kernels(const stk::mesh::Selector& selector)
{
  for (auto& kernelPair : kernels_) {
    kernelPair.second->pre_execute(timeIntegrator_, bulk_.parallel());
  }

  const int bytes_per_thread = bytes_needed_for_nodal_utility_kernels(neededData_, outputData_);
  const int bytes_per_team = 0;

  const auto& buckets = bulk_.get_buckets(stk::topology::NODE_RANK, selector);
  auto team_exec = get_team_policy(buckets.size(), bytes_per_team, bytes_per_thread);
  Kokkos::parallel_for(team_exec, [&](const sierra::nalu::TeamHandleType& team) {
    stk::mesh::Bucket& b = *buckets[team.league_rank()];

    const stk::mesh::Bucket::size_type bucketLen = b.size();
    const size_t simdBucketLen = get_num_simd_groups(bucketLen);
    std::array<stk::mesh::Entity, simdLen> simdNodes = {{}};
    NodalScratchData<DoubleType> simdInputData(team, bulk_, neededData_);
    NodalScratchData<DoubleType> simdOutputData(team, bulk_, outputData_);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, simdBucketLen), [&](const size_t k) {
      const int numSimdNodes = get_length_of_next_simd_group(k, bucketLen);

      for (int simdIndex = 0; simdIndex < numSimdNodes; ++simdIndex) {
        simdNodes[simdIndex] = b[k*simdLen + simdIndex];
      }

      gather_data(numSimdNodes, simdNodes, simdInputData);
      for (auto& kernelPair : kernels_) kernelPair.second->execute(simdInputData, simdOutputData);
      scatter_data(numSimdNodes, simdNodes, simdOutputData);
    });
  });

  for (auto& kernelPair : kernels_) {
    kernelPair.second->post_execute(timeIntegrator_, bulk_.parallel());
  }
}

void
NodalUtilityKernelManager::gather_data(
  int numSimdNodes,
  const std::array<stk::mesh::Entity, stk::simd::ndoubles>& nodes,
  NodalScratchData<DoubleType>& scratchData) const
{
  for (const FieldInfo& fieldInfo : neededData_.get_fields()) {
    unsigned scalarsDim1 = fieldInfo.scalarsDim1;
    bool isTensorField = fieldInfo.scalarsDim2 > 1;
    bool isVectorField = fieldInfo.scalarsDim1 > 1 && !isTensorField;

    const stk::mesh::FieldBase& field = *fieldInfo.field;

    for (int simdIndex = 0; simdIndex < numSimdNodes; ++simdIndex) {
      const auto* field_data = static_cast<const double*>(stk::mesh::field_data(field, nodes[simdIndex]));

      if (!isVectorField && !isTensorField) {
        stk::simd::set_data(scratchData.get_value(field), simdIndex, *field_data);
      }
      else if (isVectorField) {
        const auto& field_view = scratchData.get_1D(field);
        for (unsigned d = 0; d < scalarsDim1; ++d) {
          stk::simd::set_data(field_view(d), simdIndex, field_data[d]);
        }
      }
      else if (isTensorField) {
        unsigned scalarsDim2 = fieldInfo.scalarsDim2;
        const auto& field_view = scratchData.get_2D(field);
        for (unsigned d_outer = 0; d_outer < scalarsDim1; ++d_outer) {
          for (unsigned d_inner = 0; d_inner < scalarsDim2; ++d_inner) {
            stk::simd::set_data(field_view(d_outer, d_inner), simdIndex, field_data[d_outer * scalarsDim2 + d_inner]);
          }
        }
      }
    }
  }
}

void
NodalUtilityKernelManager::scatter_data(
  int numSimdNodes,
  const std::array<stk::mesh::Entity, stk::simd::ndoubles>& nodes,
  const NodalScratchData<DoubleType>& scratchData) const
{
  for (const FieldInfo& fieldInfo : outputData_.get_fields()) {
    const unsigned scalarsDim1 = fieldInfo.scalarsDim1;
    const bool isTensorField = fieldInfo.scalarsDim2 > 1;
    const bool isVectorField = fieldInfo.scalarsDim1 > 1 && !isTensorField;

    const stk::mesh::FieldBase& field = *fieldInfo.field;

    for (int simdIndex = 0; simdIndex < numSimdNodes; ++simdIndex) {
      auto* field_data = static_cast<double*>(stk::mesh::field_data(field, nodes[simdIndex]));

      if (!isVectorField && !isTensorField) {
        *field_data = stk::simd::get_data(scratchData.get_value(field), simdIndex);
      }
      else if (isVectorField) {
        const auto& field_view = scratchData.get_1D(field);
        for (unsigned d = 0; d < scalarsDim1; ++d) {
          field_data[d] = stk::simd::get_data(field_view(d), simdIndex);
        }
      }
      else if (isTensorField) {
        const unsigned scalarsDim2 = fieldInfo.scalarsDim2;
        const auto& field_view = scratchData.get_2D(field);
        for (unsigned d_outer = 0; d_outer < scalarsDim1; ++d_outer) {
          for (unsigned d_inner = 0; d_inner < scalarsDim2; ++d_inner) {
            field_data[d_outer * scalarsDim2 + d_inner] = stk::simd::get_data(field_view(d_outer, d_inner), simdIndex);
          }
        }
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
