/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef NodalUtilityKernelManager_h
#define NodalUtilityKernelManager_h

#include <kernel/UtilityKernelManager.h>
#include <SolverAlgorithm.h>
#include <ElemDataRequests.h>
#include <FieldTypeDef.h>
#include <nalu_make_unique.h>
#include <kernel/Kernel.h>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class NodalUtilityKernelManager : public UtilityKernelManager
{
public:
  NodalUtilityKernelManager(
    const stk::mesh::BulkData& bulk,
    const SolutionOptions& solnOpts,
    const TimeIntegrator& timeIntegrator,
    ElemDataRequests& neededData,
    ElemDataRequests& outputData);
  ~NodalUtilityKernelManager() = default;

  void execute_kernels(const stk::mesh::Selector& selector) final;
private:
  void gather_data(
    int numSimdNodes,
    const std::array<stk::mesh::Entity, stk::simd::ndoubles>& nodes,
    NodalScratchData<DoubleType>& simdInputData) const;

  void scatter_data(
    int numSimdNodes,
    const std::array<stk::mesh::Entity, stk::simd::ndoubles>& nodes,
    const NodalScratchData<DoubleType>& simdOutputData) const;

};


} // namespace nalu
} // namespace Sierra

#endif

