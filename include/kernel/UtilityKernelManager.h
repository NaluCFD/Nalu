/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef UtilityKernelManager_h
#define UtilityKernelManager_h

#include <SolverAlgorithm.h>
#include <ElemDataRequests.h>
#include <FieldTypeDef.h>
#include <nalu_make_unique.h>
#include <kernel/Kernel.h>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;


class UtilityKernelManager
{
public:
  using UtilityKernelMapType = std::map<std::string, std::unique_ptr<NodalUtilityKernel>>;

  static std::unique_ptr<UtilityKernelManager> create(
    stk::topology topo,
    const stk::mesh::BulkData& bulk,
    const SolutionOptions& solnOpts,
    const TimeIntegrator& timeIntegrator,
    ElemDataRequests& neededData,
    ElemDataRequests& outputData
  );

  UtilityKernelManager(
    const stk::mesh::BulkData& bulk,
    const SolutionOptions& solnOpts,
    const TimeIntegrator& timeIntegrator,
    ElemDataRequests& neededData,
    ElemDataRequests& outputData);

  virtual ~UtilityKernelManager() = default;

  virtual void execute_kernels(const stk::mesh::Selector& selector) {};

  template<typename T, typename... Args>
  std::pair<UtilityKernelMapType::iterator, bool> build_utility_kernel(std::string name, Args&&... args) {
    return kernels_.emplace(name, make_unique<T>(std::forward<Args>(args)...));
  }

  bool add_part(const stk::mesh::Part& part);

protected:
  const stk::mesh::BulkData& bulk_;
  const SolutionOptions& solnOpts_;
  const TimeIntegrator& timeIntegrator_;
  stk::mesh::ConstPartVector parts_;

  ElemDataRequests& neededData_;
  ElemDataRequests& outputData_;

  UtilityKernelMapType kernels_;
};

} // namespace nalu
} // namespace Sierra

#endif

