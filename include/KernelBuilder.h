/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef KernelBuilder_h
#define KernelBuilder_h

#include <Kernel.h>
#include <AssembleElemSolverAlgorithm.h>
#include <EquationSystem.h>
#include <AlgTraits.h>
#include <KernelBuilderLog.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>

#include <algorithm>
#include <tuple>

namespace sierra{
namespace nalu{
  class Realm;

  template <template <typename> class T, typename... Args>
  Kernel* build_topo_kernel(stk::topology topo, Args&&... args)
  {
    switch(topo.value()) {
      case stk::topology::HEX_8:
        return new T<AlgTraitsHex8>(std::forward<Args>(args)...);
      case stk::topology::HEX_27:
        return new T<AlgTraitsHex27>(std::forward<Args>(args)...);
      case stk::topology::TET_4:
        return new T<AlgTraitsTet4>(std::forward<Args>(args)...);
      case stk::topology::PYRAMID_5:
        return new T<AlgTraitsPyr5>(std::forward<Args>(args)...);
      case stk::topology::WEDGE_6:
        return new T<AlgTraitsWed6>(std::forward<Args>(args)...);
      case stk::topology::QUAD_4_2D:
        return new T<AlgTraitsQuad4_2D>(std::forward<Args>(args)...);
      case stk::topology::QUAD_9_2D:
        return new T<AlgTraitsQuad9_2D>(std::forward<Args>(args)...);
      case stk::topology::TRI_3_2D:
        return new T<AlgTraitsTri3_2D>(std::forward<Args>(args)...);
      default:
        return nullptr;
    }
  }

  template <template <typename> class T, typename... Args>
  bool build_topo_kernel_if_requested(
    stk::topology topo,
    EquationSystem& eqSys,
    std::vector<Kernel*>& kernelVec,
    std::string name,
    Args&&... args)
  {
    bool isCreated = false;
    KernelBuilderLog::self().add_valid_name(eqSys.eqnTypeName_,  name);
    if (eqSys.supp_alg_is_requested(name)) {
      Kernel* compKernel = build_topo_kernel<T>(topo, std::forward<Args>(args)...);
      ThrowRequire(compKernel != nullptr);
      KernelBuilderLog::self().add_built_name(eqSys.eqnTypeName_,  name);
      kernelVec.push_back(compKernel);
      isCreated = true;
    }
    return isCreated;
  }

  inline std::pair<AssembleElemSolverAlgorithm*, bool>
  build_or_add_part_to_solver_alg(
    EquationSystem& eqSys,
    stk::mesh::Part& part,
    std::map<std::string, SolverAlgorithm*>& solverAlgs)
  {
    const stk::topology topo = part.topology();
    const std::string algName = "AssembleElemSolverAlg_" + topo.name();

    auto itc = solverAlgs.find(algName);
    bool createNewAlg = itc == solverAlgs.end();
    if (createNewAlg) {
      auto* theSolverAlg = new AssembleElemSolverAlgorithm(eqSys.realm_, &part, &eqSys, topo);
      ThrowRequire(theSolverAlg != nullptr);

      NaluEnv::self().naluOutputP0() << "Created the following alg: " << algName << std::endl;
      solverAlgs.insert({algName, theSolverAlg});
    }
    else {
      auto& partVec = itc->second->partVec_;
      if (std::find(partVec.begin(), partVec.end(), &part) == partVec.end()) {
        partVec.push_back(&part);
      }
    }

    auto* theSolverAlg = dynamic_cast<AssembleElemSolverAlgorithm*>(solverAlgs.at(algName));
    ThrowRequire(theSolverAlg != nullptr);

    return {theSolverAlg, createNewAlg};
  }

} // namespace nalu
} // namespace Sierra

#endif
