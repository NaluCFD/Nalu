/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/
#ifndef SupplementalAlgorithmBuilder_h
#define SupplementalAlgorithmBuilder_h

#include <SupplementalAlgorithm.h>
#include <AssembleElemSolverAlgorithm.h>
#include <EquationSystem.h>
#include <AlgTraits.h>
#include <SupplementalAlgorithmBuilderLog.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>

#include <algorithm>
#include <tuple>

namespace sierra{
namespace nalu{
  class Realm;

  template <class T,  typename... Args>
  void append_new_supp_alg_if_requested(
    std::vector<std::string> algNameVec,
    std::vector<SupplementalAlgorithm*>& algVec,
    std::string suppAlgTypeName,
    Args&&... args)
  {
    std::string name = T::name;

    bool isRequested = std::find(algNameVec.begin(), algNameVec.end(), name) != algNameVec.end();
    if (isRequested) {
      auto* suppAlg = new T(std::forward<Args>(args)...);
      ThrowRequire(suppAlg != nullptr);
      algVec.push_back(suppAlg);
      SuppAlgBuilderLog::self().add_built_name(suppAlgTypeName, name);
    }
    SuppAlgBuilderLog::self().add_valid_name(suppAlgTypeName, name);
  }

  template <template <typename> class T, typename... Args>
  void build_topo_supp_alg_if_requested(
    stk::topology topo,
    const std::map<std::string, std::vector<std::string>>& algNameMap,
    std::vector<SupplementalAlgorithm*>& algVec,
    std::string suppAlgTypeName,
    Args&&... args)
  {
    auto it = algNameMap.find(suppAlgTypeName);
    if (it != algNameMap.end())
    {
      const auto& algNameVec = it->second;

      switch(topo.value()) {
        case stk::topology::HEX_8:
          append_new_supp_alg_if_requested<T<AlgTraitsHex8>>(algNameVec, algVec, suppAlgTypeName,
              std::forward<Args>(args)...);
          break;
        case stk::topology::TET_4:
          append_new_supp_alg_if_requested<T<AlgTraitsTet4>>(algNameVec, algVec, suppAlgTypeName,
              std::forward<Args>(args)...);
          break;
        case stk::topology::PYRAMID_5:
          append_new_supp_alg_if_requested<T<AlgTraitsPyr5>>(algNameVec, algVec, suppAlgTypeName,
              std::forward<Args>(args)...);
          break;
        case stk::topology::WEDGE_6:
          append_new_supp_alg_if_requested<T<AlgTraitsWed6>>(algNameVec, algVec, suppAlgTypeName,
              std::forward<Args>(args)...);
          break;
        default: break;
      }
    }
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

    return {theSolverAlg, createNewAlg} ;
  }

} // namespace nalu
} // namespace Sierra

#endif
