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

#include <element_promotion/ElementDescription.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>

#include <BuildTemplates.h>

#include <algorithm>
#include <tuple>

namespace sierra{
namespace nalu{
  class Realm;


  template <template <typename> class T, int order, typename... Args>
  Kernel* build_ho_kernel(int dimension, Args&&... args)
  {
    // only two topologies supported, so we can flatten some of the decision making
    if (dimension == 2) {
      return new T<AlgTraitsQuadGL<order>>(std::forward<Args>(args)...);
    }
    return new T<AlgTraitsHexGL<order>>(std::forward<Args>(args)...);
  }

  template <template <typename> class T, typename... Args>
  Kernel* build_fem_kernel(stk::topology topo, Args&&... args)
  {
    ThrowRequireMsg(topo == stk::topology::HEXAHEDRON_8, "FEM kernels only implemented for Hex8 topology");
    return new T<AlgTraitsHex8>(std::forward<Args>(args)...);
  }

  template <template <typename> class T, typename... Args>
  Kernel* build_topo_kernel(int dimension, stk::topology topo, Args&&... args)
  {
    if (!topo.is_super_topology()) {
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
    else {
      int poly_order = poly_order_from_super_topology(dimension, topo);
      switch (poly_order) {
        case 2: return build_ho_kernel<T, 2>(dimension, std::forward<Args>(args)...);
        case 3: return build_ho_kernel<T, 3>(dimension, std::forward<Args>(args)...);
        case 4: return build_ho_kernel<T, 4>(dimension, std::forward<Args>(args)...);
        case USER_POLY_ORDER: return build_ho_kernel<T, USER_POLY_ORDER>(dimension, std::forward<Args>(args)...);
        default:
          ThrowRequireMsg(false,
            "Polynomial order" + std::to_string(poly_order) + "is not supported by default.  "
            "Specify USER_POLY_ORDER and recompile to run.");
          return nullptr;
      }
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
    // dimension, in addition to topology, is necessary to distinguish the HO elements,
    const int dim = eqSys.realm_.spatialDimension_;

    bool isCreated = false;
    KernelBuilderLog::self().add_valid_name(eqSys.eqnTypeName_,  name);
    if (eqSys.supp_alg_is_requested(name)) {
      Kernel* compKernel = build_topo_kernel<T>(dim, topo, std::forward<Args>(args)...);
      ThrowRequire(compKernel != nullptr);
      KernelBuilderLog::self().add_built_name(eqSys.eqnTypeName_,  name);
      kernelVec.push_back(compKernel);
      isCreated = true;
    }
    return isCreated;
  }

  template <template <typename> class T, typename... Args>
  bool build_fem_kernel_if_requested(
    stk::topology topo,
    EquationSystem& eqSys,
    std::vector<Kernel*>& kernelVec,
    std::string name,
    Args&&... args)
  {
    bool isCreated = false;
    KernelBuilderLog::self().add_valid_name(eqSys.eqnTypeName_,  name);
    if (eqSys.supp_alg_is_requested(name)) {
      Kernel* compKernel = build_fem_kernel<T>(topo, std::forward<Args>(args)...);
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

    bool isNotNGP = !(topo == stk::topology::HEXAHEDRON_8 ||
                      topo == stk::topology::HEXAHEDRON_27 ||
                      topo == stk::topology::QUADRILATERAL_4_2D ||
                      topo == stk::topology::TRIANGLE_3_2D ||
                      topo == stk::topology::WEDGE_6 ||
                      topo == stk::topology::TETRAHEDRON_4 ||
                      topo == stk::topology::PYRAMID_5);

    auto itc = solverAlgs.find(algName);
    bool createNewAlg = itc == solverAlgs.end();
    if (createNewAlg) {
      auto* theSolverAlg = new AssembleElemSolverAlgorithm(eqSys.realm_, &part, &eqSys, stk::topology::ELEMENT_RANK, topo.num_nodes(), isNotNGP);
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
