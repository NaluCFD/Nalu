/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef EquationSystems_h
#define EquationSystems_h

#include <Enums.h>

// yaml for parsing..
#include <yaml-cpp/yaml.h>
#include<NaluParsing.h>

// stk
namespace stk{
namespace mesh{
class Part;
}
}

#include <map>
#include <string>
#include <vector>

namespace YAML {
  class Node;
}

#include<vector>
#include<string>

namespace sierra{
namespace nalu{

class Realm;
class EquationSystem;
class PostProcessingData;
class Simulation;
class AlgorithmDriver;

typedef std::vector<EquationSystem *> EquationSystemVector;

/** A collection of Equations to be solved on a Realm
 *
 *  EquationSystems holds a vector of EquationSystem instances representing the
 *  equations that are being solved in a given Realm and is responsible for the
 *  management of the solve and update of the various field quantities in a
 *  given timestep.
 *
 *  \sa EquationSystems::solve_and_update
 */
class EquationSystems
{
 public:

  EquationSystems(
    Realm &realm);
  ~EquationSystems();

  void load(const YAML::Node & node);
  
  std::string get_solver_block_name(
    const std::string eqName);

  void breadboard();

  Simulation *root();
  Realm *parent();

  // ease of access methods to particular equation system
  size_t size() {return equationSystemVector_.size();}
  EquationSystem *operator[](int i) { return equationSystemVector_[i];}
  
  void register_nodal_fields(
    const std::vector<std::string> targetNames);

  void register_edge_fields(
    const std::vector<std::string> targetNames);

  void register_element_fields(
    const std::vector<std::string> targetNames);

  void register_interior_algorithm(
    const std::vector<std::string> targetNames);

  void register_wall_bc(
    const std::string targetName,
    const WallBoundaryConditionData &wallBCData);

  void register_inflow_bc(
    const std::string targetName,
    const InflowBoundaryConditionData &inflowBCData);

  void register_open_bc(
    const std::string targetName,
    const OpenBoundaryConditionData &openBCData);

  void register_symmetry_bc(
    const std::string targetName,
    const SymmetryBoundaryConditionData &symmetryBCData);

  void register_periodic_bc(
    const std::string targetNameMaster,
    const std::string targetNameSlave,
    const PeriodicBoundaryConditionData &periodicBCData);

  void register_overset_bc(
    const OversetBoundaryConditionData &oversetBCData);

  void register_non_conformal_bc(
    const NonConformalBoundaryConditionData &nonConformalBCData);

  void register_surface_pp_algorithm(
    const PostProcessingData &theData);

  void register_initial_condition_fcn(
    stk::mesh::Part *part,
    const UserFunctionInitialConditionData &fcnIC);

  void initialize();
  void reinitialize_linear_system();
  void populate_derived_quantities();
  void initial_work();

  /** Solve and update the state of all variables for a given timestep
   *
   *  This method is responsible for executing setup actions before calling
   *  solve, performing the actual solve, updating the solution, and performing
   *  post-solve actions after the solution has been updated. To provide
   *  sufficient granularity and control of this pre- and post- solve actions,
   *  the solve method uses the following series of steps:
   *
   *  ```
   *  // Perform tasks for this timestep before any Equation system is called
   *  pre_iter_work();
   *  // Iterate over all equation systems
   *  for (auto eqsys: equationSystems_) {
   *    eqsys->pre_iter_work();
   *    eqsys->solve_and_update();
   *    eqsys->post_iter_work();
   *  }
   *  // Perform tasks after all equation systems have updated
   *  post_iter_work();
   *  ```
   *
   *  Tasks that require to be performed before any equation system is solved
   *  for needs to be registered to preIterAlgDriver_ on EquationSystems,
   *  similiary for post-solve tasks. And actions to be performed immediately
   *  before individual equation system solves need to be registered in
   *  EquationSystem::preIterAlgDriver_.
   *
   *  \sa pre_iter_work(), post_iter_work(), EquationSystem::pre_iter_work(),
   *  \sa EquationSystem::post_iter_work()
   */
  bool solve_and_update();
  double provide_system_norm();
  double provide_mean_system_norm();

  void predict_state();
  void populate_boundary_data();
  void boundary_data_to_state_data();
  void provide_output();
  void dump_eq_time();
  void pre_timestep_work();
  void post_converged_work();
  void evaluate_properties();

  /** Perform necessary setup tasks that affect all EquationSystem instances at
   *  a given timestep.
   *
   *  \sa EquationSystems::solve_and_update()
   */
  void pre_iter_work();

  /** Perform necessary actions once all EquationSystem instances have been
   * updated for the prescribed number of _outer iterations_ at a given
   * timestep.
   *
   *  \sa EquationSystems::solve_and_update()
   */
  void post_iter_work();
  
  Realm &realm_;
  std::string name_;
  int maxIterations_;

  EquationSystemVector equationSystemVector_;
  std::map<std::string, std::string> solverSpecMap_;

  /// A list of tasks to be performed before all EquationSystem::solve_and_update
  std::vector<AlgorithmDriver*> preIterAlgDriver_;

  /// A list of tasks to be performed after all EquationSystem::solve_and_update
  std::vector<AlgorithmDriver*> postIterAlgDriver_;
};

} // namespace nalu
} // namespace Sierra

#endif
