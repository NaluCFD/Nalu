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

typedef std::vector<EquationSystem *> EquationSystemVector;

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

  void register_contact_bc(
    const std::string targetName,
    const ContactBoundaryConditionData &contactBCData);

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
  void post_adapt_work();
  void populate_derived_quantities();
  void initial_work();
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
  
  Realm &realm_;
  std::string name_;
  int maxIterations_;

  EquationSystemVector equationSystemVector_;
  std::map<std::string, std::string> solverSpecMap_;
};

} // namespace nalu
} // namespace Sierra

#endif
