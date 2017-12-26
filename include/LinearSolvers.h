/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef LinearSolvers_h
#define LinearSolvers_h

#include <Enums.h>

#include <map>
#include <string>

namespace YAML {
class Node;
}

namespace sierra{
namespace nalu{

class LinearSolver;
class TpetraLinearSolverConfig;
class HypreLinearSolverConfig;
class Simulation;

/** Collection of solvers and their associated configuration
 *
 *  This class performs the following actions within a Nalu simulation:
 *
 *  - Parse the `linear_solvers` section and create a mapping of user-defined
 *    configurations.
 *  - Create solvers for specific equation system and update the mapping
 *
 */
class LinearSolvers {
public:
  LinearSolvers(Simulation& sim);
  ~LinearSolvers();
  
  /** Parse the `linear_solvers` section from Nalu input file
   */
  void load(const YAML::Node & node);

  /** Create a solver for the EquationSystem
   *
   *  @param[in] solverBlockName The name specified in the input file, e.g., solve_scalar
   *  @param[in] theEQ The type of equation
   */
  LinearSolver *create_solver(
    std::string solverBlockName,
    EquationType theEQ);
  
  Simulation *root();
  Simulation *parent();
  
  typedef std::map<EquationType, LinearSolver *> SolverMap;
  typedef std::map<std::string, TpetraLinearSolverConfig *> SolverTpetraConfigMap;
  typedef std::map<std::string, HypreLinearSolverConfig*> HypreSolverConfigMap;

  //! Mapping of solver instances to the EquationType
  SolverMap solvers_;

  //! A lookup table of solver configurations against the names provided in the
  //! input file when the `type` is `tpetra`
  SolverTpetraConfigMap solverTpetraConfig_;

  //! A lookup table of solver configurations against the names provided in the
  //! input file when `type` is `hypre` or `tpetra_hypre`
  HypreSolverConfigMap solverHypreConfig_;
  
  //! Reference to the sierra::nalu::Simulation instance
  Simulation& sim_;

private:
};

} // namespace nalu
} // namespace Sierra

#endif
