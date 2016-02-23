/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <LinearSolvers.h>
#include <LinearSolver.h>
#include <LinearSolverConfig.h>
#include <NaluEnv.h>
#include <NaluParsing.h>
#include <Simulation.h>

#include <yaml-cpp/yaml.h>

namespace sierra{
namespace nalu{

LinearSolvers::LinearSolvers(Simulation& sim) : sim_(sim) {}
LinearSolvers::~LinearSolvers()
{
  for(SolverMap::const_iterator pos=solvers_.begin(); pos!=solvers_.end(); ++pos)
    delete pos->second;
  for(SolverEpetraConfigMap::const_iterator pos=solverEpetraConfig_.begin(); pos!=solverEpetraConfig_.end(); ++pos)
    delete pos->second;
  for(SolverTpetraConfigMap::const_iterator pos=solverTpetraConfig_.begin(); pos!=solverTpetraConfig_.end(); ++pos)
    delete pos->second;
}

Simulation *LinearSolvers::root() { return &sim_; }
Simulation *LinearSolvers::parent() { return root(); }

void
LinearSolvers::load(const YAML::Node & node)
{
  const YAML::Node * nodes = node.FindValue("linear_solvers");
  if ( nodes )
  {
    for ( size_t inode = 0; inode <  nodes->size(); ++inode )
    {
      const YAML::Node & linear_solver_node = (* nodes)[inode];
      std::string solver_type = "epetra";
      get_if_present_no_default(linear_solver_node, "type", solver_type);
      if (root()->debug())
      {
        NaluEnv::self().naluOutputP0() << "solver_type= " << solver_type << std::endl;
      }
      if (solver_type == "epetra")
      {
        EpetraLinearSolverConfig * linearSolverConfig = new EpetraLinearSolverConfig();
        linearSolverConfig->load(linear_solver_node);
        solverEpetraConfig_[linearSolverConfig->name()] = linearSolverConfig;
      }
      else if (solver_type == "tpetra")
      {
        TpetraLinearSolverConfig * linearSolverConfig = new TpetraLinearSolverConfig();
        linearSolverConfig->load(linear_solver_node);
        solverTpetraConfig_[linearSolverConfig->name()] = linearSolverConfig;
      }
      else
      {
        throw std::runtime_error("unknown solver type");
      }
    }
  }
}

LinearSolver *
LinearSolvers::create_solver(
  std::string solverBlockName,
  EquationType theEQ )
{

  // provide unique name
  std::string solverName = EquationTypeMap[theEQ] + "_Solver";
  
  LinearSolver *theSolver = NULL;
  // check in epetra map...
  bool foundE = false;
  SolverEpetraConfigMap::const_iterator iterE
    = solverEpetraConfig_.find(solverBlockName);
  if (iterE != solverEpetraConfig_.end()) {
    EpetraLinearSolverConfig *linearSolverConfig = (*iterE).second;
    foundE = true;
    theSolver = new EpetraLinearSolver(solverName,
                                       linearSolverConfig,
                                       linearSolverConfig->aztec_options(),
                                       linearSolverConfig->aztec_parameters(),
                                       linearSolverConfig->use_ml(),
                                       linearSolverConfig->use_mueLu(),
                                       linearSolverConfig->ml_parameters(), this);
  }
  
  // check in tpetra map...
  bool foundT = false;
  SolverTpetraConfigMap::const_iterator iterT
    = solverTpetraConfig_.find(solverBlockName);
  if (iterT != solverTpetraConfig_.end()) {
    TpetraLinearSolverConfig *linearSolverConfig = (*iterT).second;
    foundT = true;
    theSolver = new TpetraLinearSolver(solverName,
                                       linearSolverConfig,
                                       linearSolverConfig->params(),
                                       linearSolverConfig->paramsPrecond(), this);
  }
  
  // error check; both found
  if ( foundE && foundT ) {
    throw std::runtime_error("solver name block duplicated in t and e petra; check: " + solverName);
  }
  // error check; none found
  if ( !foundE && !foundT ) {
    throw std::runtime_error("solver name block not found; error in solver creation; check: " + solverName);
  }

  // set it and return
  solvers_[theEQ] = theSolver;
  return theSolver;
}


} // namespace nalu
} // namespace Sierra
