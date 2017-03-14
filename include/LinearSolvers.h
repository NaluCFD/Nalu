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
class Simulation;

class LinearSolvers {
public:
  LinearSolvers(Simulation& sim);
  ~LinearSolvers();
  
  void load(const YAML::Node & node);
  LinearSolver *create_solver(
    std::string solverBlockName,
    EquationType theEQ);
  
  Simulation *root();
  Simulation *parent();
  
  typedef std::map<EquationType, LinearSolver *> SolverMap;
  typedef std::map<std::string, TpetraLinearSolverConfig *> SolverTpetraConfigMap;

  SolverMap solvers_;
  SolverTpetraConfigMap solverTpetraConfig_;
  
  Simulation& sim_;

private:
};

} // namespace nalu
} // namespace Sierra

#endif
