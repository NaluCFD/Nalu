/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SolverAlgorithm_h
#define SolverAlgorithm_h

#include <Algorithm.h>

#include <stk_mesh/base/Entity.hpp>
#include <vector>

namespace sierra{
namespace nalu{

class EquationSystem;
class Realm;

class SolverAlgorithm : public Algorithm
{
public:

  SolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem);
  virtual ~SolverAlgorithm() {}

  virtual void execute() = 0;
  virtual void initialize_connectivity() = 0;

protected:

  // Need to find out whether this ever gets called inside a modification cycle.
  void apply_coeff(
    const std::vector<stk::mesh::Entity> & sym_meshobj,
    std::vector<int> &scratchIds,
    std::vector<double> &scratchVals,
    const std::vector<double> &rhs,
    const std::vector<double> &lhs,
    const char *trace_tag=0);
  
  EquationSystem *eqSystem_;
};

} // namespace nalu
} // namespace Sierra

#endif
