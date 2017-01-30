/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleElemSolverAlgorithm_h
#define AssembleElemSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
class Topology;
}
}

namespace sierra{
namespace nalu{

class Realm;

class AssembleElemSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleElemSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    const stk::topology &theTopo);
  virtual ~AssembleElemSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  std::vector<double> lhs_;
  std::vector<double> rhs_;
  std::vector<int> scratchIds_;
  std::vector<double> scratchVals_;
  std::vector<stk::mesh::Entity> connectedNodes_;
};

} // namespace nalu
} // namespace Sierra

#endif
