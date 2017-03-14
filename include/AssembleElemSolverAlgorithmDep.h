/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleElemSolverAlgorithmDep_h
#define AssembleElemSolverAlgorithmDep_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;

class AssembleElemSolverAlgorithmDep : public SolverAlgorithm
{
public:

  AssembleElemSolverAlgorithmDep(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem);
  virtual ~AssembleElemSolverAlgorithmDep() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const int sizeOfSystem_;
};

} // namespace nalu
} // namespace Sierra

#endif
