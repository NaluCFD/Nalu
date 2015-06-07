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
    EquationSystem *eqSystem);
  virtual ~AssembleElemSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const int sizeOfSystem_;
};

} // namespace nalu
} // namespace Sierra

#endif
