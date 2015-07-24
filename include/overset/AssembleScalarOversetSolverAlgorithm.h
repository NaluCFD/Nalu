/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleScalarOversetSolverAlgorithm_h
#define AssembleScalarOversetSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
class FieldBase;
}
}

namespace sierra{
namespace nalu{

class Realm;

class AssembleScalarOversetSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleScalarOversetSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ);
  virtual ~AssembleScalarOversetSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();
  virtual void prepare_constraints();

  ScalarFieldType *scalarQ_;
  
  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
