/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleOversetSolverConstraintAlgorithm_h
#define AssembleOversetSolverConstraintAlgorithm_h

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

class AssembleOversetSolverConstraintAlgorithm : public SolverAlgorithm
{
public:

  AssembleOversetSolverConstraintAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    stk::mesh::FieldBase *fieldQ);
  virtual ~AssembleOversetSolverConstraintAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();
  virtual void prepare_constraints();

  // interface assumes that the correct state was provided
  stk::mesh::FieldBase *fieldQ_;
  
  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
