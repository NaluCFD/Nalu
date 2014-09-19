/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/

#ifndef AssembleScalarElemDiffSolverAlgorithm_h
#define AssembleScalarElemDiffSolverAlgorithm_h

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

class AssembleScalarElemDiffSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleScalarElemDiffSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    VectorFieldType *dqdx,
    ScalarFieldType *diffFluxCoeff,
    bool useCollcation);
  virtual ~AssembleScalarElemDiffSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

private:

  ScalarFieldType *scalarQ_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *coordinates_;

  bool useCollocation_;

};

} // namespace nalu
} // namespace Sierra

#endif
