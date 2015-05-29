/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleScalarDiffNonConformalSolverAlgorithm_h
#define AssembleScalarDiffNonConformalSolverAlgorithm_h

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

class AssembleScalarDiffNonConformalSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleScalarDiffNonConformalSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    ScalarFieldType *diffFluxCoeff);
  virtual ~AssembleScalarDiffNonConformalSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  ScalarFieldType *scalarQ_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *coordinates_;
  GenericFieldType *exposedAreaVec_;

  // options that prevail over all algorithms created
  bool robinStyle_;
  double dsFactor_;

  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;

};

} // namespace nalu
} // namespace Sierra

#endif
