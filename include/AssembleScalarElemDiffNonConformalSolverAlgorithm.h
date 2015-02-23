/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleScalarElemDiffNonConformalSolverAlgorithm_h
#define AssembleScalarElemDiffNonConformalSolverAlgorithm_h

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

class AssembleScalarElemDiffNonConformalSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleScalarElemDiffNonConformalSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    ScalarFieldType *ncNormalFlux,
    ScalarFieldType *diffFluxCoeff,
    ScalarFieldType *ncPenalty);
  virtual ~AssembleScalarElemDiffNonConformalSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  ScalarFieldType *scalarQ_;
  ScalarFieldType *ncNormalFlux_;
  ScalarFieldType *diffFluxCoeff_;
  ScalarFieldType *ncPenalty_;
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
