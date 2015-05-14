/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleScalarNonConformalSolverAlgorithm_h
#define AssembleScalarNonConformalSolverAlgorithm_h

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

class AssembleScalarNonConformalSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleScalarNonConformalSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    ScalarFieldType *ncNormalFlux,
    ScalarFieldType *ncPenalty,
    const bool normalizeByTimeScale = false);
  virtual ~AssembleScalarNonConformalSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  ScalarFieldType *scalarQ_;
  ScalarFieldType *ncNormalFlux_;
  ScalarFieldType *ncPenalty_;
  const bool normalizeByTimeScale_;
  GenericFieldType *exposedAreaVec_;

  // options that prevail over all algorithms created
  bool robinStyle_;
  double dsFactor_;

  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;

};

} // namespace nalu
} // namespace Sierra

#endif
