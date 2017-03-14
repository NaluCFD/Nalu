/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleMomentumNonConformalSolverAlgorithm_h
#define AssembleMomentumNonConformalSolverAlgorithm_h

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

class AssembleMomentumNonConformalSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleMomentumNonConformalSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    VectorFieldType *velocity,
    ScalarFieldType *diffFluxCoeff);
  virtual ~AssembleMomentumNonConformalSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  VectorFieldType *velocity_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *coordinates_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *ncMassFlowRate_;

  // options that prevail over all algorithms created
  const double eta_;
  const double includeDivU_;
  const double useCurrentNormal_;

  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
