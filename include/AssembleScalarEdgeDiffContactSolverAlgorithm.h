/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleScalarEdgeDiffContactSolverAlgorithm_h
#define AssembleScalarEdgeDiffContactSolverAlgorithm_h

#include <SolverAlgorithm.h>
#include <FieldTypeDef.h>
#include <HermitePolynomialInterpolation.h>

namespace stk {
namespace mesh {
class Part;
class FieldBase;
}
}

namespace sierra{
namespace nalu{

class Realm;

class AssembleScalarEdgeDiffContactSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleScalarEdgeDiffContactSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    VectorFieldType *dqdx,
    ScalarFieldType *diffFluxCoeff,
    const bool useHermiteInterpolation);
  virtual ~AssembleScalarEdgeDiffContactSolverAlgorithm();
  virtual void initialize_connectivity();
  virtual void execute();

  HermitePolynomialInterpolation *hermite_;
  ScalarFieldType *scalarQ_;
  VectorFieldType *dqdx_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *coordinates_;
  ScalarFieldType *haloMdot_;
  
  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;

};

} // namespace nalu
} // namespace Sierra

#endif
