/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleScalarElemOpenSolverAlgorithm_h
#define AssembleScalarElemOpenSolverAlgorithm_h

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

class AssembleScalarElemOpenSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleScalarElemOpenSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    ScalarFieldType *bcScalarQ,
    VectorFieldType *dqdx,
    ScalarFieldType *diffFluxCoeff);
  virtual ~AssembleScalarElemOpenSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  
  ScalarFieldType *scalarQ_;
  ScalarFieldType *bcScalarQ_;
  VectorFieldType *dqdx_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *velocity_;
  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  GenericFieldType *openMassFlowRate_;

};

} // namespace nalu
} // namespace Sierra

#endif
