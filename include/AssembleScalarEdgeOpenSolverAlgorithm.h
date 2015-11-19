/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleScalarEdgeOpenSolverAlgorithm_h
#define AssembleScalarEdgeOpenSolverAlgorithm_h

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

class AssembleScalarEdgeOpenSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleScalarEdgeOpenSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    ScalarFieldType *scalarQ,
    ScalarFieldType *bcScalarQ,
    VectorFieldType *dqdx,
    ScalarFieldType *diffFluxCoeff);
  virtual ~AssembleScalarEdgeOpenSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  
  ScalarFieldType *scalarQ_;
  ScalarFieldType *bcScalarQ_;
  VectorFieldType *dqdx_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *coordinates_;
  GenericFieldType *openMassFlowRate_;
};

} // namespace nalu
} // namespace Sierra

#endif
