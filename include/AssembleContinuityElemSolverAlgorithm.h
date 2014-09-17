/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleContinuityElemSolverAlgorithm_h
#define AssembleContinuityElemSolverAlgorithm_h

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

class AssembleContinuityElemSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleContinuityElemSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem);
  virtual ~AssembleContinuityElemSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  // extract fields; nodal
  VectorFieldType *velocity_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  GenericFieldType *massFlowRate_;
};

} // namespace nalu
} // namespace Sierra

#endif
