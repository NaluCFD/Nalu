/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleMomentumEdgeOpenSolverAlgorithm_h
#define AssembleMomentumEdgeOpenSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include <FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleMomentumEdgeOpenSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleMomentumEdgeOpenSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem);
  virtual ~AssembleMomentumEdgeOpenSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const double includeDivU_;

  // extract fields
  VectorFieldType *velocity_;
  GenericFieldType *dudx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *viscosity_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *openMassFlowRate_;
  VectorFieldType *velocityBc_;
};

} // namespace nalu
} // namespace Sierra

#endif
