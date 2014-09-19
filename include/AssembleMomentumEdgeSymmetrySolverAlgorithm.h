/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleMomentumEdgeSymmetrySolverAlgorithm_h
#define AssembleMomentumEdgeSymmetrySolverAlgorithm_h

#include<SolverAlgorithm.h>
#include <FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleMomentumEdgeSymmetrySolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleMomentumEdgeSymmetrySolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem);
  virtual ~AssembleMomentumEdgeSymmetrySolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const double includeDivU_;

  // extract fields
  VectorFieldType *velocity_;
  GenericFieldType *dudx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *viscosity_;
  GenericFieldType *exposedAreaVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
