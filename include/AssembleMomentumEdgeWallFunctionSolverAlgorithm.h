/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleMomentumEdgeWallFunctionSolverAlgorithm_h
#define AssembleMomentumEdgeWallFunctionSolverAlgorithm_h

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

class AssembleMomentumEdgeWallFunctionSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleMomentumEdgeWallFunctionSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem);
  virtual ~AssembleMomentumEdgeWallFunctionSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const double yplusCrit_;
  const double elog_;
  const double kappa_;

  VectorFieldType *velocity_;
  VectorFieldType *bcVelocity_;
  ScalarFieldType *density_;
  ScalarFieldType *viscosity_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *wallFrictionVelocityBip_;
  GenericFieldType *wallNormalDistanceBip_;
};

} // namespace nalu
} // namespace Sierra

#endif
