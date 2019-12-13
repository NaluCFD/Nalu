/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleMomentumElemWallFunctionProjectedSolverAlgorithm_h
#define AssembleMomentumElemWallFunctionProjectedSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
  class Part;
  class Ghosting;
}
}

namespace sierra{
namespace nalu{

class Realm;
class PointInfo;

class AssembleMomentumElemWallFunctionProjectedSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleMomentumElemWallFunctionProjectedSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    const bool &useShifted,
    std::vector<std::vector<PointInfo *> > &pointInfoVec,
    stk::mesh::Ghosting *wallFunctionGhosting);
  virtual ~AssembleMomentumElemWallFunctionProjectedSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const bool useShifted_;
  std::vector<std::vector<PointInfo *> > &pointInfoVec_;
  stk::mesh::Ghosting *wallFunctionGhosting_;

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
  
  // data structure to parallel communicate nodal data to ghosted elements
  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
