/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleMomentumMoengWallFunctionSolverAlgorithm_h
#define AssembleMomentumMoengWallFunctionSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>
#include<ABLProfileFunction.h>

namespace stk {
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;

class AssembleMomentumMoengWallFunctionSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleMomentumMoengWallFunctionSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    const bool &useShifted,
    const double &gravity,
    const double &z0,
    const double &Tref);
  virtual ~AssembleMomentumMoengWallFunctionSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();
  void compute_utau(
    const double &up, const double &zp, 
    const double &qsurf, const ABLProfileFunction *ABLProfFun, 
    double &utau );

  const bool useShifted_;
  const double z0_;
  const double Tref_;
  const double gravity_;
  const double alpha_h_;
  const double beta_m_;
  const double beta_h_;
  const double gamma_m_;
  const double gamma_h_;
  const double kappa_;

  VectorFieldType *velocity_;
  VectorFieldType *bcVelocity_;
  VectorFieldType *coordinates_;
  ScalarFieldType *bcHeatFlux_;
  ScalarFieldType *density_;
  ScalarFieldType *specificHeat_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *wallNormalDistanceBip_;
  GenericFieldType *tauWallBip_;
};

} // namespace nalu
} // namespace Sierra

#endif
