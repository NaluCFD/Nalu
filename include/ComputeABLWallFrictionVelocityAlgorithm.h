/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeABLWallFrictionVelocityAlgorithm_h
#define ComputeABLWallFrictionVelocityAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>
#include<ABLProfileFunction.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeABLWallFrictionVelocityAlgorithm : public Algorithm
{
public:

  ComputeABLWallFrictionVelocityAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const bool &useShifted,
    const double &gravity,
    const double &z0,
    const double &Tref);
  virtual ~ComputeABLWallFrictionVelocityAlgorithm();

  void execute();

  void zero_nodal_fields();

  void compute_utau(
      const double &up, const double &zp,
      const double &qsurf, const ABLProfileFunction *ABLProfFun,
      double &utau);
  
  void normalize_nodal_fields();

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
  const int maxIteration_;
  const double tolerance_;

  VectorFieldType *velocity_;
  VectorFieldType *bcVelocity_;
  ScalarFieldType *bcHeatFlux_;
  ScalarFieldType *density_;
  ScalarFieldType *specificHeat_;
  VectorFieldType *coordinates_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *wallFrictionVelocityBip_;
  GenericFieldType *wallNormalDistanceBip_;
  ScalarFieldType *assembledWallNormalDistance_;
  ScalarFieldType *assembledWallArea_;
};

} // namespace nalu
} // namespace Sierra

#endif
