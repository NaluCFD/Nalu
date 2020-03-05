/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SurfaceForceAndMomentWallFunctionAlgorithm_h
#define SurfaceForceAndMomentWallFunctionAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class SurfaceForceAndMomentWallFunctionAlgorithm : public Algorithm
{
public:

  SurfaceForceAndMomentWallFunctionAlgorithm(
    Realm &realm,
    stk::mesh::PartVector &partVec,
    const std::string &outputFileName,
    const std::vector<double > &parameters,
    const bool &useShifted,
    ScalarFieldType *assembledArea);
  ~SurfaceForceAndMomentWallFunctionAlgorithm();

  void execute();

  void pre_work();

  void cross_product(
    double *force, double *cross, double *rad);

  const std::string &outputFileName_;
  const std::vector<double > &parameters_;
  const bool useShifted_;
  const double yplusCrit_;
  const double elog_;
  const double kappa_;

  ScalarFieldType *assembledArea_;
  VectorFieldType *coordinates_;
  VectorFieldType *velocity_;
  ScalarFieldType *pressure_;
  VectorFieldType *pressureForce_;
  ScalarFieldType *tauWall_;
  ScalarFieldType *yplus_;
  VectorFieldType *bcVelocity_;
  ScalarFieldType *density_;
  ScalarFieldType *viscosity_;
  GenericFieldType *wallFrictionVelocityBip_;
  GenericFieldType *wallNormalDistanceBip_;
  GenericFieldType *exposedAreaVec_;

  const int w_;
};

} // namespace nalu
} // namespace Sierra

#endif
