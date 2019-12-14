/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SurfaceForceAndMomentWallFunctionProjectedAlgorithm_h
#define SurfaceForceAndMomentWallFunctionProjectedAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

// forwards
namespace stk {
  namespace mesh {
    class Ghosting;
  }
}

namespace sierra{
namespace nalu{

class Realm;
class PointInfo;

class SurfaceForceAndMomentWallFunctionProjectedAlgorithm : public Algorithm
{
public:

  SurfaceForceAndMomentWallFunctionProjectedAlgorithm(
    Realm &realm,
    stk::mesh::PartVector &partVec,
    const std::string &outputFileName,
    const int &frequency_,
    const std::vector<double > &parameters,
    const bool &useShifted,
    ScalarFieldType *assembledArea,
    std::vector<std::vector<PointInfo *> > &pointInfoVec,
    stk::mesh::Ghosting *wallFunctionGhosting);
  ~SurfaceForceAndMomentWallFunctionProjectedAlgorithm();

  void execute();
  
  void pre_work();

  void cross_product(
    double *force, double *cross, double *rad);

  void error_check();

  const std::string &outputFileName_;
  const int &frequency_;
  const std::vector<double > &parameters_;
  const bool useShifted_;
  std::vector<std::vector<PointInfo *> > &pointInfoVec_;  
  stk::mesh::Ghosting *wallFunctionGhosting_;
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

  // sanity check for like partVec_ and pointInfoVec_ size
  bool errorCheckProcessed_;

  // data structure to parallel communicate nodal data to ghosted elements
  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
