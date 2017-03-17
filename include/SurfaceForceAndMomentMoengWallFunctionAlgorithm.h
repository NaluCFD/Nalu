/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SurfaceForceAndMomentMoengWallFunctionAlgorithm_h
#define SurfaceForceAndMomentMoengWallFunctionAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class SurfaceForceAndMomentMoengWallFunctionAlgorithm : public Algorithm
{
public:

  SurfaceForceAndMomentMoengWallFunctionAlgorithm(
    Realm &realm,
    stk::mesh::PartVector &partVec,
    const std::string &outputFileName,
    const int &frequency_,
    const std::vector<double > &parameters,
    const bool &useShifted);
  ~SurfaceForceAndMomentMoengWallFunctionAlgorithm();

  void execute();

  void pre_work();

  void cross_product(
    double *force, double *cross, double *rad);

  const std::string &outputFileName_;
  const int &frequency_;
  const std::vector<double > &parameters_;
  const bool useShifted_;

  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  VectorFieldType *pressureForce_;
  ScalarFieldType *tauWall_;
  ScalarFieldType *yplus_;
  ScalarFieldType *density_;
  ScalarFieldType *viscosity_;
  GenericFieldType *tauWallBip_;
  GenericFieldType *wallNormalDistanceBip_;
  GenericFieldType *exposedAreaVec_;
  ScalarFieldType *assembledArea_;

  const int w_;
};

} // namespace nalu
} // namespace Sierra

#endif
