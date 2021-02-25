/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SixDofSurfaceForceAndMomentAlgorithm_h
#define SixDofSurfaceForceAndMomentAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>
#include<MeshMotionInfo.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class SixDofSurfaceForceAndMomentAlgorithm : public Algorithm
{
public:

  SixDofSurfaceForceAndMomentAlgorithm(
    Realm &realm,
    MeshMotionInfo *motion, 
    stk::mesh::PartVector &partVec,
    const bool &useShifted,
    ScalarFieldType *assembledArea);
  ~SixDofSurfaceForceAndMomentAlgorithm();

  void execute();

  void pre_work();

  void cross_product(
    double *force, double *cross, double *rad);

  const bool useShifted_;
  const double includeDivU_;

  MeshMotionInfo *motion_; 
  ScalarFieldType *assembledArea_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  ScalarFieldType *viscosity_;
  GenericFieldType *dudx_;
  GenericFieldType *exposedAreaVec_;

  const int w_;

};

} // namespace nalu
} // namespace Sierra

#endif
