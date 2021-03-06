/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeWallModelTurbDissipationWallAlgorithm_h
#define ComputeWallModelTurbDissipationWallAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeWallModelTurbDissipationWallAlgorithm : public Algorithm
{
public:

  ComputeWallModelTurbDissipationWallAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  virtual ~ComputeWallModelTurbDissipationWallAlgorithm();

  void execute();

  const double kappa_;

  GenericFieldType *exposedAreaVec_;
  GenericFieldType *wallFrictionVelocityBip_;
  GenericFieldType *wallNormalDistanceBip_;
  ScalarFieldType *epsBc_;
  ScalarFieldType *assembledWallArea_;
};

} // namespace nalu
} // namespace Sierra

#endif
