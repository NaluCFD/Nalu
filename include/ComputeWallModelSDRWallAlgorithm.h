/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeWallModelSDRWallAlgorithm_h
#define ComputeWallModelSDRWallAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeWallModelSDRWallAlgorithm : public Algorithm
{
public:

  ComputeWallModelSDRWallAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const bool &useShifted);
  virtual ~ComputeWallModelSDRWallAlgorithm();

  void execute();

  const bool useShifted_;
  const double sqrtBetaStar_;
  const double kappa_;

  VectorFieldType *coordinates_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *wallFrictionVelocityBip_;
  ScalarFieldType *sdrBc_;
  ScalarFieldType *assembledWallArea_;
};

} // namespace nalu
} // namespace Sierra

#endif
