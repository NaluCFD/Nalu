/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeLowReynoldsSDRWallAlgorithm_h
#define ComputeLowReynoldsSDRWallAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeLowReynoldsSDRWallAlgorithm : public Algorithm
{
public:

  ComputeLowReynoldsSDRWallAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const bool &useShifted);
  virtual ~ComputeLowReynoldsSDRWallAlgorithm();

  void execute();

  const bool useShifted_;
  const double betaOne_;
  const double wallFactor_;

  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  ScalarFieldType *viscosity_;
  GenericFieldType *exposedAreaVec_;
  ScalarFieldType *sdrBc_;
  ScalarFieldType *assembledWallArea_;
};

} // namespace nalu
} // namespace Sierra

#endif
