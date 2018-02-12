/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeTurbKineticEnergyWallFunctionAlgorithm_h
#define ComputeTurbKineticEnergyWallFunctionAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeTurbKineticEnergyWallFunctionAlgorithm : public Algorithm
{
public:

  ComputeTurbKineticEnergyWallFunctionAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);

  virtual ~ComputeTurbKineticEnergyWallFunctionAlgorithm();

  void execute();

  void zero_nodal_fields();

  void normalize_nodal_fields();

  const double cMu_;

  ScalarFieldType *turbKineticEnergy_;
  ScalarFieldType *bcTurbKineticEnergy_;
  ScalarFieldType *bcAssembledTurbKineticEnergy_;
  ScalarFieldType *assembledWallArea_;
  GenericFieldType *wallFrictionVelocityBip_;
  GenericFieldType *exposedAreaVec_;

};

} // namespace nalu
} // namespace Sierra

#endif
