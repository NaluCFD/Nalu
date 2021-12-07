/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMdotVofElemOpenAlgorithm_h
#define ComputeMdotVofElemOpenAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;
class SolutionOptions;

class ComputeMdotVofElemOpenAlgorithm : public Algorithm
{
public:

  ComputeMdotVofElemOpenAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const SolutionOptions &solnOpts);
  ~ComputeMdotVofElemOpenAlgorithm();

  void execute();

  VectorFieldType *velocityRTM_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  ScalarFieldType *interfaceCurvature_;
  ScalarFieldType *surfaceTension_;
  ScalarFieldType *vof_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *openMassFlowRate_;
  GenericFieldType *openVolumeFlowRate_;
  GenericFieldType *dynamicPressure_;
  ScalarFieldType *pressureBc_;

  const bool shiftMdot_;
  const bool shiftedGradOp_;
  const double penaltyFac_;
};

} // namespace nalu
} // namespace Sierra

#endif
