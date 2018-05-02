/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMdotElemOpenPenaltyAlgorithm_h
#define ComputeMdotElemOpenPenaltyAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeMdotElemOpenPenaltyAlgorithm : public Algorithm
{
public:

  ComputeMdotElemOpenPenaltyAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  ~ComputeMdotElemOpenPenaltyAlgorithm();

  void execute();

  VectorFieldType *velocityRTM_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *openMassFlowRate_;
  ScalarFieldType *pressureBc_;

  const double interpTogether_;
  const double om_interpTogether_;
  const bool shiftMdot_;
  const bool shiftedGradOp_;
  const double stabFac_;
};

} // namespace nalu
} // namespace Sierra

#endif
