/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMdotVofElemAlgorithm_h
#define ComputeMdotVofElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;
class SolutionOptions;

class ComputeMdotVofElemAlgorithm : public Algorithm
{
public:

  ComputeMdotVofElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const SolutionOptions &solnOpts);
  ~ComputeMdotVofElemAlgorithm();

  void execute();

  const bool meshMotion_;

  // extract fields; nodal
  VectorFieldType *velocityRTM_;
  VectorFieldType *Gpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  ScalarFieldType *interfaceCurvature_;
  ScalarFieldType *surfaceTension_;
  ScalarFieldType *vof_;
  GenericFieldType *massFlowRate_;
  GenericFieldType *volumeFlowRate_;
  std::array<double, 3> gravity_;

  const bool shiftMdot_;
  const bool shiftPoisson_;
};

} // namespace nalu
} // namespace Sierra

#endif
