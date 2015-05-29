/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleContinuityNonConformalPenaltyAlgorithm_h
#define AssembleContinuityNonConformalPenaltyAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class AssembleContinuityNonConformalPenaltyAlgorithm : public Algorithm
{
public:

  AssembleContinuityNonConformalPenaltyAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *ncNormalFlux,
    ScalarFieldType *ncPenalty,
    ScalarFieldType *ncArea,
    const bool useShifted);
  ~AssembleContinuityNonConformalPenaltyAlgorithm();

  void execute();

  ScalarFieldType *ncNormalFlux_;
  ScalarFieldType *ncPenalty_;
  ScalarFieldType *ncArea_;
  const bool useShifted_;
  const bool meshMotion_;

  VectorFieldType *velocityRTM_;
  ScalarFieldType *density_;
  VectorFieldType *coordinates_;
  GenericFieldType *exposedAreaVec_;  
};

} // namespace nalu
} // namespace Sierra

#endif
