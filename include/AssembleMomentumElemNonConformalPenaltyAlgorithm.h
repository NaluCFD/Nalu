/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleMomentumElemNonConformalPenaltyAlgorithm_h
#define AssembleMomentumElemNonConformalPenaltyAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class AssembleMomentumElemNonConformalPenaltyAlgorithm : public Algorithm
{
public:

  AssembleMomentumElemNonConformalPenaltyAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    VectorFieldType *velocity,
    VectorFieldType *ncNormalFlux,
    ScalarFieldType *ncPenalty,
    ScalarFieldType *ncArea,
    ScalarFieldType *diffFluxCoeff);
  ~AssembleMomentumElemNonConformalPenaltyAlgorithm();

  void execute();

  const double includeDivU_;

  VectorFieldType *velocity_;
  VectorFieldType *ncNormalFlux_;
  ScalarFieldType *ncPenalty_;
  ScalarFieldType *ncArea_;
  ScalarFieldType *diffFluxCoeff_;

  VectorFieldType *coordinates_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *massFlowRate_;
};

} // namespace nalu
} // namespace Sierra

#endif
