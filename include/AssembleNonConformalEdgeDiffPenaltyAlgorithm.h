/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNonConformalEdgeDiffPenaltyAlgorithm_h
#define AssembleNonConformalEdgeDiffPenaltyAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNonConformalEdgeDiffPenaltyAlgorithm : public Algorithm
{
public:

  AssembleNonConformalEdgeDiffPenaltyAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *scalarQ,
    VectorFieldType *GjQ,
    ScalarFieldType *ncNormalFlux,
    ScalarFieldType *ncPenalty,
    ScalarFieldType *ncArea,
    ScalarFieldType *diffFluxCoeff);
  ~AssembleNonConformalEdgeDiffPenaltyAlgorithm();

  void execute();
  void zero_fields();
  void assemble_and_normalize();

  ScalarFieldType *scalarQ_;
  VectorFieldType *GjQ_;
  ScalarFieldType *ncNormalFlux_;
  ScalarFieldType *ncPenalty_;
  ScalarFieldType *ncArea_;
  ScalarFieldType *diffFluxCoeff_;

  VectorFieldType *coordinates_;
  GenericFieldType *exposedAreaVec_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
