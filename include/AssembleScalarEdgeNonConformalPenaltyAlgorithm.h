/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleScalarEdgeNonConformalPenaltyAlgorithm_h
#define AssembleScalarEdgeNonConformalPenaltyAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class AssembleScalarEdgeNonConformalPenaltyAlgorithm : public Algorithm
{
public:

  AssembleScalarEdgeNonConformalPenaltyAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *scalarQ,
    VectorFieldType *GjQ,
    ScalarFieldType *ncNormalFlux,
    ScalarFieldType *ncPenalty,
    ScalarFieldType *ncArea,
    ScalarFieldType *diffFluxCoeff);
  ~AssembleScalarEdgeNonConformalPenaltyAlgorithm();

  void execute();

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
