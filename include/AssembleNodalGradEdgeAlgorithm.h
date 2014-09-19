/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalGradEdgeAlgorithm_h
#define AssembleNodalGradEdgeAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;
class AssembleNodalGradEdgeAlgorithm : public Algorithm
{
public:

  AssembleNodalGradEdgeAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *scalarQ,
    VectorFieldType *dqdx);
  virtual ~AssembleNodalGradEdgeAlgorithm() {}

  virtual void execute();

  ScalarFieldType *scalarQ_;
  VectorFieldType *dqdx_;
  VectorFieldType *edgeAreaVec_;
  ScalarFieldType *dualNodalVolume_;

};

} // namespace nalu
} // namespace Sierra

#endif
