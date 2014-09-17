/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalGradElemBoundaryAlgorithm_h
#define AssembleNodalGradElemBoundaryAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalGradElemBoundaryAlgorithm : public Algorithm
{
public:
  AssembleNodalGradElemBoundaryAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *scalarQ,
    VectorFieldType *dqdx,
    const bool useShifted = false);
  virtual ~AssembleNodalGradElemBoundaryAlgorithm() {}

  virtual void execute();

  ScalarFieldType *scalarQ_;
  VectorFieldType *dqdx_;
  const bool useShifted_;

};

} // namespace nalu
} // namespace Sierra

#endif
