/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/

#ifndef AssembleNodalGradUElemBoundaryAlgorithm_h
#define AssembleNodalGradUElemBoundaryAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalGradUElemBoundaryAlgorithm : public Algorithm
{
public:
  AssembleNodalGradUElemBoundaryAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    VectorFieldType *vectorQ,
    GenericFieldType *dqdx,
    const bool useShifted = false);
  virtual ~AssembleNodalGradUElemBoundaryAlgorithm() {}

  virtual void execute();

  VectorFieldType *vectorQ_;
  GenericFieldType *dqdx_;

  const bool useShifted_;

};

} // namespace nalu
} // namespace Sierra

#endif
