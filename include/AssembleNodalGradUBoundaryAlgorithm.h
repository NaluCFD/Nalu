/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/

#ifndef AssembleNodalGradUBoundaryAlgorithm_h
#define AssembleNodalGradUBoundaryAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalGradUBoundaryAlgorithm : public Algorithm
{
public:
  AssembleNodalGradUBoundaryAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    VectorFieldType *vectorQ,
    GenericFieldType *dqdx,
    const bool useShifted);
  virtual ~AssembleNodalGradUBoundaryAlgorithm() {}

  virtual void execute();

  VectorFieldType *vectorQ_;
  GenericFieldType *dqdx_;

  const bool useShifted_;

};

} // namespace nalu
} // namespace Sierra

#endif
