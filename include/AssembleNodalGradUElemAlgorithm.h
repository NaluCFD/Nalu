/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalGradUElemAlgorithm_h
#define AssembleNodalGradUElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;
class AssembleNodalGradUElemAlgorithm : public Algorithm
{
public:

  AssembleNodalGradUElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    VectorFieldType *velocity,
    GenericFieldType *dudx,
    const bool useShifted = false);
  virtual ~AssembleNodalGradUElemAlgorithm() {}

  virtual void execute();

  VectorFieldType *vectorQ_;
  GenericFieldType *dqdx_;
  const bool useShifted_;

};

} // namespace nalu
} // namespace Sierra

#endif
