/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalGradPAWElemAlgorithm_h
#define AssembleNodalGradPAWElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;
class AssembleNodalGradPAWElemAlgorithm : public Algorithm
{
public:

  AssembleNodalGradPAWElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *pressure,
    VectorFieldType *dpdx);
  virtual ~AssembleNodalGradPAWElemAlgorithm() {}

  virtual void execute();

  ScalarFieldType *pressure_;
  ScalarFieldType *density_;
  VectorFieldType *dpdx_;
  VectorFieldType *areaWeight_;
  const bool useShifted_;
};

} // namespace nalu
} // namespace Sierra

#endif
