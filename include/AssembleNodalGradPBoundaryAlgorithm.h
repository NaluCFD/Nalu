/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalGradPBoundaryAlgorithm_h
#define AssembleNodalGradPBoundaryAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalGradPBoundaryAlgorithm : public Algorithm
{
public:
  AssembleNodalGradPBoundaryAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *pressure,
    VectorFieldType *dpdx,
    const bool useShifted);
  virtual ~AssembleNodalGradPBoundaryAlgorithm() {}

  virtual void execute();

  ScalarFieldType *pressure_;
  VectorFieldType *dpdx_;
  const bool useShifted_;
};

} // namespace nalu
} // namespace Sierra

#endif
