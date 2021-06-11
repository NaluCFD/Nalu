/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalGradPAWBoundaryAlgorithm_h
#define AssembleNodalGradPAWBoundaryAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalGradPAWBoundaryAlgorithm : public Algorithm
{
public:
  AssembleNodalGradPAWBoundaryAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *pressure,
    VectorFieldType *dpdx,
    const std::string bcPressureName = "pressure");
  virtual ~AssembleNodalGradPAWBoundaryAlgorithm() {}

  virtual void execute();

  ScalarFieldType *pressure_;
  VectorFieldType *dpdx_;
  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  ScalarFieldType *bcPressure_;
  VectorFieldType *areaWeight_;
  const bool useShifted_;
};

} // namespace nalu
} // namespace Sierra

#endif
