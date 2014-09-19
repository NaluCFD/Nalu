/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/

#ifndef AssembleNodalGradUEdgeAlgorithm_h
#define AssembleNodalGradUEdgeAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;
class AssembleNodalGradUEdgeAlgorithm : public Algorithm
{
public:

  AssembleNodalGradUEdgeAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    VectorFieldType *velocity,
    GenericFieldType *dudx);
  virtual ~AssembleNodalGradUEdgeAlgorithm() {}

  virtual void execute();

  VectorFieldType *velocity_;
  GenericFieldType *dudx_;
  VectorFieldType *edgeAreaVec_;
  ScalarFieldType *dualNodalVolume_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
