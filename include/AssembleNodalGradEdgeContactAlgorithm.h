/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalGradEdgeContactAlgorithm_h
#define AssembleNodalGradEdgeContactAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class FieldBase;
}
}

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalGradEdgeContactAlgorithm : public Algorithm
{
public:
  AssembleNodalGradEdgeContactAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *scalarQ,
    VectorFieldType *dqdx);
  virtual ~AssembleNodalGradEdgeContactAlgorithm() {}

  virtual void execute();

  ScalarFieldType *scalarQ_;
  VectorFieldType *dqdx_;
  ScalarFieldType *dualNodalVolume_;

  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;

};

} // namespace nalu
} // namespace Sierra

#endif
