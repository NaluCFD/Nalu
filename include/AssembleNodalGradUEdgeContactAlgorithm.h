/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalGradUEdgeContactAlgorithm_h
#define AssembleNodalGradUEdgeContactAlgorithm_h

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

class AssembleNodalGradUEdgeContactAlgorithm : public Algorithm
{
public:
  AssembleNodalGradUEdgeContactAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    VectorFieldType *velocity,
    GenericFieldType *dudx);
  virtual ~AssembleNodalGradUEdgeContactAlgorithm() {}

  virtual void execute();

  VectorFieldType *velocity_;
  GenericFieldType *dudx_;
  ScalarFieldType *dualNodalVolume_;

  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;

};

} // namespace nalu
} // namespace Sierra

#endif
