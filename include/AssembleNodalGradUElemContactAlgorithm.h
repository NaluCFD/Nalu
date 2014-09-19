/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalGradUElemContactAlgorithm_h
#define AssembleNodalGradUElemContactAlgorithm_h

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

class AssembleNodalGradUElemContactAlgorithm : public Algorithm
{
public:
  AssembleNodalGradUElemContactAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    VectorFieldType *vectorQ,
    GenericFieldType *dqdx,
    VectorFieldType *haloQ);
  virtual ~AssembleNodalGradUElemContactAlgorithm() {}

  virtual void execute();

  VectorFieldType *vectorQ_;
  GenericFieldType *dqdx_;
  VectorFieldType *haloQ_;
  ScalarFieldType *dualNodalVolume_;

  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;

  void populate_halo_state();
  void add_elem_gradq();

};

} // namespace nalu
} // namespace Sierra

#endif
