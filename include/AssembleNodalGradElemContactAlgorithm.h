/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalGradElemContactAlgorithm_h
#define AssembleNodalGradElemContactAlgorithm_h

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

class AssembleNodalGradElemContactAlgorithm : public Algorithm
{
public:
  AssembleNodalGradElemContactAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *scalarQ,
    VectorFieldType *dqdx,
    ScalarFieldType *haloQ);
  virtual ~AssembleNodalGradElemContactAlgorithm() {}

  virtual void execute();

  ScalarFieldType *scalarQ_;
  VectorFieldType *dqdx_;
  ScalarFieldType *haloQ_;
  ScalarFieldType *dualNodalVolume_;

  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;

  void populate_halo_state();
  void add_elem_gradq();

};

} // namespace nalu
} // namespace Sierra

#endif
