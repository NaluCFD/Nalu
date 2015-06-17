/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalGradUNonConformalAlgorithm_h
#define AssembleNodalGradUNonConformalAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalGradUNonConformalAlgorithm : public Algorithm
{
public:

  AssembleNodalGradUNonConformalAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    VectorFieldType *vectorQ,
    GenericFieldType *dqdx);

  ~AssembleNodalGradUNonConformalAlgorithm();

  void execute();

  VectorFieldType *vectorQ_;
  GenericFieldType *dqdx_;
  
  ScalarFieldType *dualNodalVolume_;
  GenericFieldType *exposedAreaVec_;

  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
