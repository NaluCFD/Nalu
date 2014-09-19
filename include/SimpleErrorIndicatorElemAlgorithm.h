/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SimpleErrorIndicatorElemAlgorithm_h
#define SimpleErrorIndicatorElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class SimpleErrorIndicatorElemAlgorithm : public Algorithm
{
public:

  SimpleErrorIndicatorElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part);
  ~SimpleErrorIndicatorElemAlgorithm();

  void execute();

  // extract fields; nodal
  VectorFieldType *velocity_;
  VectorFieldType *coordinates_;
  GenericFieldType *dudx_;
  GenericFieldType *errorIndicatorField_;

};

} // namespace nalu
} // namespace Sierra

#endif
