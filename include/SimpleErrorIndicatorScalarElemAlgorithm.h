/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SimpleErrorIndicatorScalarElemAlgorithm_h
#define SimpleErrorIndicatorScalarElemAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class SimpleErrorIndicatorScalarElemAlgorithm : public Algorithm
{
public:

  SimpleErrorIndicatorScalarElemAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *scalarQ,
    VectorFieldType *dqdx
  );
  ~SimpleErrorIndicatorScalarElemAlgorithm();

  void execute();

  // extract fields; nodal
  ScalarFieldType *velocity_;
  VectorFieldType *coordinates_;
  VectorFieldType *dudx_;
  GenericFieldType *errorIndicatorField_;
  int nUnk_;
};

} // namespace nalu
} // namespace Sierra

#endif
