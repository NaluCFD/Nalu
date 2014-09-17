/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef EnthalpyPmrSrcNodeSuppAlg_h
#define EnthalpyPmrSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class EnthalpyPmrSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  EnthalpyPmrSrcNodeSuppAlg(
    Realm &realm);

  virtual ~EnthalpyPmrSrcNodeSuppAlg() {}

  virtual void setup() {}

  virtual void elem_execute(
    const int &numScvIntPoints,
    const int &numScsIntPoints,
    double *lhs,
    double *rhs,
    stk::mesh::Entity elem) {}

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  ScalarFieldType *divRadFlux_;
  ScalarFieldType *dualNodalVolume_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
