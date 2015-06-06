/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef EnthalpyPressureWorkNodeSuppAlg_h
#define EnthalpyPressureWorkNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class EnthalpyPressureWorkNodeSuppAlg : public SupplementalAlgorithm
{
public:

  EnthalpyPressureWorkNodeSuppAlg(
    Realm &realm);

  virtual ~EnthalpyPressureWorkNodeSuppAlg() {}

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *dpdx_;
  VectorFieldType *velocity_;
  ScalarFieldType *dualNodalVolume_;
  const int nDim_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
