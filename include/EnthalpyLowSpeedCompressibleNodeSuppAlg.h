/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef EnthalpyLowSpeedCompressibleNodeSuppAlg_h
#define EnthalpyLowSpeedCompressibleNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class EnthalpyLowSpeedCompressibleNodeSuppAlg : public SupplementalAlgorithm
{
public:

  EnthalpyLowSpeedCompressibleNodeSuppAlg(
    Realm &realm);

  virtual ~EnthalpyLowSpeedCompressibleNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  ScalarFieldType *pressureN_;
  ScalarFieldType *pressureNp1_;
  ScalarFieldType *dualNodalVolume_;
  double dt_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
