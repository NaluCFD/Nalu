/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SteadyThermalContactSrcNodeSuppAlg_h
#define SteadyThermalContactSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class SteadyThermalContactSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  SteadyThermalContactSrcNodeSuppAlg(
    Realm &realm);

  virtual ~SteadyThermalContactSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *coordinates_;
  ScalarFieldType *dualNodalVolume_;
  double a_;
  double k_;
  double pi_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
