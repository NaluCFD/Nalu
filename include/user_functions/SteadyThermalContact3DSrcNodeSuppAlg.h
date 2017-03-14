/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SteadyThermalContact3DSrcNodeSuppAlg_h
#define SteadyThermalContact3DSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class SteadyThermalContact3DSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  SteadyThermalContact3DSrcNodeSuppAlg(
    Realm &realm);

  virtual ~SteadyThermalContact3DSrcNodeSuppAlg() {}

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
