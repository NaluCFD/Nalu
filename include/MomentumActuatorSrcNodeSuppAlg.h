/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumActuatorSrcNodeSuppAlg_h
#define MomentumActuatorSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class MomentumActuatorSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  MomentumActuatorSrcNodeSuppAlg(
    Realm &realm);

  virtual ~MomentumActuatorSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *actuatorSrc_;
  ScalarFieldType *actuatorSrcLHS_;
  ScalarFieldType *dualNodalVolume_;
  int nDim_;
};

} // namespace nalu
} // namespace Sierra

#endif
