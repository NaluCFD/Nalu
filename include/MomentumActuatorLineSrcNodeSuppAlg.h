/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumActuatorLineSrcNodeSuppAlg_h
#define MomentumActuatorLineSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class MomentumActuatorLineSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  MomentumActuatorLineSrcNodeSuppAlg(
    Realm &realm);

  virtual ~MomentumActuatorLineSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *actuatorLineSrc_;
  ScalarFieldType *actuatorLineSrcLHS_;
  ScalarFieldType *dualNodalVolume_;
  int nDim_;
};

} // namespace nalu
} // namespace Sierra

#endif
