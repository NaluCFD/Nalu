/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumMassBackwardEulerNodeSuppAlg_h
#define MomentumMassBackwardEulerNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class MomentumMassBackwardEulerNodeSuppAlg : public SupplementalAlgorithm
{
public:

  MomentumMassBackwardEulerNodeSuppAlg(
    Realm &realm);

  virtual ~MomentumMassBackwardEulerNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *velocityN_;
  VectorFieldType *velocityNp1_;
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  VectorFieldType *dpdx_;
  ScalarFieldType *dualNodalVolume_;

  double dt_;
  int nDim_;
};

} // namespace nalu
} // namespace Sierra

#endif
