/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ContinuityMassBackwardEulerNodeSuppAlg_h
#define ContinuityMassBackwardEulerNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ContinuityMassBackwardEulerNodeSuppAlg : public SupplementalAlgorithm
{
public:

  ContinuityMassBackwardEulerNodeSuppAlg(
    Realm &realm);

  virtual ~ContinuityMassBackwardEulerNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  ScalarFieldType *dualNodalVolume_;
  double dt_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
