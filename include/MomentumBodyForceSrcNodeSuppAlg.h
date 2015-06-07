/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumBodyForceSrcNodeSuppAlg_h
#define MomentumBodyForceSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class MomentumBodyForceSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  MomentumBodyForceSrcNodeSuppAlg(
    Realm &realm,
    const std::vector<double> theParams);

  virtual ~MomentumBodyForceSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  std::vector<double> params_;
  ScalarFieldType *dualNodalVolume_;
  int nDim_;

};

} // namespace nalu
} // namespace Sierra

#endif
