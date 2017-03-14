/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumCoriolisSrcNodeSuppAlg_h
#define MomentumCoriolisSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <CoriolisSrc.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class MomentumCoriolisSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  MomentumCoriolisSrcNodeSuppAlg(
    Realm &realm);

  virtual ~MomentumCoriolisSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  ScalarFieldType *densityNp1_;
  VectorFieldType *velocityNp1_;
  ScalarFieldType *dualNodalVolume_;
  
  CoriolisSrc cor_;
};

} // namespace nalu
} // namespace Sierra

#endif
