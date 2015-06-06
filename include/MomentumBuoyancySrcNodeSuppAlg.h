/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumBuoyancySrcNodeSuppAlg_h
#define MomentumBuoyancySrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class MomentumBuoyancySrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  MomentumBuoyancySrcNodeSuppAlg(
    Realm &realm);

  virtual ~MomentumBuoyancySrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  ScalarFieldType *densityNp1_;
  ScalarFieldType *dualNodalVolume_;
  int nDim_;
  double rhoRef_;
  std::vector<double> gravity_;

};

} // namespace nalu
} // namespace Sierra

#endif
