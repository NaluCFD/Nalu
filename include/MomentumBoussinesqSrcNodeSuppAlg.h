/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumBoussinesqSrcNodeSuppAlg_h
#define MomentumBoussinesqSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class MomentumBoussinesqSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  MomentumBoussinesqSrcNodeSuppAlg(
    Realm &realm);

  virtual ~MomentumBoussinesqSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  ScalarFieldType *temperature_;
  ScalarFieldType *dualNodalVolume_;
  double tRef_;
  double rhoRef_;
  double beta_;
  int nDim_;
  std::array<double, 3> gravity_ = {0.0, 0.0, 0.0};

};

} // namespace nalu
} // namespace Sierra

#endif
