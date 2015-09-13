/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SteadyTaylorVortexMomentumSrcNodeSuppAlg_h
#define SteadyTaylorVortexMomentumSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class SteadyTaylorVortexMomentumSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  SteadyTaylorVortexMomentumSrcNodeSuppAlg(
    Realm &realm);

  virtual ~SteadyTaylorVortexMomentumSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *coordinates_;
  ScalarFieldType *dualNodalVolume_;

  const int nDim_;
  const double unot_;
  const double a_;
  const double visc_;
  const double pi_;

  std::vector<double> srcXi_;
};

} // namespace nalu
} // namespace Sierra

#endif
