/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SteadyTaylorVortexMixFracSrcNodeSuppAlg_h
#define SteadyTaylorVortexMixFracSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class SteadyTaylorVortexMixFracSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  SteadyTaylorVortexMixFracSrcNodeSuppAlg(
    Realm &realm);

  virtual ~SteadyTaylorVortexMixFracSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *coordinates_;
  ScalarFieldType *dualNodalVolume_;
  const double rhoP_;
  const double rhoS_;
  const double unot_;
  const double vnot_;
  const double znot_;
  const double pnot_;
  const double visc_;
  const double a_;
  const double amf_;
  const double Sc_;  
  const double pi_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
