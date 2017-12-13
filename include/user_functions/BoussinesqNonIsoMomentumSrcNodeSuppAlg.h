/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef BoussinesqNonIsoMomentumSrcNodeSuppAlg_h
#define BoussinesqNonIsoMomentumSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class BoussinesqNonIsoMomentumSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  BoussinesqNonIsoMomentumSrcNodeSuppAlg(
    Realm &realm);

  virtual ~BoussinesqNonIsoMomentumSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *coordinates_;
  ScalarFieldType *dualNodalVolume_;

  const double visc_;
  const double Cp;

  double beta;
  double rhoRef;
  double TRef;
  std::vector<double> gravity;
};

} // namespace nalu
} // namespace Sierra

#endif
