/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef RadTransIsoScatteringNodeSuppAlg_h
#define RadTransIsoScatteringNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class RadiativeTransportEquationSystem;
class Realm;

//--------------------------------------------------------------------
class RadTransIsoScatteringNodeSuppAlg : public SupplementalAlgorithm
{
public:

  RadTransIsoScatteringNodeSuppAlg(
    Realm &realm,
    RadiativeTransportEquationSystem *radEqSystem);

  virtual ~RadTransIsoScatteringNodeSuppAlg() {}

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
 
  const RadiativeTransportEquationSystem *radEqSystem_;
  ScalarFieldType *scalarFlux_;
  ScalarFieldType *scattering_;
  ScalarFieldType *dualNodalVolume_;

  const double invPi_;

};

} // namespace nalu
} // namespace Sierra

#endif
