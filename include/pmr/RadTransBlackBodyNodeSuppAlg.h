/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef RadTransBlackBodyNodeSuppAlg_h
#define RadTransBlackBodyNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class RadiativeTransportEquationSystem;
class Realm;

class RadTransBlackBodyNodeSuppAlg : public SupplementalAlgorithm
{
public:

  RadTransBlackBodyNodeSuppAlg(
      Realm &realm,
      RadiativeTransportEquationSystem *radEqSystem);

  virtual ~RadTransBlackBodyNodeSuppAlg() {}

  virtual void setup();

  virtual void elem_execute(
      const int &numScvIntPoints,
      const int &numScsIntPoints,
      double *lhs,
      double *rhs,
      stk::mesh::Entity elem) {}

  virtual void node_execute(
      double *lhs,
      double *rhs,
      stk::mesh::Entity node);
 
  const RadiativeTransportEquationSystem *radEqSystem_;
  ScalarFieldType *intensity_;
  ScalarFieldType *absorption_;
  ScalarFieldType *scattering_;
  ScalarFieldType *temperature_;
  ScalarFieldType *dualNodalVolume_;

  const double invPi_;
  const double sb_;

};

} // namespace nalu
} // namespace Sierra

#endif
