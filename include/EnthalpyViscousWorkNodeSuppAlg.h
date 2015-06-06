/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef EnthalpyViscousWorkNodeSuppAlg_h
#define EnthalpyViscousWorkNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class EnthalpyViscousWorkNodeSuppAlg : public SupplementalAlgorithm
{
public:

  EnthalpyViscousWorkNodeSuppAlg(
    Realm &realm);

  virtual ~EnthalpyViscousWorkNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  GenericFieldType *dudx_;
  ScalarFieldType *viscosity_;
  ScalarFieldType *dualNodalVolume_;
  const double includeDivU_;
  const int nDim_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
