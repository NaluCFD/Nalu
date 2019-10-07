/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef TemperaturePmrSrcNodeSuppAlg_h
#define TemperaturePmrSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class TemperaturePmrSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  TemperaturePmrSrcNodeSuppAlg(
    Realm &realm);

  virtual ~TemperaturePmrSrcNodeSuppAlg() {}

  virtual void setup() {}

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  ScalarFieldType *divRadFlux_;
  ScalarFieldType *divRadFluxLin_;
  ScalarFieldType *dualNodalVolume_;
};

} // namespace nalu
} // namespace Sierra

#endif
