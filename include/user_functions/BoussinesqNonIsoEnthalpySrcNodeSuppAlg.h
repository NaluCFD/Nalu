/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef BoussinesqNonIsoEnthalpySrcNodeSuppAlg_h
#define BoussinesqNonIsoEnthalpySrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class BoussinesqNonIsoEnthalpySrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  BoussinesqNonIsoEnthalpySrcNodeSuppAlg(
    Realm &realm);

  virtual ~BoussinesqNonIsoEnthalpySrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *coordinates_;
  ScalarFieldType *dualNodalVolume_;
};

} // namespace nalu
} // namespace Sierra

#endif
