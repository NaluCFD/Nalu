/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef VariableDensityContinuitySrcNodeSuppAlg_h
#define VariableDensityContinuitySrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class VariableDensityContinuitySrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  VariableDensityContinuitySrcNodeSuppAlg(
    Realm &realm);

  virtual ~VariableDensityContinuitySrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *coordinates_;
  ScalarFieldType *dualNodalVolume_;
  const double unot_;
  const double vnot_;
  const double wnot_;
  const double znot_;
  const double rhoP_;
  const double rhoS_;
  const double a_;
  const double amf_; 
  const double pi_;
  double projTimeScale_;
};

} // namespace nalu
} // namespace Sierra

#endif
