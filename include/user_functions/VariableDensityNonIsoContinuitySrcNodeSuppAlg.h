/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef VariableDensityNonIsoContinuitySrcNodeSuppAlg_h
#define VariableDensityNonIsoContinuitySrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class VariableDensityNonIsoContinuitySrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  VariableDensityNonIsoContinuitySrcNodeSuppAlg(
    Realm &realm);

  virtual ~VariableDensityNonIsoContinuitySrcNodeSuppAlg() {}

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
  const double hnot_;
  const double a_;
  const double ah_;
  const double Pref_;
  const double MW_;
  const double R_; 
  const double Tref_;
  const double Cp_;
  const double Pr_;
  const double pi_;
  double projTimeScale_;
};

} // namespace nalu
} // namespace Sierra

#endif
