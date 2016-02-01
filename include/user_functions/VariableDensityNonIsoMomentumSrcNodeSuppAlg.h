/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef VariableDensityNonIsoMomentumSrcNodeSuppAlg_h
#define VariableDensityNonIsoMomentumSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class VariableDensityNonIsoMomentumSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  VariableDensityNonIsoMomentumSrcNodeSuppAlg(
    Realm &realm);

  virtual ~VariableDensityNonIsoMomentumSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *coordinates_;
  ScalarFieldType *dualNodalVolume_;

  const int nDim_;
  const double unot_;
  const double vnot_;
  const double wnot_;
  const double pnot_;
  const double hnot_;
  const double a_;
  const double ah_;
  const double visc_;
  const double Pref_;
  const double MW_;
  const double R_;
  const double Tref_;
  const double Cp_;
  const double pi_;
  const double twoThirds_;
  double rhoRef_;
  double gx_;
  double gy_;
  double gz_;
  
  // space for source terms
  std::vector<double> srcXi_;
};

} // namespace nalu
} // namespace Sierra

#endif
