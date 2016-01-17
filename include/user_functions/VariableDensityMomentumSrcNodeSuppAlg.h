/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef VariableDensityMomentumSrcNodeSuppAlg_h
#define VariableDensityMomentumSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class VariableDensityMomentumSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  VariableDensityMomentumSrcNodeSuppAlg(
    Realm &realm);

  virtual ~VariableDensityMomentumSrcNodeSuppAlg() {}

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
  const double znot_;
  const double a_;
  const double amf_;
  const double visc_;
  const double rhoP_;
  const double rhoS_;
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
