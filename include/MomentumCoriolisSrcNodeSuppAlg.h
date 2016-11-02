/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumCoriolisSrcNodeSuppAlg_h
#define MomentumCoriolisSrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class MomentumCoriolisSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  MomentumCoriolisSrcNodeSuppAlg(
    Realm &realm);

  virtual ~MomentumCoriolisSrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  void cross_product(std::vector<double> u, std::vector<double> v, std::vector<double> cross);

  ScalarFieldType *densityNp1_;
  ScalarFieldType *dualNodalVolume_;
  VectorFieldType *velocityNp1_;
  int nDim_;
  double earthAngularVelocity_;
  double latitude_;
  double sinphi_;
  double cosphi_;
  double corfac_;
  double Jxy_, Jxz_, Jyz_;
  double pi_;
  std::vector<double> eastVector_;
  std::vector<double> northVector_;
  std::vector<double> upVector_;

};

} // namespace nalu
} // namespace Sierra

#endif
