/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumBoussinesqRASrcNodeSuppAlg_h
#define MomentumBoussinesqRASrcNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class MomentumBoussinesqRASrcNodeSuppAlg : public SupplementalAlgorithm
{
public:

  MomentumBoussinesqRASrcNodeSuppAlg(
    Realm &realm);

  virtual ~MomentumBoussinesqRASrcNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);

  double compute_alpha(double delta_t);
  double update_average(double avg, double newVal);

  ScalarFieldType *temperature_;
  ScalarFieldType *raTemperature_;
  std::string raName_;
  ScalarFieldType *dualNodalVolume_;
  double rhoRef_;
  double beta_;
  int nDim_;

  std::vector<double> gravity_;
};

} // namespace nalu
} // namespace Sierra

#endif
