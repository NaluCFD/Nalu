/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MeshDisplacementMassBackwardEulerNodeSuppAlg_h
#define MeshDisplacementMassBackwardEulerNodeSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;

class MeshDisplacementMassBackwardEulerNodeSuppAlg : public SupplementalAlgorithm
{
public:

  MeshDisplacementMassBackwardEulerNodeSuppAlg(
    Realm &realm);

  virtual ~MeshDisplacementMassBackwardEulerNodeSuppAlg() {}

  virtual void setup();

  virtual void node_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity node);
  
  VectorFieldType *displacementN_;
  VectorFieldType *displacementNp1_;
  ScalarFieldType *density_;
  ScalarFieldType *dualNodalVolume_;

  double dt_;
  int nDim_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
