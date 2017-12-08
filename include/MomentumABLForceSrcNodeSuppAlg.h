/******************************************************************************/
/* This software is released under the license detailed in the file, LICENSE, */
/* located in the top-level Nalu directory structure.                         */
/******************************************************************************/

#ifndef MOMENTUMABLFORCESRCNODESUPPALG_H
#define MOMENTUMABLFORCESRCNODESUPPALG_H

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra {
namespace nalu {

class Realm;
class ABLForcingAlgorithm;

class MomentumABLForceSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:
  MomentumABLForceSrcNodeSuppAlg(Realm&, ABLForcingAlgorithm*);

  virtual ~MomentumABLForceSrcNodeSuppAlg() {}

  virtual void setup() {}

  virtual void node_execute(double*, double*, stk::mesh::Entity);

private:
  MomentumABLForceSrcNodeSuppAlg();
  MomentumABLForceSrcNodeSuppAlg(const MomentumABLForceSrcNodeSuppAlg&);

  //! Pointer to ABL Forcing Algorithm object
  ABLForcingAlgorithm* ablSrc_;

  //! Pointer to the mesh coordinates
  VectorFieldType* coords_;

  //! Pointer to the dual volume of the mesh
  ScalarFieldType* dualNodalVolume_;

  //! Spatial dimension of the computational mesh
  const int nDim_;
};
}
}

#endif /* MOMENTUMABLFORCESRCNODESUPPALG_H */
