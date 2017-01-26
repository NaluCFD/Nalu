/******************************************************************************/
/* This software is released under the license detailed in the file, LICENSE, */
/* located in the top-level Nalu directory structure.                         */
/******************************************************************************/

#ifndef ENTHALPYABLSRCNODESUPPALG_H
#define ENTHALPYABLSRCNODESUPPALG_H

#include "SupplementalAlgorithm.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/Entity.hpp>

namespace sierra {
namespace nalu {

// Forward declarations
class Realm;
class ABLForcingAlgorithm;

class EnthalpyABLSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:
  EnthalpyABLSrcNodeSuppAlg(Realm&, ABLForcingAlgorithm*);

  virtual ~EnthalpyABLSrcNodeSuppAlg() {}

  virtual void setup() {}

  virtual void node_execute(double*, double*, stk::mesh::Entity);

private:
  EnthalpyABLSrcNodeSuppAlg();
  EnthalpyABLSrcNodeSuppAlg(const EnthalpyABLSrcNodeSuppAlg&);

  //! Pointer to ABL Forcing algorithm object
  ABLForcingAlgorithm* ablSrc_;

  //! Pointer to mesh coordinates
  VectorFieldType* coords_;

  //! Pointer to dual volume of the mesh
  ScalarFieldType* dualNodalVolume_;

  //! Pointer to the density of the mesh
  ScalarFieldType* density_;

  //! Pointer to specific heat
  ScalarFieldType* specificHeat_;

  const int nDim_;
};

} // namespace nalu
} // namespace sierra

#endif /* ENTHALPYABLSRCNODESUPPALG_H */
