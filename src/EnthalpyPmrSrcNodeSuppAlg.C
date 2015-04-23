/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <EnthalpyPmrSrcNodeSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <TimeIntegrator.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// EnthalpyPmrSrcNodeSuppAlg - base class for algorithm
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyPmrSrcNodeSuppAlg::EnthalpyPmrSrcNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    divRadFlux_(NULL),
    dualNodalVolume_(NULL)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  divRadFlux_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "div_radiative_heat_flux");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
EnthalpyPmrSrcNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // explicit coupling for now...
  const double divQ       = *stk::mesh::field_data(*divRadFlux_, node );
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );

  rhs[0] -= divQ*dualVolume;
  lhs[0] += 0.0;
}

} // namespace nalu
} // namespace Sierra
