/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <EnthalpyPressureWorkNodeSuppAlg.h>
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
// EnthalpyPressureWorkNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
EnthalpyPressureWorkNodeSuppAlg::EnthalpyPressureWorkNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    dpdx_(NULL),
    velocity_(NULL),
    dualNodalVolume_(NULL),
    nDim_(realm_.meta_data().spatial_dimension())
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  dpdx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx");
  velocity_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
EnthalpyPressureWorkNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // pressure work
  const double *dpdx    = stk::mesh::field_data(*dpdx_, node );
  const double *velocity = stk::mesh::field_data(*velocity_, node );
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );

  double uDotGp = 0.0;
  for ( int j = 0; j < nDim_; ++j )
    uDotGp += velocity[j]*dpdx[j];

  rhs[0] += uDotGp*dualVolume;
  lhs[0] += 0.0;
}

} // namespace nalu
} // namespace Sierra
