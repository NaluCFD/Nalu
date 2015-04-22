/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ContinuityMassBackwardEulerNodeSuppAlg.h>
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
// ContinuityMassBackwardEulerNodeSuppAlg - lumped mass drho/dt; BE
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ContinuityMassBackwardEulerNodeSuppAlg::ContinuityMassBackwardEulerNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    densityN_(NULL),
    densityNp1_(NULL),
    dualNodalVolume_(NULL),
    dt_(0.0)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityN_ = &(density->field_of_state(stk::mesh::StateN));
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
ContinuityMassBackwardEulerNodeSuppAlg::setup()
{
  dt_ = realm_.timeIntegrator_->get_time_step();
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
ContinuityMassBackwardEulerNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // deal with lumped mass matrix
  const double projTimeScale = dt_;
  const double rhoN       = *stk::mesh::field_data(*densityN_, node );
  const double rhoNp1     = *stk::mesh::field_data(*densityNp1_, node );
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  rhs[0] -= (rhoNp1 - rhoN)*dualVolume/dt_/projTimeScale;
  lhs[0] += 0.0;
}

} // namespace nalu
} // namespace Sierra
