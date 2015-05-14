/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ContinuityMassBDF2NodeSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>

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
// ContinuityMassBDF2NodeSuppAlg - lumped mass drho/dt; BDF2
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ContinuityMassBDF2NodeSuppAlg::ContinuityMassBDF2NodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    densityNm1_(NULL),
    densityN_(NULL),
    densityNp1_(NULL),
    dualNodalVolume_(NULL),
    dt_(0.0),
    gamma1_(0.0),
    gamma2_(0.0),
    gamma3_(0.0)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNm1_ = &(density->field_of_state(stk::mesh::StateNM1));
  densityN_ = &(density->field_of_state(stk::mesh::StateN));
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
ContinuityMassBDF2NodeSuppAlg::setup()
{
  dt_ = realm_.get_time_step();
  gamma1_ = realm_.get_gamma1();
  gamma2_ = realm_.get_gamma2();
  gamma3_ = realm_.get_gamma3();
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
ContinuityMassBDF2NodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // deal with lumped mass matrix
  const double projTimeScale = dt_/gamma1_;
  const double rhoNm1     = *stk::mesh::field_data(*densityNm1_, node);
  const double rhoN       = *stk::mesh::field_data(*densityN_, node);
  const double rhoNp1     = *stk::mesh::field_data(*densityNp1_, node);
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node);

  rhs[0] -= (gamma1_*rhoNp1 + gamma2_*rhoN + gamma3_*rhoNm1)*dualVolume/dt_/projTimeScale;
  lhs[0] += 0.0;
}

} // namespace nalu
} // namespace Sierra
