/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <HeatCondMassBackwardEulerNodeSuppAlg.h>
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
// HeatCondMassBackwardEulerNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
HeatCondMassBackwardEulerNodeSuppAlg::HeatCondMassBackwardEulerNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    temperatureN_(NULL),
    temperatureNp1_(NULL),
    density_(NULL),
    specificHeat_(NULL),
    dualNodalVolume_(NULL),
    dt_(0.0)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  ScalarFieldType *temperature = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");
  temperatureN_ = &(temperature->field_of_state(stk::mesh::StateN));
  temperatureNp1_ = &(temperature->field_of_state(stk::mesh::StateNP1));
  density_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  specificHeat_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_heat");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
HeatCondMassBackwardEulerNodeSuppAlg::setup()
{
  dt_ = realm_.timeIntegrator_->get_time_step();
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
HeatCondMassBackwardEulerNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // deal with lumped mass matrix
  const double tN         = *stk::mesh::field_data(*temperatureN_, node);
  const double tNp1       = *stk::mesh::field_data(*temperatureNp1_, node);
  const double rho     = *stk::mesh::field_data(*density_, node);
  const double cp     = *stk::mesh::field_data(*specificHeat_, node);
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node);
  const double lhsTime = rho*cp*dualVolume/dt_;
  rhs[0] -= lhsTime*(tNp1 - tN);
  lhs[0] += lhsTime;
}

} // namespace nalu
} // namespace Sierra
