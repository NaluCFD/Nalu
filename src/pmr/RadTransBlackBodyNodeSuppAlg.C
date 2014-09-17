/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <pmr/RadTransBlackBodyNodeSuppAlg.h>
#include <pmr/RadiativeTransportEquationSystem.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <TimeIntegrator.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
// stk_io
#include <stk_io/StkMeshIoBroker.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// RadTransBlackBodyNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
RadTransBlackBodyNodeSuppAlg::RadTransBlackBodyNodeSuppAlg(
  Realm &realm,
  RadiativeTransportEquationSystem *radEqSystem)
  : SupplementalAlgorithm(realm),
    radEqSystem_(radEqSystem),
    intensity_(NULL),
    absorption_(NULL),
    scattering_(NULL),
    temperature_(NULL),
    dualNodalVolume_(NULL),
    invPi_(1.0/std::acos(-1.0)),
    sb_(radEqSystem->get_stefan_boltzmann())
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.fixture_->meta_data();
  absorption_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "absorption_coefficient");
  scattering_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "scattering_coefficient");
  temperature_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
RadTransBlackBodyNodeSuppAlg::setup()
{
  intensity_ = radEqSystem_->get_intensity();
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
RadTransBlackBodyNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // ((mu+beta)*I - mu*sigma*T^4/pi)*dVol = 0.0
  const double intensity  = *stk::mesh::field_data(*intensity_, node);
  const double mu = *stk::mesh::field_data(*absorption_, node);
  const double beta = *stk::mesh::field_data(*scattering_, node);
  const double temp       = *stk::mesh::field_data(*temperature_, node);
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node);
  rhs[0] -= ((mu+beta)*intensity - mu*sb_*temp*temp*temp*temp*invPi_)*dualVolume;
  lhs[0] += (mu+beta)*dualVolume;
}

} // namespace nalu
} // namespace Sierra
