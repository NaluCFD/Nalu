/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <TurbKineticEnergyKsgsNodeSourceSuppAlg.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SupplementalAlgorithm.h>
#include <TimeIntegrator.h>
#include <stk_mesh/base/Field.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// TurbKineticEnergyKsgsNodeSourceSuppAlg - Ksgs LES source term algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TurbKineticEnergyKsgsNodeSourceSuppAlg::TurbKineticEnergyKsgsNodeSourceSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    tkeNp1_(NULL),
    densityNp1_(NULL),
    tvisc_(NULL),
    dudx_(NULL),
    dualNodalVolume_(NULL),
    cEps_(realm_.get_turb_model_constant(TM_cEps)),
    tkeProdLimitRatio_(realm_.get_turb_model_constant(TM_tkeProdLimitRatio)),
    nDim_(realm_.meta_data().spatial_dimension())
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  ScalarFieldType *tke = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_ke");
  tkeNp1_ = &(tke->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  tvisc_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity");
  dudx_ = meta_data.get_field<GenericFieldType>(stk::topology::NODE_RANK, "dudx");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
TurbKineticEnergyKsgsNodeSourceSuppAlg::setup()
{
  // could extract user-based values for cEps_ and tkeProdLimitRatio_
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
TurbKineticEnergyKsgsNodeSourceSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  const double tke        = *stk::mesh::field_data(*tkeNp1_, node );
  const double rho        = *stk::mesh::field_data(*densityNp1_, node );
  const double tvisc      = *stk::mesh::field_data(*tvisc_, node );
  const double *dudx      =  stk::mesh::field_data(*dudx_, node );
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );

  // filter
  double filter = std::pow(dualVolume, 1.0/nDim_);

  int nDim = nDim_;
  double Pk = 0.0;
  for ( int i = 0; i < nDim; ++i ) {
    const int offSet = nDim*i;
    for ( int j = 0; j < nDim; ++j ) {
      Pk += dudx[offSet+j]*(dudx[offSet+j] + dudx[nDim*j+i]);
    }
  }
  Pk *= tvisc;

  double Dk = cEps_*rho*std::pow(tke, 1.5)/filter;

  if ( Pk > tkeProdLimitRatio_*Dk )
    Pk = tkeProdLimitRatio_*Dk;

  rhs[0] += (Pk - Dk)*dualVolume;
  lhs[0] += 1.5*cEps_*rho*std::sqrt(tke)/filter*dualVolume;
}

} // namespace nalu
} // namespace Sierra
