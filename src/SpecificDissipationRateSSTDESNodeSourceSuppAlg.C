/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <SpecificDissipationRateSSTDESNodeSourceSuppAlg.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SupplementalAlgorithm.h>
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
// SpecificDissipationRateSSTDESNodeSourceSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SpecificDissipationRateSSTDESNodeSourceSuppAlg::SpecificDissipationRateSSTDESNodeSourceSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    sigmaWTwo_(realm.get_turb_model_constant(TM_sigmaWTwo)),
    betaStar_(realm.get_turb_model_constant(TM_betaStar)),
    betaOne_(realm.get_turb_model_constant(TM_betaOne)),
    betaTwo_(realm.get_turb_model_constant(TM_betaTwo)),
    gammaOne_(realm.get_turb_model_constant(TM_gammaOne)),
    gammaTwo_(realm.get_turb_model_constant(TM_gammaTwo)),
    sdrNp1_(NULL),
    tkeNp1_(NULL),
    densityNp1_(NULL),
    fOneBlend_(NULL),
    tvisc_(NULL),
    dkdx_(NULL),
    dwdx_(NULL),
    dualNodalVolume_(NULL),
    maxLengthScale_(NULL),
    tkeProdLimitRatio_(realm.get_turb_model_constant(TM_tkeProdLimitRatio)),
    nDim_(realm_.meta_data().spatial_dimension()),
    cDESke_(realm.get_turb_model_constant(TM_cDESke)),
    cDESkw_(realm.get_turb_model_constant(TM_cDESkw))
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  ScalarFieldType *sdr = meta_data.get_field<double>(stk::topology::NODE_RANK, "specific_dissipation_rate");
  sdrNp1_ = &(sdr->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *tke = meta_data.get_field<double>(stk::topology::NODE_RANK, "turbulent_ke");
  tkeNp1_ = &(tke->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *density = meta_data.get_field<double>(stk::topology::NODE_RANK, "density");
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  fOneBlend_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "sst_f_one_blending");
  tvisc_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "turbulent_viscosity");
  dudx_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "dudx");
  dkdx_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "dkdx");
  dwdx_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "dwdx");
  maxLengthScale_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "sst_max_length_scale");
  dualNodalVolume_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "dual_nodal_volume");
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
SpecificDissipationRateSSTDESNodeSourceSuppAlg::setup()
{
  // could extract user-based values
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
SpecificDissipationRateSSTDESNodeSourceSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  const double sdr        = *stk::mesh::field_data(*sdrNp1_, node );
  const double tke        = *stk::mesh::field_data(*tkeNp1_, node );
  const double rho        = *stk::mesh::field_data(*densityNp1_, node );
  const double fOneBlend  = *stk::mesh::field_data(*fOneBlend_, node );
  const double tvisc      = *stk::mesh::field_data(*tvisc_, node );
  const double *dudx      =  stk::mesh::field_data(*dudx_, node );
  const double *dkdx      =  stk::mesh::field_data(*dkdx_, node );
  const double *dwdx      =  stk::mesh::field_data(*dwdx_, node );
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double maxLengthScale = *stk::mesh::field_data(*maxLengthScale_, node );

  int nDim = nDim_;
  double Pk = 0.0;
  double crossDiff = 0.0;
  for ( int i = 0; i < nDim; ++i ) {
    crossDiff += dkdx[i]*dwdx[i];
    const int offSet = nDim*i;
    for ( int j = 0; j < nDim; ++j ) {
      Pk += dudx[offSet+j]*(dudx[offSet+j] + dudx[nDim*j+i]);
    }
  }
  Pk *= tvisc;

  // blend cDES  
  const double cDES = fOneBlend*cDESkw_ + (1.0-fOneBlend)*cDESke_;

  const double sqrtTke = std::sqrt(tke);
  const double lSST = sqrtTke/betaStar_/sdr;
  // prevent divide by zero (possible at resolved tke bcs = 0.0)
  const double lDES = std::max(1.0e-16, std::min(lSST, cDES*maxLengthScale));
  
  double Dk = rho*tke*sqrtTke/lDES;

  if ( Pk > tkeProdLimitRatio_*Dk )
    Pk = tkeProdLimitRatio_*Dk;

  // start the blending and constants
  const double om_fOneBlend = 1.0 - fOneBlend;
  const double beta = fOneBlend*betaOne_ + om_fOneBlend*betaTwo_;
  const double gamma = fOneBlend*gammaOne_ + om_fOneBlend*gammaTwo_;
  const double sigmaD = 2.0*om_fOneBlend*sigmaWTwo_;

  // Pw includes 1/tvisc scaling; tvisc may be zero at a dirichlet low Re approach (clip)
  const double Pw = gamma*rho*Pk/std::max(tvisc, 1.0e-16);
  const double Dw = beta*rho*sdr*sdr;
  const double Sw = sigmaD*rho*crossDiff/sdr;

  rhs[0] += (Pw - Dw + Sw)*dualVolume;
  lhs[0] += (2.0*beta*rho*sdr + std::max(Sw/sdr,0.0))*dualVolume;
}

} // namespace nalu
} // namespace Sierra
