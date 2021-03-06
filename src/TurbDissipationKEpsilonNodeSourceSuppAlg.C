/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <TurbDissipationKEpsilonNodeSourceSuppAlg.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SupplementalAlgorithm.h>
#include <TimeIntegrator.h>

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

// c++
#include<limits>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// TurbDissipationKEpsilonNodeSourceSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TurbDissipationKEpsilonNodeSourceSuppAlg::TurbDissipationKEpsilonNodeSourceSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    cEpsOne_(realm.get_turb_model_constant(TM_cEpsOne)),
    cEpsTwo_(realm.get_turb_model_constant(TM_cEpsTwo)),
    epsNp1_(NULL),
    tkeNp1_(NULL),
    densityNp1_(NULL),
    tvisc_(NULL),
    dualNodalVolume_(NULL),
    tkeProdLimitRatio_(realm.get_turb_model_constant(TM_tkeProdLimitRatio)),
    includeDivU_(realm_.get_divU()),
    twoThirds_(2.0/3.0),
    nDim_(realm_.meta_data().spatial_dimension())
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  ScalarFieldType *eps = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_dissipation");
  epsNp1_ = &(eps->field_of_state(stk::mesh::StateNP1));
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
TurbDissipationKEpsilonNodeSourceSuppAlg::setup()
{
  // could extract user-based values
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
TurbDissipationKEpsilonNodeSourceSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  const double eps        = *stk::mesh::field_data(*epsNp1_, node );
  const double tke        = *stk::mesh::field_data(*tkeNp1_, node );
  const double rho        = *stk::mesh::field_data(*densityNp1_, node );
  const double tvisc      = *stk::mesh::field_data(*tvisc_, node );
  const double *dudx      =  stk::mesh::field_data(*dudx_, node );
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );

  int nDim = nDim_;
  double Pk = 0.0;
  double divU = 0.0;
  for ( int i = 0; i < nDim; ++i ) {
    const int offSet = nDim*i;
    divU += dudx[offSet+i]*includeDivU_;
    for ( int j = 0; j < nDim; ++j ) {
      Pk += dudx[offSet+j]*(dudx[offSet+j] + dudx[nDim*j+i]);
    }
  }
  Pk = tvisc*(Pk-twoThirds_*divU*divU) - twoThirds_*rho*tke*divU;
  
  // clip
  Pk = std::max(0.0, Pk);
  const double tkeC = std::max(std::numeric_limits<double>::min(), tke);

  const double Dk = rho*eps;
  if ( Pk > tkeProdLimitRatio_*Dk )
    Pk = tkeProdLimitRatio_*Dk;

  rhs[0] += eps/tkeC*(cEpsOne_*Pk - cEpsTwo_*Dk)*dualVolume;
  lhs[0] += cEpsTwo_*rho*eps/tkeC*dualVolume;
}

} // namespace nalu
} // namespace Sierra
