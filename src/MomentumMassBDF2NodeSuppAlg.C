/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <MomentumMassBDF2NodeSuppAlg.h>
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
// MomentumMassBDF2NodeSuppAlg - lumped mass BDF2
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MomentumMassBDF2NodeSuppAlg::MomentumMassBDF2NodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    velocityNm1_(NULL),
    velocityN_(NULL),
    velocityNp1_(NULL),
    densityNm1_(NULL),
    densityN_(NULL),
    densityNp1_(NULL),
    dpdx_(NULL),
    dualNodalVolume_(NULL),
    dt_(0.0),
    nDim_(1),
    gamma1_(0.0),
    gamma2_(0.0),
    gamma3_(0.0)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  VectorFieldType *velocity = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  velocityNm1_ = &(velocity->field_of_state(stk::mesh::StateNM1));
  velocityN_ = &(velocity->field_of_state(stk::mesh::StateN));
  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNm1_ = &(density->field_of_state(stk::mesh::StateNM1));
  densityN_ = &(density->field_of_state(stk::mesh::StateN));
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  dpdx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  nDim_ = meta_data.spatial_dimension();

}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumMassBDF2NodeSuppAlg::setup()
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
MomentumMassBDF2NodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  // deal with lumped mass matrix (diagonal matrix)
  const double *uNm1      =  stk::mesh::field_data(*velocityNm1_, node);
  const double *uN        =  stk::mesh::field_data(*velocityN_, node);
  const double *uNp1      =  stk::mesh::field_data(*velocityNp1_, node);
  const double rhoNm1     = *stk::mesh::field_data(*densityNm1_, node);
  const double rhoN       = *stk::mesh::field_data(*densityN_, node);
  const double rhoNp1     = *stk::mesh::field_data(*densityNp1_, node);
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node);
  const double *dpdx = stk::mesh::field_data(*dpdx_, node);

  const double lhsfac = gamma1_*rhoNp1*dualVolume/dt_;
  const int nDim = nDim_;
  for ( int i = 0; i < nDim; ++i ) {
    rhs[i] += -(gamma1_*rhoNp1*uNp1[i] + gamma2_*rhoN*uN[i] + gamma3_*rhoNm1*uNm1[i])*dualVolume/dt_
                    - dpdx[i]*dualVolume;
    const int row = i*nDim;
    lhs[row+i] += lhsfac;
  }
}

} // namespace nalu
} // namespace Sierra
