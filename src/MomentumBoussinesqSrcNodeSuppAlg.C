/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <MomentumBoussinesqSrcNodeSuppAlg.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <SolutionOptions.h>
#include <SupplementalAlgorithm.h>

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
// MomentumBoussinesqSrcNodeSuppAlg - -rho*beta*(T-Tref) g
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MomentumBoussinesqSrcNodeSuppAlg::MomentumBoussinesqSrcNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    temperature_(NULL),
    dualNodalVolume_(NULL),
    tRef_(298.0),
    rhoRef_(1.0),
    beta_(1.0),
    nDim_(1)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  temperature_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "temperature");
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  // extract user parameters from solution options
  tRef_ = realm_.solutionOptions_->referenceTemperature_;
  rhoRef_ = realm_.solutionOptions_->referenceDensity_;
  beta_ = realm_.solutionOptions_->thermalExpansionCoeff_;
  nDim_ = meta_data.spatial_dimension();
  gravity_.resize(nDim_);
  gravity_ = realm_.solutionOptions_->gravity_;
}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumBoussinesqSrcNodeSuppAlg::setup()
{
  // all set up in constructor
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumBoussinesqSrcNodeSuppAlg::node_execute(
  double */*lhs*/,
  double *rhs,
  stk::mesh::Entity node)
{
  // no lhs contribution; all rhs source term; density should be constant...
  const double temperature = *stk::mesh::field_data(*temperature_, node );
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );
  const double fac = -rhoRef_*beta_*(temperature - tRef_)*dualVolume;
  for ( int i = 0; i < nDim_; ++i ) {
    rhs[i] += fac*gravity_[i];
  }
}

} // namespace nalu
} // namespace Sierra
