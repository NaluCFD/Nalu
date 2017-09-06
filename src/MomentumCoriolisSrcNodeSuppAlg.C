/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <MomentumCoriolisSrcNodeSuppAlg.h>
#include <CoriolisSrc.h>
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
// MomentumCoriolisSrcNodeSuppAlg - base class for algorithm
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MomentumCoriolisSrcNodeSuppAlg::MomentumCoriolisSrcNodeSuppAlg(
  Realm &realm)
  : SupplementalAlgorithm(realm),
    densityNp1_(NULL),
    velocityNp1_(NULL),
    dualNodalVolume_(NULL),
    cor_(*realm.solutionOptions_)
{
  
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  VectorFieldType *velocity = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

}

//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumCoriolisSrcNodeSuppAlg::setup()
{
  // all set up in constructor
}

//--------------------------------------------------------------------------
//-------- node_execute ----------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumCoriolisSrcNodeSuppAlg::node_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity node)
{
  const double *uNp1 = stk::mesh::field_data(*velocityNp1_, node );


  // calculate the velocity vector in east-north-up coordinates
  const double ue = cor_.eastVector_[0]*uNp1[0] + cor_.eastVector_[1]*uNp1[1] + cor_.eastVector_[2]*uNp1[2];
  const double un = cor_.northVector_[0]*uNp1[0] + cor_.northVector_[1]*uNp1[1] + cor_.northVector_[2]*uNp1[2];
  const double uu = cor_.upVector_[0]*uNp1[0] + cor_.upVector_[1]*uNp1[1] + cor_.upVector_[2]*uNp1[2];

  // calculate acceleration in east-north-up coordinates
  const double ae = cor_.corfac_ * (un*cor_.sinphi_ - uu*cor_.cosphi_);
  const double an = -cor_.corfac_*ue*cor_.sinphi_;
  const double au = cor_.corfac_*ue*cor_.cosphi_;
  
  // calculate acceleration in model x-y-z coordinates
  const double ax = ae*cor_.eastVector_[0] + an*cor_.northVector_[0] + au*cor_.upVector_[0];
  const double ay = ae*cor_.eastVector_[1] + an*cor_.northVector_[1] + au*cor_.upVector_[1];
  const double az = ae*cor_.eastVector_[2] + an*cor_.northVector_[2] + au*cor_.upVector_[2];

  const double rhoNp1     = *stk::mesh::field_data(*densityNp1_, node );
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );

  const double fac2 = rhoNp1*dualVolume;
  rhs[0] += fac2*ax;
  rhs[1] += fac2*ay;
  rhs[2] += fac2*az;

  // Only the off-diagonal LHS entries are non-zero
  lhs[1] += fac2*cor_.Jxy_;
  lhs[2] += fac2*cor_.Jxz_;
  lhs[3] -= fac2*cor_.Jxy_; // Jyx = - Jxy
  lhs[5] += fac2*cor_.Jyz_;
  lhs[6] -= fac2*cor_.Jxz_; // Jzx = - Jxz
  lhs[7] -= fac2*cor_.Jyz_; // Jzy = - Jyz

}

} // namespace nalu
} // namespace Sierra
