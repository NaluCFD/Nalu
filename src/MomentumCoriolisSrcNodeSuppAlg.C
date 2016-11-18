/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <MomentumCoriolisSrcNodeSuppAlg.h>
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
    nDim_(1),
    earthAngularVelocity_(0.0),
    latitude_(0.0)
{

  pi_ = std::acos(-1.0);

  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  VectorFieldType *velocity = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  // Do we need density for the compressible form of the Coriolis term?
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  dualNodalVolume_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  // extract user parameters from solution options
  earthAngularVelocity_ = realm_.solutionOptions_->earthAngularVelocity_;
  latitude_ = realm_.solutionOptions_->latitude_ * pi_ / 180.0;
  nDim_ = meta_data.spatial_dimension();
  if (nDim_ != 3) 
    throw std::runtime_error("MomentumCoriolisSrcNodeSuppAlg: nDim != 3");
  eastVector_.resize(nDim_);
  northVector_.resize(nDim_);
  upVector_.resize(nDim_);
  eastVector_ = realm_.solutionOptions_->eastVector_;
  northVector_ = realm_.solutionOptions_->northVector_;

  // normalize the east and north vectors
  double magE = std::sqrt(eastVector_[0]*eastVector_[0]+eastVector_[1]*eastVector_[1]+eastVector_[2]*eastVector_[2]);
  double magN = std::sqrt(northVector_[0]*northVector_[0]+northVector_[1]*northVector_[1]+northVector_[2]*northVector_[2]);
  for (int i=0; i<nDim_; ++i) {
    eastVector_[i] /= magE;
    northVector_[i] /= magN;
  }

  // calculate the 'up' unit vector
  cross_product(eastVector_, northVector_, upVector_);

  // some factors that do not change
  sinphi_ = std::sin(latitude_);
  cosphi_ = std::cos(latitude_);
  corfac_ = 2.0*earthAngularVelocity_;

  // Jacobian entries
  Jxy_ = corfac_ * (eastVector_[0]*northVector_[1]-northVector_[0]*eastVector_[1])*sinphi_
                 + (upVector_[0]*eastVector_[1]-eastVector_[0]*upVector_[1])*cosphi_;
  Jxz_ = corfac_ * (eastVector_[0]*northVector_[2]-northVector_[0]*eastVector_[2])*sinphi_
                 + (upVector_[0]*eastVector_[2]-eastVector_[0]*upVector_[2])*cosphi_;
  Jyz_ = corfac_ * (eastVector_[1]*northVector_[2]-northVector_[1]*eastVector_[2])*sinphi_
                 + (upVector_[1]*eastVector_[2]-eastVector_[1]*upVector_[2])*cosphi_;

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
  const double ue = eastVector_[0]*uNp1[0] + eastVector_[1]*uNp1[1] + eastVector_[2]*uNp1[2];
  const double un = northVector_[0]*uNp1[0] + northVector_[1]*uNp1[1] + northVector_[2]*uNp1[2];
  const double uu = upVector_[0]*uNp1[0] + upVector_[1]*uNp1[1] + upVector_[2]*uNp1[2];

  // calculate acceleration in east-north-up coordinates
  const double ae = corfac_ * (un*sinphi_ - uu*cosphi_);
  const double an = -corfac_*ue*sinphi_;
  const double au = corfac_*ue*cosphi_;
  
  // calculate acceleration in model x-y-z coordinates
  const double ax = ae*eastVector_[0] + an*northVector_[0] + au*upVector_[0];
  const double ay = ae*eastVector_[1] + an*northVector_[1] + au*upVector_[1];
  const double az = ae*eastVector_[2] + an*northVector_[2] + au*upVector_[2];

  const double rhoNp1     = *stk::mesh::field_data(*densityNp1_, node );
  const double dualVolume = *stk::mesh::field_data(*dualNodalVolume_, node );

  const double fac2 = rhoNp1*dualVolume;
  rhs[0] += fac2*ax;
  rhs[1] += fac2*ay;
  rhs[2] += fac2*az;

  // Only the off-diagonal LHS entries are non-zero
  lhs[1] += fac2*Jxy_;
  lhs[2] += fac2*Jxz_;
  lhs[3] -= fac2*Jxy_; // Jyx = - Jxy
  lhs[5] += fac2*Jyz_;
  lhs[6] -= fac2*Jxz_; // Jzx = - Jxz
  lhs[7] -= fac2*Jyz_; // Jzy = - Jyz

}

//--------------------------------------------------------------------------
//-------- cross_product ----------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumCoriolisSrcNodeSuppAlg::cross_product(
  std::vector<double> u, std::vector<double> v, std::vector<double> cross)
{
  cross[0] =   u[1]*v[2] - u[2]*v[1];
  cross[1] = -(u[0]*v[2] - u[2]*v[0]);
  cross[2] =   u[0]*v[1] - u[1]*v[0];
}

} // namespace nalu
} // namespace Sierra
