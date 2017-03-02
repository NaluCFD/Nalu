/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <ContinuityAdvElemSuppAlg.h>
#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <master_element/MasterElement.h>

// template and scratch space
#include <BuildTemplates.h>
#include <ScratchViews.h>

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
// ContinuityAdvElemSuppAlg - CMM (BDF2) for continuity equation ()p-dof)
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
template<typename AlgTraits>
ContinuityAdvElemSuppAlg<AlgTraits>::ContinuityAdvElemSuppAlg(
  Realm &realm,
  ElemDataRequests& dataPreReqs)
  : SupplementalAlgorithm(realm),
    bulkData_(&realm.bulk_data()),
    velocityRTM_(NULL),
    Gpdx_(NULL),
    pressure_(NULL),
    densityNp1_(NULL),
    coordinates_(NULL),
    projTimeScale_(1.0),
    meshMotion_(realm_.does_mesh_move()),
    shiftMdot_(realm_.get_cvfem_shifted_mdot()),
    shiftPoisson_(realm_.get_cvfem_shifted_poisson()),
    reducedSensitivities_(realm_.get_cvfem_reduced_sens_poisson()),
    interpTogether_(realm_.get_mdot_interp()),
    om_interpTogether_(1.0-interpTogether_),
    lrscv_(realm.get_surface_master_element(AlgTraits::topo_)->adjacentNodes())
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  if ( meshMotion_ )
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  Gpdx_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx");
  pressure_ = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure");
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  coordinates_ = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  MasterElement *meSCS = realm.get_surface_master_element(AlgTraits::topo_);
  if ( shiftMdot_ )
    meSCS->shifted_shape_fcn(&v_shape_function_(0,0));
  else
    meSCS->shape_fcn(&v_shape_function_(0,0));

  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields and data
  dataPreReqs.add_gathered_nodal_field(*coordinates_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityRTM_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*pressure_, 1);
  dataPreReqs.add_gathered_nodal_field(*Gpdx_, AlgTraits::nDim_);
  dataPreReqs.add_master_element_call(SCS_AREAV);

  // TODO: Integrate checks for shiftPoisson and reducedSensitivities
  dataPreReqs.add_master_element_call(SCS_GRAD_OP);
  dataPreReqs.add_master_element_call(SCS_SHIFTED_GRAD_OP);

  // Implementation details: code is designed to manage the following
  // When shiftPoisson_ is TRUE, reducedSensitivities_ is enforced to be TRUE
  // However, shiftPoisson_ can be FALSE while reducedSensitivities_ is TRUE
}


//--------------------------------------------------------------------------
//-------- setup -----------------------------------------------------------
//--------------------------------------------------------------------------
template<typename AlgTraits>
void
ContinuityAdvElemSuppAlg<AlgTraits>::setup()
{
  const double dt = realm_.get_time_step();
  const double gamma1 = realm_.get_gamma1();
  projTimeScale_ = dt/gamma1;
}

//--------------------------------------------------------------------------
//-------- elem_execute ----------------------------------------------------
//--------------------------------------------------------------------------
template<typename AlgTraits>
void
ContinuityAdvElemSuppAlg<AlgTraits>::element_execute(
  double *lhs,
  double *rhs,
  stk::mesh::Entity /*element*/,
  ScratchViews& scratchViews)
{
  SharedMemView<double*>& densityNp1_view = scratchViews.get_scratch_view_1D(*densityNp1_);
  SharedMemView<double*>& pressure_view = scratchViews.get_scratch_view_1D(*pressure_);

  SharedMemView<double**>& velocity_view = scratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<double**>& Gpdx_view = scratchViews.get_scratch_view_2D(*Gpdx_);

  SharedMemView<double**>& scs_areav = scratchViews.scs_areav;
  // SharedMemView<double***>& dndx = scratchViews.dndx;
  // SharedMemView<double***>& dndx_shifted = scratchViews.dndx_shifted;

  SharedMemView<double***>& dndx = shiftPoisson_?
    scratchViews.dndx_shifted : scratchViews.dndx;
  SharedMemView<double***>& dndx_lhs = (!shiftPoisson_ && reducedSensitivities_)?
    scratchViews.dndx_shifted : scratchViews.dndx;

  for (int ip=0; ip < AlgTraits::numScsIp_; ++ip) {
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];

    const int rowL = il * AlgTraits::nodesPerElement_;
    const int rowR = ir * AlgTraits::nodesPerElement_;

    double rhoIp = 0.0;
    for (int j=0; j<AlgTraits::nDim_; ++j) {
      v_uIp_(j) = 0.0;
      v_rho_uIp_(j) = 0.0;
      v_Gpdx_Ip_(j) = 0.0;
      v_dpdxIp_(j) = 0.0;
    }

    for (int ic=0; ic < AlgTraits::nodesPerElement_; ++ic) {
      const double r = v_shape_function_(ip, ic);
      const double nodalPressure = pressure_view(ic);
      const double nodalRho = densityNp1_view(ic);

      rhoIp += r * nodalRho;
      double lhsfac = 0.0;

      for (int j=0; j<AlgTraits::nDim_; j++) {
        v_Gpdx_Ip_(j) += r * Gpdx_view(ic, j);
        v_uIp_(j)     += r * velocity_view(ic, j);
        v_rho_uIp_(j) += r * nodalRho * velocity_view(ic, j);
        v_dpdxIp_(j)  += dndx(ip, ic, j) * nodalPressure;
        lhsfac += -dndx_lhs(ip, ic, j) * scs_areav(ip, j);
      }

      lhs[rowL + ic] += lhsfac;
      lhs[rowR + ic] -= lhsfac;
    }

    // assemble mdot
    double mdot = 0.0;
    for (int j=0; j<AlgTraits::nDim_; j++) {
      mdot += (interpTogether_ * v_rho_uIp_(j) +
               om_interpTogether_ * rhoIp * v_uIp_(j) -
               projTimeScale_ * ( v_dpdxIp_(j) - v_Gpdx_Ip_(j))) * scs_areav(ip,j);
    }

    // residuals
    rhs[il] -= mdot / projTimeScale_;
    rhs[ir] += mdot / projTimeScale_;
  }
}

INSTANTIATE_SUPPLEMENTAL_ALGORITHM(ContinuityAdvElemSuppAlg);

} // namespace nalu
} // namespace Sierra
