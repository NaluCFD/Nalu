/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "nso/MomentumNSOKeElemKernel.h"
#include "AlgTraits.h"
#include "master_element/MasterElement.h"
#include "SolutionOptions.h"
#include "TimeIntegrator.h"

// template and scratch space
#include "BuildTemplates.h"
#include "ScratchViews.h"

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra {
namespace nalu {

template<typename AlgTraits>
MomentumNSOKeElemKernel<AlgTraits>::MomentumNSOKeElemKernel(
  const stk::mesh::BulkData& bulkData,
  SolutionOptions& solnOpts,
  VectorFieldType* ,
  GenericFieldType* Gju,
  const double fourthFac,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    Gju_(Gju),
    lrscv_(sierra::nalu::get_surface_master_element(AlgTraits::topo_)->adjacentNodes()),
    fourthFac_(fourthFac)
{
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  velocityNp1_ = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "velocity");
  densityNp1_ = metaData.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "density");
  pressure_ = metaData.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "pressure");

  if (solnOpts.does_mesh_move())
    velocityRTM_ = metaData.get_field<VectorFieldType>(
      stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = metaData.get_field<VectorFieldType>(
      stk::topology::NODE_RANK, "velocity");

  pressure_ = metaData.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "pressure");

  coordinates_ = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  Gjp_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx");

  MasterElement *meSCS = sierra::nalu::get_surface_master_element(AlgTraits::topo_);
  meSCS->shape_fcn(&v_shape_function_(0,0));

  // add master elements
  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields
  dataPreReqs.add_gathered_nodal_field(*Gju_, AlgTraits::nDim_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*coordinates_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityNp1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityRTM_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*Gjp_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_,1);
  dataPreReqs.add_gathered_nodal_field(*pressure_,1);

  // master element data
  dataPreReqs.add_master_element_call(SCS_AREAV);
  dataPreReqs.add_master_element_call(SCS_GRAD_OP);
  dataPreReqs.add_master_element_call(SCS_GIJ);
}

template<typename AlgTraits>
void
MomentumNSOKeElemKernel<AlgTraits>::execute(
  SharedMemView<double**>& lhs,
  SharedMemView<double *>& rhs,
  stk::mesh::Entity ,
  ScratchViews& scratchViews)
{
  SharedMemView<double***>& v_Gju = scratchViews.get_scratch_view_3D(*Gju_);
  SharedMemView<double**>& v_uNp1 = scratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<double**>& v_velocityRTM = scratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<double**>& v_Gjp = scratchViews.get_scratch_view_2D(*Gjp_);
  SharedMemView<double*>& v_rhoNp1 = scratchViews.get_scratch_view_1D(*densityNp1_);
  SharedMemView<double*>& v_pressure = scratchViews.get_scratch_view_1D(*pressure_);

  SharedMemView<double**>& v_scs_areav = scratchViews.scs_areav;
  SharedMemView<double***>& v_dndx = scratchViews.dndx;
  SharedMemView<double***>& v_gijUpper = scratchViews.gijUpper;
  SharedMemView<double***>& v_gijLower = scratchViews.gijLower;

  // compute nodal ke
  for ( int n = 0; n < AlgTraits::nodesPerElement_; ++n ) {
    double ke = 0.0;
    for ( int j = 0; j < AlgTraits::nDim_; ++j )
      ke += v_uNp1(n,j)*v_uNp1(n,j)/2.0;
    v_ke_(n) = ke;
  }

  for ( int ip = 0; ip < AlgTraits::numScsIp_; ++ip ) {

    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];

    // save off some offsets
    const int ilNdim = il*AlgTraits::nDim_;
    const int irNdim = ir*AlgTraits::nDim_;

    // zero out vector that prevail over all components
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      v_rhoVrtmScs_(i) = 0.0;
      v_uNp1Scs_(i) = 0.0;
      v_dpdxScs_(i) = 0.0;
      v_GjpScs_(i) = 0.0;
      v_dkedxScs_(i) = 0.0;
    }
    double rhoScs = 0.0;

    // determine scs values of interest
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {

      // save off shape function
      const double r = v_shape_function_(ip,ic);

      // compute scs derivatives and flux derivative
      const double pressureIC = v_pressure(ic);
      const double rhoIC = v_rhoNp1(ic);
      const double keIC = v_ke_(ic);
      rhoScs += r*rhoIC;
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        const double dnj = v_dndx(ip,ic,j);
        const double vrtmj = v_velocityRTM(ic,j);
        const double uNp1 = v_uNp1(ic,j);
        const double Gjp = v_Gjp(ic,j);

        v_rhoVrtmScs_(j) += r*rhoIC*vrtmj;
        v_uNp1Scs_(j) += r*uNp1;
        v_dpdxScs_(j) += pressureIC*dnj;
        v_GjpScs_(j) += r*Gjp;
        v_dkedxScs_(j) += keIC*dnj;
      }
    }

    // form ke residual (based on fine scale momentum residual used in Pstab)
    double keResidual = 0.0;
    for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
      keResidual += v_uNp1Scs_(j)*(v_dpdxScs_(j) - v_GjpScs_(j))/rhoScs/2.0;
    }

    // denominator for nu as well as terms for "upwind" nu
    double gUpperMagGradQ = 0.0;
    double rhoVrtmiGLowerRhoVrtmj = 0.0;
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      const double dkedxScsi = v_dkedxScs_(i);
      const double rhoVrtmi = v_rhoVrtmScs_(i);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        gUpperMagGradQ += dkedxScsi*v_gijUpper(ip,i,j)*v_dkedxScs_(j);
        rhoVrtmiGLowerRhoVrtmj += rhoVrtmi*v_gijLower(ip,i,j)*v_rhoVrtmScs_(j);
      }
    }

    // construct nu from ke residual
    const double nuResidual = rhoScs*std::sqrt((keResidual*keResidual)/(gUpperMagGradQ+small_));

    // construct nu from first-order-like approach; SNL-internal write-up (eq 209)
    // for now, only include advection as full set of terms is too diffuse
    const double nuFirstOrder = std::sqrt(rhoVrtmiGLowerRhoVrtmj);

    // limit based on first order; Cupw_ is a fudge factor similar to Guermond's approach
    const double nu = std::min(Cupw_*nuFirstOrder, nuResidual);

    // assemble each component
    for ( int k = 0; k < AlgTraits::nDim_; ++k ) {

      const int indexL = ilNdim + k;
      const int indexR = irNdim + k;

      double gijFac = 0.0;
      for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {

        // save off shape function
        const double r = v_shape_function_(ip,ic);

        // find the row
        const int icNdim = ic*AlgTraits::nDim_;

        // save of some variables
        const double ukNp1 = v_uNp1(ic,k);

        // NSO diffusion-like term; -nu*gUpper*dQ/dxj*ai (residual below)
        double lhsfac = 0.0;
        for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
          const double axi = v_scs_areav(ip,i);
          for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
            const double dnxj = v_dndx(ip,ic,j);
            const double fac = v_gijUpper(ip,i,j)*dnxj*axi;
            const double facGj = r*v_gijUpper(ip,i,j)*v_Gju(ic,k,j)*axi;
            gijFac += fac*ukNp1 - facGj*fourthFac_;
            lhsfac += -fac;
          }
        }

        // no coupling between components
        lhs(indexL,icNdim+k) += nu*lhsfac;
        lhs(indexR,icNdim+k) -= nu*lhsfac;
      }

      // residual; left and right
      const double residualNSO = -nu*gijFac;
      rhs[indexL] -= residualNSO;
      rhs[indexR] += residualNSO;
    }
  }
}

INSTANTIATE_SUPPLEMENTAL_ALGORITHM(MomentumNSOKeElemKernel);

}  // nalu
}  // sierra
