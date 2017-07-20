/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "nso/ScalarNSOElemKernel.h"
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
ScalarNSOElemKernel<AlgTraits>::ScalarNSOElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ScalarFieldType* scalarQ,
  VectorFieldType* Gjq,
  ScalarFieldType* diffFluxCoeff,
  const double fourthFac,
  const double altResFac,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    diffFluxCoeff_(diffFluxCoeff),
    Gjq_(Gjq),
    lrscv_(sierra::nalu::get_surface_master_element(AlgTraits::topo_)->adjacentNodes()),
    fourthFac_(fourthFac),
    altResFac_(altResFac),
    om_altResFac_(1.0 - altResFac),
    shiftedGradOp_(solnOpts.get_shifted_grad_op(scalarQ->name()))
{
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  ScalarFieldType *density = metaData.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "density");

  scalarQN_ = &(scalarQ->field_of_state(stk::mesh::StateN));
  scalarQNp1_ = &(scalarQ->field_of_state(stk::mesh::StateNP1));
  if (scalarQ->number_of_states() == 2)
    scalarQNm1_ = scalarQN_;
  else
    scalarQNm1_ = &(scalarQ->field_of_state(stk::mesh::StateNM1));

  densityN_ = &(density->field_of_state(stk::mesh::StateN));
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  if (density->number_of_states() == 2)
    densityNm1_ = densityN_;
  else
    densityNm1_ = &(density->field_of_state(stk::mesh::StateNM1));

  if (solnOpts.does_mesh_move())
    velocityRTM_ = metaData.get_field<VectorFieldType>(
      stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = metaData.get_field<VectorFieldType>(
      stk::topology::NODE_RANK, "velocity");

  coordinates_ = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  MasterElement *meSCS = sierra::nalu::get_surface_master_element(AlgTraits::topo_);
  meSCS->shape_fcn(&v_shape_function_(0,0));

  // add master elements
  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*velocityRTM_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*Gjq_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*scalarQNm1_, 1);
  dataPreReqs.add_gathered_nodal_field(*scalarQN_, 1);
  dataPreReqs.add_gathered_nodal_field(*scalarQNp1_, 1);

  dataPreReqs.add_gathered_nodal_field(*densityNm1_,1);
  dataPreReqs.add_gathered_nodal_field(*densityN_,1);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_,1);
  dataPreReqs.add_gathered_nodal_field(*diffFluxCoeff_,1);

  // master element data
  dataPreReqs.add_master_element_call(SCS_AREAV, CURRENT_COORDINATES);
  if ( shiftedGradOp_ )
    dataPreReqs.add_master_element_call(SCS_SHIFTED_GRAD_OP, CURRENT_COORDINATES);
  else
    dataPreReqs.add_master_element_call(SCS_GRAD_OP, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(SCS_GIJ, CURRENT_COORDINATES);
}

template<typename AlgTraits>
void
ScalarNSOElemKernel<AlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
  dt_ = timeIntegrator.get_time_step();
  gamma1_ = timeIntegrator.get_gamma1();
  gamma2_ = timeIntegrator.get_gamma2();
  gamma3_ = timeIntegrator.get_gamma3(); // gamma3 may be zero
}

template<typename AlgTraits>
void
ScalarNSOElemKernel<AlgTraits>::execute(
  SharedMemView<double**>& lhs,
  SharedMemView<double *>& rhs,
  ScratchViews& scratchViews)
{
  double w_dqdxScs    [AlgTraits::nDim_];
  double w_rhoVrtmScs [AlgTraits::nDim_];

  SharedMemView<double**>& v_Gjq = scratchViews.get_scratch_view_2D(*Gjq_);
  SharedMemView<double**>& v_velocityRTM = scratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<double*>& v_qNm1 = scratchViews.get_scratch_view_1D(*scalarQNm1_);
  SharedMemView<double*>& v_qN = scratchViews.get_scratch_view_1D(*scalarQN_);
  SharedMemView<double*>& v_qNp1 = scratchViews.get_scratch_view_1D(*scalarQNp1_);
  SharedMemView<double*>& v_rhoNm1 = scratchViews.get_scratch_view_1D(*densityNm1_);
  SharedMemView<double*>& v_rhoN = scratchViews.get_scratch_view_1D(*densityN_);
  SharedMemView<double*>& v_rhoNp1 = scratchViews.get_scratch_view_1D(*densityNp1_);
  SharedMemView<double*>& v_diffFluxCoeff = scratchViews.get_scratch_view_1D(*diffFluxCoeff_);

  SharedMemView<double**>& v_scs_areav = scratchViews.get_me_views(CURRENT_COORDINATES).scs_areav;
  SharedMemView<double***>& v_dndx = shiftedGradOp_ 
    ? scratchViews.get_me_views(CURRENT_COORDINATES).dndx_shifted 
    : scratchViews.get_me_views(CURRENT_COORDINATES).dndx;
  SharedMemView<double***>& v_gijUpper = scratchViews.get_me_views(CURRENT_COORDINATES).gijUpper;
  SharedMemView<double***>& v_gijLower = scratchViews.get_me_views(CURRENT_COORDINATES).gijLower;

  for ( int ip = 0; ip < AlgTraits::numScsIp_; ++ip ) {

    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];

    // zero out; scalar
    double qNm1Scs = 0.0;
    double qNScs = 0.0;
    double qNp1Scs = 0.0;
    double rhoNm1Scs = 0.0;
    double rhoNScs = 0.0;
    double rhoNp1Scs = 0.0;
    double dFdxAdv = 0.0;
    double dFdxDiff = 0.0;
    double dFdxCont = 0.0;

    // zero out vector
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      w_dqdxScs[i] = 0.0;
      w_rhoVrtmScs[i] = 0.0;
    }

    // determine scs values of interest
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      // save off shape function
      const double r = v_shape_function_(ip,ic);

      // time term; scalar q
      qNm1Scs += r*v_qNm1(ic);
      qNScs += r*v_qN(ic);
      qNp1Scs += r*v_qNp1(ic);

      // time term, density
      rhoNm1Scs += r*v_rhoNm1(ic);
      rhoNScs += r*v_rhoN(ic);
      rhoNp1Scs += r*v_rhoNp1(ic);

      // compute scs derivatives and flux derivative
      const double qIC = v_qNp1(ic);
      const double rhoIC = v_rhoNp1(ic);
      const double diffFluxCoeffIC = v_diffFluxCoeff(ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        const double dnj = v_dndx(ip,ic,j);
        const double vrtmj = v_velocityRTM(ic,j);
        w_dqdxScs[j] += qIC*dnj;
        w_rhoVrtmScs[j] += r*rhoIC*vrtmj;
        dFdxAdv += rhoIC*vrtmj*qIC*dnj;
        dFdxDiff += diffFluxCoeffIC*v_Gjq(ic,j)*dnj;
        dFdxCont += rhoIC*vrtmj*dnj;
      }
    }

    // full continuity residual
    const double contRes = (gamma1_*rhoNp1Scs + gamma2_*rhoNScs + gamma3_*rhoNm1Scs)/dt_ + dFdxCont;

    // compute residual for NSO; linearized first
    double residualAlt = dFdxAdv - qNp1Scs*dFdxCont;
    for ( int j = 0; j < AlgTraits::nDim_; ++j )
      residualAlt -= w_rhoVrtmScs[j]*w_dqdxScs[j];

    // compute residual for NSO; pde-based second
    const double time = (gamma1_*rhoNp1Scs*qNp1Scs + gamma2_*rhoNScs*qNScs + gamma3_*rhoNm1Scs*qNm1Scs)/dt_;
    const double residualPde = time + dFdxAdv - dFdxDiff - contRes*qNp1Scs*nonConservedForm_;

    // final form
    const double residual = residualAlt*altResFac_ + residualPde*om_altResFac_;

    // denominator for nu as well as terms for "upwind" nu
    double gUpperMagGradQ = 0.0;
    double rhoVrtmiGLowerRhoVrtmj = 0.0;
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      const double dqdxScsi = w_dqdxScs[i];
      const double rhoVrtmi = w_rhoVrtmScs[i];
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        gUpperMagGradQ += dqdxScsi*v_gijUpper(ip,i,j)*w_dqdxScs[j];
        rhoVrtmiGLowerRhoVrtmj += rhoVrtmi*v_gijLower(ip,i,j)*w_rhoVrtmScs[j];
      }
    }

    // construct nu from residual
    const double nuResidual = std::sqrt((residual*residual)/(gUpperMagGradQ+small_));

    // construct nu from first-order-like approach; SNL-internal write-up (eq 209)
    // for now, only include advection as full set of terms is too diffuse
    const double nuFirstOrder = std::sqrt(rhoVrtmiGLowerRhoVrtmj);

    // limit based on first order; Cupw_ is a fudge factor similar to Guermond's approach
    const double nu = std::min(Cupw_*nuFirstOrder, nuResidual);

    double gijFac = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {

      // save off shape function
      const double r = v_shape_function_(ip,ic);

      // save of some variables
      const double qIC = v_qNp1(ic);

      // NSO diffusion-like term; -nu*gUpper*(dQ/dxj - Gjq)*ai (residual below)
      double lhsfac = 0.0;
      for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
        const double axi = v_scs_areav(ip,i);
        for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
          const double dnxj = v_dndx(ip,ic,j);
          const double fac = v_gijUpper(ip,i,j)*dnxj*axi;
          const double facGj = r*v_gijUpper(ip,i,j)*v_Gjq(ic,j)*axi;
          gijFac += fac*qIC - facGj*fourthFac_;
          lhsfac += -fac;
        }
      }

      lhs(il,ic) += nu*lhsfac;
      lhs(ir,ic) -= nu*lhsfac;
    }

    // residual; left and right
    const double residualNSO = -nu*gijFac;
    rhs(il) -= residualNSO;
    rhs(ir) += residualNSO;
  }
}

INSTANTIATE_KERNEL(ScalarNSOElemKernel);

}  // nalu
}  // sierra
