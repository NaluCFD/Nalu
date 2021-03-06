/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/VolumeOfFluidSucvNsoElemKernel.h"
#include "AlgTraits.h"
#include "master_element/MasterElement.h"
#include "TimeIntegrator.h"
#include "SolutionOptions.h"

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
VolumeOfFluidSucvNsoElemKernel<AlgTraits>::VolumeOfFluidSucvNsoElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ScalarFieldType* vof,
  const double sucvFac,
  const double nsoFac,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    sucvFac_(sucvFac),
    nsoFac_(nsoFac),
    lrscv_(sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_)->adjacentNodes())
{
  // save off fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  
  vofN_ = &(vof->field_of_state(stk::mesh::StateN));
  vofNp1_ = &(vof->field_of_state(stk::mesh::StateNP1));
  if (vof->number_of_states() == 2)
    vofNm1_ = vofN_;
  else
    vofNm1_ = &(vof->field_of_state(stk::mesh::StateNM1));
  
  std::string velocity_name = solnOpts.does_mesh_move() ? "velocity_rtm" : "velocity";
  velocityRTM_ = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, velocity_name);
  coordinates_ = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_);

  // compute shape function
  get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCS->shape_fcn(ptr);}, v_shape_function_);

  // add master elements
  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*velocityRTM_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*vofNm1_, 1);
  dataPreReqs.add_gathered_nodal_field(*vofN_, 1);
  dataPreReqs.add_gathered_nodal_field(*vofNp1_, 1);

  // master element data
  dataPreReqs.add_master_element_call(SCS_AREAV, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(SCS_GRAD_OP, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(SCS_GIJ, CURRENT_COORDINATES);
  
  // correction to gij
  gijFac_ = (AlgTraits::nodesPerElement_ == 4 || AlgTraits::nodesPerElement_ == 8) ? 4.0 : 1.0;
}

template<typename AlgTraits>
VolumeOfFluidSucvNsoElemKernel<AlgTraits>::~VolumeOfFluidSucvNsoElemKernel()
{}

template<typename AlgTraits>
void
VolumeOfFluidSucvNsoElemKernel<AlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
  dt_ = timeIntegrator.get_time_step();
  gamma1_ = timeIntegrator.get_gamma1();
  gamma2_ = timeIntegrator.get_gamma2();
  gamma3_ = timeIntegrator.get_gamma3(); // gamma3 may be zero
}

template<typename AlgTraits>
void
VolumeOfFluidSucvNsoElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>& lhs,
  SharedMemView<DoubleType *>&rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  NALU_ALIGNED DoubleType w_uIp[AlgTraits::nDim_];
  NALU_ALIGNED DoubleType w_dvofdxIp[AlgTraits::nDim_];

  SharedMemView<DoubleType*>& v_vofNm1 = scratchViews.get_scratch_view_1D(
    *vofNm1_);
  SharedMemView<DoubleType*>& v_vofN = scratchViews.get_scratch_view_1D(
    *vofN_);
  SharedMemView<DoubleType*>& v_vofNp1 = scratchViews.get_scratch_view_1D(
    *vofNp1_);
  SharedMemView<DoubleType**>& v_vrtm = scratchViews.get_scratch_view_2D(*velocityRTM_);

  SharedMemView<DoubleType**>& v_scs_areav = scratchViews.get_me_views(CURRENT_COORDINATES).scs_areav;
  SharedMemView<DoubleType***>& v_dndx = scratchViews.get_me_views(CURRENT_COORDINATES).dndx;
  SharedMemView<DoubleType***>& v_gijUpper = scratchViews.get_me_views(CURRENT_COORDINATES).gijUpper;
  SharedMemView<DoubleType***>& v_gijLower = scratchViews.get_me_views(CURRENT_COORDINATES).gijLower;

  for ( int ip = 0; ip < AlgTraits::numScsIp_; ++ip ) {

    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];
    
    // zero out; scalar
    DoubleType vofNm1Ip = 0.0;
    DoubleType vofNIp = 0.0;
    DoubleType vofNp1Ip = 0.0;
    
    // zero out; vector
    for (int j = 0; j < AlgTraits::nDim_; ++j) {
      w_uIp[j] = 0.0;
      w_dvofdxIp[j] = 0.0;
    }

    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      // save off shape function
      const DoubleType r = v_shape_function_(ip,ic);
      vofNm1Ip += r*v_vofNm1(ic);
      vofNIp += r*v_vofN(ic);
      vofNp1Ip += r*v_vofNp1(ic);
      const DoubleType vofIc = v_vofNp1(ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        w_uIp[j] += r*v_vrtm(ic,j);
        w_dvofdxIp[j] += vofIc*v_dndx(ip,ic,j);
      }
    }

    // form uj*nj; fold in uigijuj
    DoubleType ujaj = 0.0;
    DoubleType uigijuj = 0.0;
    DoubleType gUpperMagGradQ = 0.0;
    DoubleType uigLijuj = 0.0;
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      const DoubleType ai = v_scs_areav(ip,i);
      const DoubleType ui = w_uIp[i];
      const DoubleType dvofdxIpi = w_dvofdxIp[i];

      ujaj += ui*ai;
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        uigijuj += ui*v_gijLower(ip,i,j)*w_uIp[j];
        uigLijuj += ui*v_gijLower(ip,i,j)*w_uIp[j];
        gUpperMagGradQ += dvofdxIpi*v_gijUpper(ip,i,j)*w_dvofdxIp[j];
      }
    }

    // construct tau (from flow-aligned time scale) and nu (from residual)
    const DoubleType tau = 1.0/(stk::math::sqrt((2.0/dt_)*(2.0/dt_) + gijFac_*uigijuj));
    
    // form mass term residual
    const DoubleType mass = (gamma1_*vofNp1Ip + gamma2_*vofNIp + gamma3_*vofNm1Ip)/dt_;

    DoubleType uidVofdxi = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      // save of some variables
      const DoubleType vofIc = v_vofNp1(ic);
      
      // SUCV diffusion-like term; -tau*uidvof/dxi*ujnj (residual below)
      DoubleType lhsfac = 0.0;
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        const DoubleType ujdNj = w_uIp[j]*v_dndx(ip,ic,j);
        uidVofdxi += ujdNj*vofIc;
        lhsfac += -ujdNj;
      }
      lhs(il,ic) += tau*ujaj*lhsfac*sucvFac_;
      lhs(ir,ic) -= tau*ujaj*lhsfac*sucvFac_;	  
    }
    
    // form residual
    const DoubleType residual = mass+uidVofdxi;

    // full sucv residual
    const DoubleType residualSucv = -tau*ujaj*residual*sucvFac_;
    
    // residual; left and right
    rhs(il) -= residualSucv;
    rhs(ir) += residualSucv; 

    // construct and nu (from residual)
    const DoubleType nuResidual = stk::math::sqrt((residual*residual)/(gUpperMagGradQ/gijFac_+small_));

    // construct nu from first-order-like approach; SNL-internal write-up (eq 209)
    const DoubleType nuFirstOrder = stk::math::sqrt(uigLijuj);

    // limit based on first order; Cupw_ is a fudge factor similar to Guermond's approach
    const DoubleType nu = stk::math::min(Cupw_*nuFirstOrder, nuResidual);

    DoubleType gijFac = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {

      // save of some variables
      const DoubleType vofIc = v_vofNp1(ic);

      // NSO diffusion-like term; -nu*gUpper*dvof/dxj*ai (residual below)
      DoubleType lhsfac = 0.0;
      for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
        const DoubleType axi = v_scs_areav(ip,i);
        for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
          const DoubleType dnxj = v_dndx(ip,ic,j);
          const DoubleType fac = v_gijUpper(ip,i,j)*dnxj*axi;
          gijFac += fac*vofIc;
          lhsfac += -fac;
        }
      }
      
      lhs(il,ic) += nu*lhsfac;
      lhs(ir,ic) -= nu*lhsfac;
    }

    // residual; left and right
    const DoubleType residualNSO = -nu*gijFac;
    rhs(il) -= residualNSO;
    rhs(ir) += residualNSO;   
  }
  
}

INSTANTIATE_KERNEL(VolumeOfFluidSucvNsoElemKernel);

}  // nalu
}  // sierra
