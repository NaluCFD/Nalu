/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/VolumeOfFluidScsUpwAdvElemKernel.h"
#include "AlgTraits.h"
#include "master_element/MasterElement.h"
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
VolumeOfFluidScsUpwAdvElemKernel<AlgTraits>::VolumeOfFluidScsUpwAdvElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ScalarFieldType* vof,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    hoUpwind_(solnOpts.get_upw_factor(vof->name())),
    useLimiter_(solnOpts.primitive_uses_limiter(vof->name())),
    lrscv_(sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_)->adjacentNodes())
{
  // save off fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  
  vofNp1_ = &(vof->field_of_state(stk::mesh::StateNP1));
  dvofdx_ = metaData.get_field<double>(stk::topology::NODE_RANK, "dvofdx");

  coordinates_ = metaData.get_field<double>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  volumeFlowRate_ = metaData.get_field<double>(stk::topology::ELEMENT_RANK, "volume_flow_rate_scs");

  MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_);

  // compute shape function
  get_scs_shape_fn_data<AlgTraits>([&](double* ptr){meSCS->shape_fcn(ptr);}, v_shape_function_);

  // add master elements
  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*vofNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*dvofdx_, AlgTraits::nDim_);
  dataPreReqs.add_element_field(*volumeFlowRate_, AlgTraits::numScsIp_);

  NaluEnv::self().naluOutputP0() << "VolumeOfFluidScsUpwAdvElemKernel hoUpwind/limit: " 
                                 << hoUpwind_ << "/" << useLimiter_ << std::endl;
}

template<typename AlgTraits>
VolumeOfFluidScsUpwAdvElemKernel<AlgTraits>::~VolumeOfFluidScsUpwAdvElemKernel()
{}

template<typename AlgTraits>
void
VolumeOfFluidScsUpwAdvElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>& lhs,
  SharedMemView<DoubleType *>&rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  DoubleType w_coordIp[AlgTraits::nDim_];

  SharedMemView<DoubleType**>& v_coordinates = scratchViews.get_scratch_view_2D(*coordinates_);
  SharedMemView<DoubleType*>& v_vofNp1 = scratchViews.get_scratch_view_1D(*vofNp1_);
  SharedMemView<DoubleType**>& v_dvofdx = scratchViews.get_scratch_view_2D(*dvofdx_);
  SharedMemView<DoubleType*>& v_vdot = scratchViews.get_scratch_view_1D(*volumeFlowRate_);

  //============================================
  // SCS contribution; vdot*alpha*nj*dS
  //============================================
  for ( int ip = 0; ip < AlgTraits::numScsIp_; ++ip ) {
    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];

    // save off vdot
    const DoubleType vdot = v_vdot(ip);

    // zero out values of interest for this ip
    for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
      w_coordIp[j] = 0.0;
    }

    // compute ip property and
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const DoubleType r = v_shape_function_(ip,ic);
      for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
        w_coordIp[i] += r*v_coordinates(ic,i);
      }
    }
    
    // left and right extrapolation
    DoubleType dqL = 0.0;
    DoubleType dqR = 0.0;
    for(int j = 0; j < AlgTraits::nDim_; ++j ) {
      const DoubleType dxjL = w_coordIp[j] - v_coordinates(il,j);
      const DoubleType dxjR = v_coordinates(ir,j) - w_coordIp[j];
      dqL += dxjL*v_dvofdx(il,j);
      dqR += dxjR*v_dvofdx(ir,j);
    }

    // add limiter if appropriate
    DoubleType limitL = 1.0;
    DoubleType limitR = 1.0;
    if ( useLimiter_ ) {
      const DoubleType dq = v_vofNp1(ir) - v_vofNp1(il);
      const DoubleType dqMl = 2.0*2.0*dqL - dq;
      const DoubleType dqMr = 2.0*2.0*dqR - dq;
      limitL = van_leer(dqMl, dq);
      limitR = van_leer(dqMr, dq);
    }

    // extrapolated; for now limit (along edge is fine)
    const DoubleType vofIpL = v_vofNp1(il) + dqL*hoUpwind_*limitL;
    const DoubleType vofIpR = v_vofNp1(ir) - dqR*hoUpwind_*limitR;

    // upwind
    const DoubleType vofUpwind = stk::math::if_then_else(vdot > 0,
                                                         vofIpL,
                                                         vofIpR);

    // advection
    const DoubleType adv = vdot*vofUpwind;

    // right hand side; L and R
    rhs(il) -= adv;
    rhs(ir) += adv;

    // upwind; left node
    const DoubleType alhsfacL = 0.5*(vdot+stk::math::abs(vdot));
    lhs(il,il) += alhsfacL;
    lhs(ir,il) -= alhsfacL;

    // upwind; right node
    const DoubleType alhsfacR = 0.5*(vdot-stk::math::abs(vdot));
    lhs(il,il) -= alhsfacR;
    lhs(ir,il) += alhsfacR;
  }
}

template<class AlgTraits>
DoubleType
VolumeOfFluidScsUpwAdvElemKernel<AlgTraits>::van_leer(
  const DoubleType &dqm,
  const DoubleType &dqp)
{
  DoubleType limit = (2.0*(dqm*dqp+stk::math::abs(dqm*dqp))) /
    ((dqm+dqp)*(dqm+dqp)+small_);
  return limit;
}

INSTANTIATE_KERNEL(VolumeOfFluidScsUpwAdvElemKernel);

}  // nalu
}  // sierra
