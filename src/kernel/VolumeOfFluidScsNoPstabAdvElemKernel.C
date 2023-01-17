/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/VolumeOfFluidScsNoPstabAdvElemKernel.h"
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
VolumeOfFluidScsNoPstabAdvElemKernel<AlgTraits>::VolumeOfFluidScsNoPstabAdvElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ScalarFieldType* vof,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    lrscv_(sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_)->adjacentNodes())
{
  // save off fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  
  vofNp1_ = &(vof->field_of_state(stk::mesh::StateNP1));

  coordinates_ = metaData.get_field<double>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  std::string velocity_name = solnOpts.does_mesh_move() ? "velocity_rtm" : "velocity";
  velocityRTM_ = metaData.get_field<double>(stk::topology::NODE_RANK, velocity_name);

  MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_);

  // compute shape function
  get_scs_shape_fn_data<AlgTraits>([&](double* ptr){meSCS->shape_fcn(ptr);}, v_shape_function_);

  // add master elements
  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*vofNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*velocityRTM_, AlgTraits::nDim_);
  dataPreReqs.add_master_element_call(SCS_AREAV, CURRENT_COORDINATES);
}

template<typename AlgTraits>
VolumeOfFluidScsNoPstabAdvElemKernel<AlgTraits>::~VolumeOfFluidScsNoPstabAdvElemKernel()
{}

template<typename AlgTraits>
void
VolumeOfFluidScsNoPstabAdvElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>& lhs,
  SharedMemView<DoubleType *>&rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  NALU_ALIGNED DoubleType w_velocityRTMIp[AlgTraits::nDim_];

  SharedMemView<DoubleType*>& v_vofNp1 = scratchViews.get_scratch_view_1D(*vofNp1_);
  SharedMemView<DoubleType**>& v_velocityRTM = scratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<DoubleType**>& v_scs_areav = scratchViews.get_me_views(CURRENT_COORDINATES).scs_areav;

  //============================================
  // SCS contribution; vdot*alpha*nj*dS
  //============================================
  for ( int ip = 0; ip < AlgTraits::numScsIp_; ++ip ) {
    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];

    // zero out values of interest for this ip
    for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
      w_velocityRTMIp[j] = 0.0;
    }

    // compute ip quantities
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const DoubleType r = v_shape_function_(ip,ic);
      for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
        w_velocityRTMIp[i] += r*v_velocityRTM(ic,i);
      }
    }

    // compute vdot
    DoubleType vdot = 0.0;
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      vdot += w_velocityRTMIp[i]*v_scs_areav(ip,i);
    }

    // advection and diffusion
    DoubleType vofAdv = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {

      // advection
      const DoubleType lhsfacAdv = v_shape_function_(ip,ic)*vdot;
      vofAdv += lhsfacAdv*v_vofNp1(ic);

      // lhs; il then ir
      lhs(il,ic) += lhsfacAdv;
      lhs(ir,ic) -= lhsfacAdv;
    }

    // rhs; il then ir
    rhs(il) -= vofAdv;
    rhs(ir) += vofAdv;
  }
}

INSTANTIATE_KERNEL(VolumeOfFluidScsNoPstabAdvElemKernel);

}  // nalu
}  // sierra
