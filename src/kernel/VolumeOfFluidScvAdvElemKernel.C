/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/VolumeOfFluidScvAdvElemKernel.h"
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
VolumeOfFluidScvAdvElemKernel<AlgTraits>::VolumeOfFluidScvAdvElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ScalarFieldType* vof,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap())
{
  // save off fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  
  vofNp1_ = &(vof->field_of_state(stk::mesh::StateNP1));

  std::string velocity_name = solnOpts.does_mesh_move() ? "velocity_rtm" : "velocity";
  velocityRTM_ = metaData.get_field<double>(stk::topology::NODE_RANK, velocity_name);
  coordinates_ = metaData.get_field<double>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);

  // compute shape function
  get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shape_fcn(ptr);}, v_shape_function_);

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*velocityRTM_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*vofNp1_, 1);
  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(SCV_GRAD_OP, CURRENT_COORDINATES);
}

template<typename AlgTraits>
VolumeOfFluidScvAdvElemKernel<AlgTraits>::~VolumeOfFluidScvAdvElemKernel()
{}

template<typename AlgTraits>
void
VolumeOfFluidScvAdvElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>& lhs,
  SharedMemView<DoubleType *>&rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  NALU_ALIGNED DoubleType w_uIp[AlgTraits::nDim_];

  SharedMemView<DoubleType*>& v_vofNp1 = scratchViews.get_scratch_view_1D(
    *vofNp1_);
  SharedMemView<DoubleType**>& v_vrtm = scratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<DoubleType*>& v_scv_volume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;
  SharedMemView<DoubleType***>& v_dndx = scratchViews.get_me_views(CURRENT_COORDINATES).dndx_scv;

  for ( int ip = 0; ip < AlgTraits::numScvIp_; ++ip ) {

    // nearest node to ip
    const int nearestNode = ipNodeMap_[ip];
    
    // zero out; scalar
    DoubleType vofNp1Ip = 0.0;
    
    // zero out; vector
    for (int j = 0; j < AlgTraits::nDim_; ++j) {
      w_uIp[j] = 0.0;
    }

    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      // save off shape function
      const DoubleType r = v_shape_function_(ip,ic);
      vofNp1Ip += r*v_vofNp1(ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        w_uIp[j] += r*v_vrtm(ic,j);
      }
    }

    const DoubleType scV = v_scv_volume(ip);
    
    // Galerkin first
    DoubleType uidVofdxi = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      DoubleType lhsAdv = 0.0;
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        lhsAdv += w_uIp[j]*v_dndx(ip,ic,j);
      }
      uidVofdxi += lhsAdv*v_vofNp1(ic);
      lhs(nearestNode,ic) += lhsAdv*scV;
    }
    
    // assemble rhs
    rhs(nearestNode) -= uidVofdxi*scV;
  }
}

INSTANTIATE_KERNEL(VolumeOfFluidScvAdvElemKernel);

}  // nalu
}  // sierra
