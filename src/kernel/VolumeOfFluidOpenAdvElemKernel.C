/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/VolumeOfFluidOpenAdvElemKernel.h"
#include "EquationSystem.h"
#include "master_element/MasterElement.h"
#include "PecletFunction.h"
#include "SolutionOptions.h"
#include "BuildTemplates.h"

// template and scratch space
#include "ScratchViews.h"

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra {
namespace nalu {

template<typename BcAlgTraits>
VolumeOfFluidOpenAdvElemKernel<BcAlgTraits>::VolumeOfFluidOpenAdvElemKernel(
  const stk::mesh::MetaData &metaData,
  const SolutionOptions &solnOpts,
  EquationSystem *eqSystem,
  ScalarFieldType *scalarQ,
  ScalarFieldType *bcScalarQ,
  ElemDataRequests &faceDataPreReqs,
  ElemDataRequests &elemDataPreReqs)
  : Kernel(),
    scalarQ_(scalarQ),
    bcScalarQ_(bcScalarQ),
    meSCS_(sierra::nalu::MasterElementRepo::get_surface_master_element(BcAlgTraits::elemTopo_))
{
  // save off fields
  coordinates_ = metaData.get_field<double>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  openVolumeFlowRate_ = metaData.get_field<double>(metaData.side_rank(), "open_volume_flow_rate");
  
  // extract master elements
  MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(BcAlgTraits::faceTopo_);
  
  // add master elements
  faceDataPreReqs.add_cvfem_face_me(meFC);
 
  // fields and data; face and then element
  faceDataPreReqs.add_coordinates_field(*coordinates_, BcAlgTraits::nDim_, CURRENT_COORDINATES);
  faceDataPreReqs.add_gathered_nodal_field(*scalarQ_, 1);
  faceDataPreReqs.add_gathered_nodal_field(*bcScalarQ_, 1);
  faceDataPreReqs.add_face_field(*openVolumeFlowRate_, BcAlgTraits::numFaceIp_);

  // never shift evaluation
  get_face_shape_fn_data<BcAlgTraits>([&](double* ptr){meFC->shape_fcn(ptr);}, vf_adv_shape_function_);
}

template<typename BcAlgTraits>
VolumeOfFluidOpenAdvElemKernel<BcAlgTraits>::~VolumeOfFluidOpenAdvElemKernel()
{}

template<typename BcAlgTraits>
void
VolumeOfFluidOpenAdvElemKernel<BcAlgTraits>::execute(
  SharedMemView<DoubleType**> &lhs,
  SharedMemView<DoubleType *> &rhs,
  ScratchViews<DoubleType> &faceScratchViews,
  ScratchViews<DoubleType> &elemScratchViews,
  int elemFaceOrdinal)
{  
  const int *face_node_ordinals = meSCS_->side_node_ordinals(elemFaceOrdinal);
 
  // face
  SharedMemView<DoubleType*>& vf_scalarQ = faceScratchViews.get_scratch_view_1D(*scalarQ_);
  SharedMemView<DoubleType*>& vf_bcScalarQ = faceScratchViews.get_scratch_view_1D(*bcScalarQ_);
  SharedMemView<DoubleType*>& vf_openVolumeFlowRate = faceScratchViews.get_scratch_view_1D(*openVolumeFlowRate_);
  
  for (int ip=0; ip < BcAlgTraits::numFaceIp_; ++ip) {
    
    const int nearestNode = meSCS_->ipNodeMap(elemFaceOrdinal)[ip]; // "Right"
        
    // interpolate to bip
    DoubleType qIp = 0.0;
    DoubleType qIpEntrain = 0.0;
    for ( int ic = 0; ic < BcAlgTraits::nodesPerFace_; ++ic ) {
      const DoubleType rAdv = vf_adv_shape_function_(ip,ic);
      qIp += rAdv*vf_scalarQ(ic);
      qIpEntrain += rAdv*vf_bcScalarQ(ic);
    }
    
    const DoubleType tvdot = vf_openVolumeFlowRate(ip);

    // account for both cases, i.e., leaving or entering to avoid hard loop over SIMD length
    const DoubleType vdotLeaving = stk::math::if_then_else(tvdot > 0, 1.0, 0.0);
    const DoubleType om_vdotLeaving = 1.0 - vdotLeaving;

    //================================
    // flow is leaving
    //================================
    
    // central is simply qIp*vdot
    rhs(nearestNode) -= tvdot*qIp*vdotLeaving;

    // central part
    const DoubleType fac = tvdot*vdotLeaving;
    for ( int ic = 0; ic < BcAlgTraits::nodesPerFace_; ++ic ) {
      const int nn = face_node_ordinals[ic];
      lhs(nearestNode,nn) += vf_adv_shape_function_(ip,ic)*fac;
    }

    //================================
    // flow is entering
    //================================

    // extrainment; advect in from specified value
    rhs(nearestNode) -= tvdot*qIpEntrain*om_vdotLeaving;    
  }
}

INSTANTIATE_KERNEL_FACE_ELEMENT(VolumeOfFluidOpenAdvElemKernel);

}  // nalu
}  // sierra
