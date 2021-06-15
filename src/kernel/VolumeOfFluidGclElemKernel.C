/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/VolumeOfFluidGclElemKernel.h"
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
VolumeOfFluidGclElemKernel<AlgTraits>::VolumeOfFluidGclElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ScalarFieldType* vof,
  ElemDataRequests& dataPreReqs,
  const bool lumpedMass)
  : Kernel(),
    lumpedMass_(lumpedMass),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap())
{
  // save off fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();

  vofNp1_ = &(vof->field_of_state(stk::mesh::StateNP1));
  divV_ = metaData.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "div_mesh_velocity");
  coordinates_ = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);

  // compute shape function
  if ( lumpedMass_ )
    get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shifted_shape_fcn(ptr);}, v_shape_function_);
  else
    get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shape_fcn(ptr);}, v_shape_function_);

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*divV_, 1);
  dataPreReqs.add_gathered_nodal_field(*vofNp1_, 1);
  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
}

template<typename AlgTraits>
VolumeOfFluidGclElemKernel<AlgTraits>::~VolumeOfFluidGclElemKernel()
{}

template<typename AlgTraits>
void
VolumeOfFluidGclElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>&lhs,
  SharedMemView<DoubleType *>&rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  SharedMemView<DoubleType*>& v_vofNp1 = scratchViews.get_scratch_view_1D(
    *vofNp1_);
  SharedMemView<DoubleType*>& v_divV = scratchViews.get_scratch_view_1D(
    *divV_);
  SharedMemView<DoubleType*>& v_scv_volume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;

  for (int ip=0; ip < AlgTraits::numScvIp_; ++ip) {
    const int nearestNode = ipNodeMap_[ip];

    DoubleType divVIp = 0.0;
    DoubleType vofIp = 0.0;
    for (int ic=0; ic < AlgTraits::nodesPerElement_; ++ic) {
      const DoubleType r = v_shape_function_(ip, ic);
      divVIp += r*v_divV(ic);
      vofIp += r*v_vofNp1(ic);
    }

    const DoubleType scV = v_scv_volume(ip);

    rhs(nearestNode) -= vofIp*divVIp*scV;

    // Compute LHS
    const DoubleType lhsFac = divVIp*scV;
    for (int ic=0; ic < AlgTraits::nodesPerElement_; ++ic) {
      lhs(nearestNode,ic) += v_shape_function_(ip, ic)*lhsFac;
    }
  }
}

INSTANTIATE_KERNEL(VolumeOfFluidGclElemKernel);

}  // nalu
}  // sierra
