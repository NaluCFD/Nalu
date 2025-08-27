/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/MomentumGclElemKernel.h"
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
MomentumGclElemKernel<AlgTraits>::MomentumGclElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ElemDataRequests& dataPreReqs,
  const bool lumpedMass)
  : Kernel(),
    lumpedMass_(lumpedMass),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap())
{
  // save off fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();

  VectorFieldType *velocity = metaData.get_field<double>(stk::topology::NODE_RANK, "velocity");
  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));

  ScalarFieldType *density = metaData.get_field<double>(stk::topology::NODE_RANK, "density");
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));

  divV_ = metaData.get_field<double>(stk::topology::NODE_RANK, "div_mesh_velocity");
  velocityNp1_ = metaData.get_field<double>(stk::topology::NODE_RANK, "velocity");
  coordinates_ = metaData.get_field<double>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);

  // compute shape function
  if ( lumpedMass_ )
    get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shifted_shape_fcn(ptr);}, v_shape_function_);
  else
    get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shape_fcn(ptr);}, v_shape_function_);

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_gathered_nodal_field(*densityNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*divV_, 1);
  dataPreReqs.add_gathered_nodal_field(*velocityNp1_, AlgTraits::nDim_);
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
}

template<typename AlgTraits>
MomentumGclElemKernel<AlgTraits>::~MomentumGclElemKernel()
{}

template<typename AlgTraits>
void
MomentumGclElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>&lhs,
  SharedMemView<DoubleType *>&rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  DoubleType w_uIp [AlgTraits::nDim_];

  SharedMemView<DoubleType**>& v_velocityNp1 = scratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<DoubleType*>& v_densityNp1 = scratchViews.get_scratch_view_1D(
    *densityNp1_);
  SharedMemView<DoubleType*>& v_divV = scratchViews.get_scratch_view_1D(
    *divV_);
  SharedMemView<DoubleType*>& v_scv_volume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;

  for (int ip=0; ip < AlgTraits::numScvIp_; ++ip) {
    const int nearestNode = ipNodeMap_[ip];

    DoubleType rhoIp = 0.0;
    DoubleType divVIp = 0.0;
    for (int j=0; j < AlgTraits::nDim_; j++) {
      w_uIp[j] = 0.0;
    }

    for (int ic=0; ic < AlgTraits::nodesPerElement_; ++ic) {
      const DoubleType r = v_shape_function_(ip, ic);
      rhoIp += r*v_densityNp1(ic);
      divVIp += r*v_divV(ic);
      for (int j=0; j < AlgTraits::nDim_; j++) {
        w_uIp[j] += r*v_velocityNp1(ic, j);
      }
    }
    
    const DoubleType scV = v_scv_volume(ip);
    const int nnNdim = nearestNode*AlgTraits::nDim_;
    for (int i=0; i < AlgTraits::nDim_; ++i) {
      rhs(nnNdim+i) -= rhoIp*w_uIp[i]*divVIp*scV;
    }

    // Compute LHS
    for (int ic=0; ic < AlgTraits::nodesPerElement_; ++ic) {
      const int icNdim = ic*AlgTraits::nDim_;
      const DoubleType r = v_shape_function_(ip, ic);
      const DoubleType lhsfac = r*rhoIp*divVIp*scV;

      for (int j=0; j<AlgTraits::nDim_; ++j) {
        const int indexNN = nnNdim + j;
        lhs(indexNN,icNdim+j) += lhsfac;
      }
    }
  }
}

INSTANTIATE_KERNEL(MomentumGclElemKernel);

}  // nalu
}  // sierra
