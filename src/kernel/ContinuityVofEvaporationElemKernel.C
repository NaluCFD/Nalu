/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/ContinuityVofEvaporationElemKernel.h"
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
ContinuityVofEvaporationElemKernel<AlgTraits>::ContinuityVofEvaporationElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    rhoL_(1.0),    // cgs g/cm^3
    rhoG_(1.2e-3), // g/cm^3
    jm_(solnOpts.evapJm_), // cgs g/cm^2/s
    m_(solnOpts.evapM_),
    n_(solnOpts.evapN_),
    c_(solnOpts.evapC_),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap())
{
  // save off fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  
  vofNp1_ = metaData.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "volume_of_fluid");
  coordinates_ = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);

  // compute shape function
  get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shape_fcn(ptr);}, v_shape_function_);

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*vofNp1_, 1);
  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(SCV_GRAD_OP, CURRENT_COORDINATES);
}

template<typename AlgTraits>
ContinuityVofEvaporationElemKernel<AlgTraits>::~ContinuityVofEvaporationElemKernel()
{}

template<typename AlgTraits>
void
ContinuityVofEvaporationElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>& lhs,
  SharedMemView<DoubleType *>&rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  NALU_ALIGNED DoubleType w_gradVof[AlgTraits::nDim_];

  SharedMemView<DoubleType*>& v_vofNp1 = scratchViews.get_scratch_view_1D(
    *vofNp1_);
  SharedMemView<DoubleType*>& v_scv_volume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;
  SharedMemView<DoubleType***>& v_dndx = scratchViews.get_me_views(CURRENT_COORDINATES).dndx_scv;

  for ( int ip = 0; ip < AlgTraits::numScvIp_; ++ip ) {

    // nearest node to ip
    const int nearestNode = ipNodeMap_[ip];
    
    // zero out; vector
    for (int j = 0; j < AlgTraits::nDim_; ++j) {
      w_gradVof[j] = 0.0;
    }

    // compute scv quantities
    DoubleType vofNp1Ip = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      vofNp1Ip += v_shape_function_(ip,ic)*v_vofNp1(ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        w_gradVof[j] += vofNp1Ip*v_dndx(ip,ic,j);
      }
    }

    // magnitude of gradient
    DoubleType magGradVof = 0.0;
    for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
      magGradVof += w_gradVof[j]*w_gradVof[j];
    }
    magGradVof = stk::math::sqrt(magGradVof);

    // compute source term
    const DoubleType sM = c_*stk::math::pow(vofNp1Ip,n_)*stk::math::pow((1.0-vofNp1Ip),m_)*magGradVof*jm_;
 
    // assemble rhs (no RHS)
    rhs(nearestNode) -= c_*sM*(1.0/rhoG_ - 1.0/rhoL_)*v_scv_volume(ip);
  }
}

INSTANTIATE_KERNEL(ContinuityVofEvaporationElemKernel);

}  // nalu
}  // sierra
