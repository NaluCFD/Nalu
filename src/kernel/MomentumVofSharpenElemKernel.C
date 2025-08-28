/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/MomentumVofSharpenElemKernel.h"
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
MomentumVofSharpenElemKernel<AlgTraits>::MomentumVofSharpenElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  VectorFieldType* velocity,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    lrscv_(sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_)->adjacentNodes()),
    cAlpha_(solnOpts.vofCalpha_),
    densDiff_(solnOpts.vofDensityPhaseOne_ - solnOpts.vofDensityPhaseTwo_)
{
  // save off fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  
  coordinates_ = metaData.get_field<double>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  std::string velocity_name = solnOpts.does_mesh_move() ? "velocity_rtm" : "velocity";
  velocityRTM_ = metaData.get_field<double>(stk::topology::NODE_RANK, velocity_name);
  vof_ = metaData.get_field<double>(stk::topology::NODE_RANK, "volume_of_fluid");
  interfaceNormal_ = metaData.get_field<double>(stk::topology::NODE_RANK, "interface_normal");

  MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_);

  NaluEnv::self().naluOutputP0() << "In MomentumVofSharpenElemKernel; topo is: " << AlgTraits::topo_ << std::endl;

  get_scs_shape_fn_data<AlgTraits>([&](double* ptr){meSCS->shape_fcn(ptr);}, v_shape_function_);

  // add master elements
  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*velocityNp1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityRTM_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*vof_, 1);
  dataPreReqs.add_gathered_nodal_field(*interfaceNormal_, AlgTraits::nDim_);
  dataPreReqs.add_master_element_call(SCS_AREAV, CURRENT_COORDINATES);
}

template<typename AlgTraits>
MomentumVofSharpenElemKernel<AlgTraits>::~MomentumVofSharpenElemKernel()
{}

template<typename AlgTraits>
void
MomentumVofSharpenElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>& lhs,
  SharedMemView<DoubleType *>&rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  DoubleType w_uIp[AlgTraits::nodesPerElement_];
  DoubleType w_magU[AlgTraits::nodesPerElement_];

  SharedMemView<DoubleType**>& v_uNp1 = scratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<DoubleType**>& v_vrtm = scratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<DoubleType*>& v_vof = scratchViews.get_scratch_view_1D(*vof_);
  SharedMemView<DoubleType**>& v_interfaceNormal = scratchViews.get_scratch_view_2D(*interfaceNormal_);
  SharedMemView<DoubleType**>& v_scs_areav = scratchViews.get_me_views(CURRENT_COORDINATES).scs_areav;

  // determine nodal velocity magnitude
  for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
    DoubleType magU = 0.0;
    for ( int k = 0; k < AlgTraits::nDim_; ++k ) {
      magU += v_vrtm(ic,k)*v_vrtm(ic,k);
    }
    w_magU[ic] = stk::math::sqrt(magU);
  }

  for ( int ip = 0; ip < AlgTraits::numScsIp_; ++ip ) {

    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];

    // save off some offsets
    const int ilNdim = il*AlgTraits::nDim_;
    const int irNdim = ir*AlgTraits::nDim_;

    // zero
    for ( int k = 0; k < AlgTraits::nDim_; ++k ) {
      w_uIp[k] = 0.0;
    }
    
    // interpolate to ips
    DoubleType sharpenIp = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      DoubleType magU = 0.0;
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        magU += v_vrtm(ic,j)*v_vrtm(ic,j);
      }
      const DoubleType r = v_shape_function_(ip,ic);

      // save off ic
      const DoubleType vofIc = v_vof(ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        const DoubleType ucj = cAlpha_*w_magU[ic]*v_interfaceNormal(ic,j);
        sharpenIp += r*ucj*vofIc*(1.0-vofIc)*v_scs_areav(ip,j);
        w_uIp[j] += r*v_uNp1(ic,j);
      }
    }

    // correct sharpen
    sharpenIp *= densDiff_;

    // assemble rhs
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      const int indexL = ilNdim + i;
      const int indexR = irNdim + i;
      rhs(indexL) -= w_uIp[i]*sharpenIp;
      rhs(indexR) += w_uIp[i]*sharpenIp;
    }

    // manage LHS, ui*ucj*alpha*(1-alpha)*(rho1- rho2)*njdS
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const DoubleType lhsfacAdv = v_shape_function_(ip,ic)*sharpenIp;
      const int icNdim = ic*AlgTraits::nDim_;
      for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
        const int indexL = ilNdim + i;
        const int indexR = irNdim + i;
        lhs(indexL,icNdim+i) += lhsfacAdv;
        lhs(indexR,icNdim+i) -= lhsfacAdv;        
      }
    }
    
  }
}

INSTANTIATE_KERNEL(MomentumVofSharpenElemKernel);

}  // nalu
}  // sierra
