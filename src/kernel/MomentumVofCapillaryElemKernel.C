/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/MomentumVofCapillaryElemKernel.h"
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

template<class AlgTraits>
MomentumVofCapillaryElemKernel<AlgTraits>::MomentumVofCapillaryElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  VectorFieldType* velocity,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    lrscv_(sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_)->adjacentNodes()),
    shiftedGradOp_(solnOpts.get_shifted_grad_op(velocity->name()))
{
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  coordinates_ = metaData.get_field<double>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  interfaceNormal_ = metaData.get_field<double>(stk::topology::NODE_RANK, "interface_normal");
  vof_ = metaData.get_field<double>(stk::topology::NODE_RANK, "volume_of_fluid");
  surfaceTension_ = metaData.get_field<double>(stk::topology::NODE_RANK, "surface_tension");

  MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_);

  get_scs_shape_fn_data<AlgTraits>([&](double* ptr){meSCS->shape_fcn(ptr);}, v_shape_function_);

  // add master elements
  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields and data; mdot not gathered as element data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*velocityNp1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*interfaceNormal_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*vof_, 1);
  dataPreReqs.add_gathered_nodal_field(*surfaceTension_, 1);
  dataPreReqs.add_master_element_call(SCS_AREAV, CURRENT_COORDINATES);
  if ( shiftedGradOp_ )
    dataPreReqs.add_master_element_call(SCS_SHIFTED_GRAD_OP, CURRENT_COORDINATES);
  else
    dataPreReqs.add_master_element_call(SCS_GRAD_OP, CURRENT_COORDINATES);
}

template<class AlgTraits>
MomentumVofCapillaryElemKernel<AlgTraits>::~MomentumVofCapillaryElemKernel()
{}

template<typename AlgTraits>
void
MomentumVofCapillaryElemKernel<AlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
  dt_ = timeIntegrator.get_time_step();
}

template<class AlgTraits>
void
MomentumVofCapillaryElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>& lhs,
  SharedMemView<DoubleType *>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  DoubleType w_nxIp[AlgTraits::nDim_];
  DoubleType w_dvofdxIp[AlgTraits::nDim_];

  SharedMemView<DoubleType**>& v_uNp1 = scratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<DoubleType**>& v_interfaceNormal = scratchViews.get_scratch_view_2D(*interfaceNormal_);
  SharedMemView<DoubleType*>& v_vof = scratchViews.get_scratch_view_1D(*vof_);
  SharedMemView<DoubleType*>& v_surfaceTension = scratchViews.get_scratch_view_1D(*surfaceTension_);

  SharedMemView<DoubleType**>& v_scs_areav = scratchViews.get_me_views(CURRENT_COORDINATES).scs_areav;
  SharedMemView<DoubleType***>& v_dndx = shiftedGradOp_
    ? scratchViews.get_me_views(CURRENT_COORDINATES).dndx_shifted 
    : scratchViews.get_me_views(CURRENT_COORDINATES).dndx;

  for ( int ip = 0; ip < AlgTraits::numScsIp_; ++ip ) {

    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];

    // save off some offsets
    const int ilNdim = il*AlgTraits::nDim_;
    const int irNdim = ir*AlgTraits::nDim_;

    // compute scs point values; zero
    DoubleType sigmaIp = 0.0;
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      w_nxIp[i] = 0.0;
      w_dvofdxIp[i] = 0.0;
    }

    // now, accumulate
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const DoubleType r = v_shape_function_(ip,ic);
      sigmaIp += r*v_surfaceTension(ic);
      DoubleType vofIc = v_vof(ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        w_nxIp[j] += r*v_interfaceNormal(ic,j);
        w_dvofdxIp[j] += vofIc*v_dndx(ip,ic,j);
      }
    }
    
    // form diffusional flux coefficient
    DoubleType vofMagIp = 0.0;
    for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
      vofMagIp += w_dvofdxIp[j]*w_dvofdxIp[j];
    }
    const DoubleType muIp = sigmaIp*dt_*stk::math::sqrt(vofMagIp); 
      
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      
      const int icNdim = ic*AlgTraits::nDim_;
      
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        
        const DoubleType axj = v_scs_areav(ip,j);
        const DoubleType dndxj = v_dndx(ip,ic,j);
        
        for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
          
          // matrix entries
          const int indexL = ilNdim + i;
          const int indexR = irNdim + i;
          
          const DoubleType uxi = v_uNp1(ic,i);
          const DoubleType nxi = w_nxIp[i];
          const DoubleType om_nxinxi = 1.0-nxi*nxi;
          
          // -mu*dui/dxj*Aj(1.0-nini)
          DoubleType lhsfac = -muIp*dndxj*axj*om_nxinxi;
          // lhs; il then ir
          lhs(indexL,icNdim+i) += lhsfac;
          lhs(indexR,icNdim+i) -= lhsfac;
          // rhs; il then ir
          rhs(indexL) -= lhsfac*uxi;
          rhs(indexR) += lhsfac*uxi;

          // now we need the -nx*ny*Fy - nx*nz*Fz part
          for ( int l = 0; l < AlgTraits::nDim_; ++l ) {
            
            if ( i != l ) {
              const DoubleType nxinxl = nxi*w_nxIp[l];
              const DoubleType uxl = v_uNp1(ic,l);
              
              // +ni*nl*mu*dul/dxj*Aj
              lhsfac = muIp*dndxj*axj*nxinxl;
              lhs(indexL,icNdim+l) += lhsfac;
              lhs(indexR,icNdim+l) -= lhsfac;
              rhs(indexL) -= lhsfac*uxl;              
              rhs(indexR) += lhsfac*uxl;              
            }
          }
        }
      }
    }
  }
}

INSTANTIATE_KERNEL(MomentumVofCapillaryElemKernel);

}  // nalu
}  // sierra
