/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/MomentumMassFemKernel.h"
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
MomentumMassFemKernel<AlgTraits>::MomentumMassFemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  VectorFieldType* velocity,
  ScalarFieldType* density,
  ElemDataRequests& dataPreReqs)
  : Kernel()
{
  // Save of required fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  coordinates_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  Gjp_ =  metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx");
  
  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  velocityN_ = &(velocity->field_of_state(stk::mesh::StateN));
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  densityN_ = &(density->field_of_state(stk::mesh::StateN));
  if (velocity->number_of_states() == 2) {
    velocityNm1_ = velocityN_;
    densityNm1_ = densityN_;
  }
  else {
    velocityNm1_ = &(velocity->field_of_state(stk::mesh::StateNM1));
    densityNm1_ = &(density->field_of_state(stk::mesh::StateNM1));
  }
  
  // extract master element
  MasterElement *meFEM = sierra::nalu::MasterElementRepo::get_fem_master_element(AlgTraits::topo_);
  
  // copy ip weights into our 1-d view
  for ( int k = 0; k < AlgTraits::numGp_; ++k )
    v_ip_weight_[k] = meFEM->weights_[k];
  
  // master element, shape function is shifted consistently
  if ( solnOpts.get_shifted_grad_op(velocity->name()) )
    get_fem_shape_fn_data<AlgTraits>([&](double* ptr){meFEM->shifted_shape_fcn(ptr);}, v_shape_function_);
  else
    get_fem_shape_fn_data<AlgTraits>([&](double* ptr){meFEM->shape_fcn(ptr);}, v_shape_function_);
  
  // add FEM master element
  dataPreReqs.add_fem_volume_me(meFEM);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*velocityNp1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityN_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityNm1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*Gjp_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityN_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityNm1_, 1);
  
  dataPreReqs.add_master_element_call(FEM_DET_J, CURRENT_COORDINATES);
}

template<typename AlgTraits>
MomentumMassFemKernel<AlgTraits>::~MomentumMassFemKernel()
{
  // does nothing
}

template<typename AlgTraits>
void
MomentumMassFemKernel<AlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
  dt_ = timeIntegrator.get_time_step();
  gamma1_ = timeIntegrator.get_gamma1();
  gamma2_ = timeIntegrator.get_gamma2();
  gamma3_ = timeIntegrator.get_gamma3(); // gamma3 may be zero
}

template<typename AlgTraits>
void
MomentumMassFemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType**>& lhs,
  SharedMemView<DoubleType*>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  SharedMemView<DoubleType**>& v_uNp1 = scratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<DoubleType**>& v_uN = scratchViews.get_scratch_view_2D(*velocityN_);
  SharedMemView<DoubleType**>& v_uNm1 = scratchViews.get_scratch_view_2D(*velocityNm1_);
  SharedMemView<DoubleType**>& v_Gp = scratchViews.get_scratch_view_2D(*Gjp_);
  SharedMemView<DoubleType*>& v_rhoNp1 = scratchViews.get_scratch_view_1D(*densityNp1_);
  SharedMemView<DoubleType*>& v_rhoN = scratchViews.get_scratch_view_1D(*densityN_);
  SharedMemView<DoubleType*>& v_rhoNm1 = scratchViews.get_scratch_view_1D(*densityNm1_);
  
  SharedMemView<DoubleType*>& v_det_j = scratchViews.get_me_views(CURRENT_COORDINATES).det_j_fem;

  for ( int ip = 0; ip < AlgTraits::numGp_; ++ip ) {
  
    // zero out ip values
    DoubleType rhoNp1Ip = 0.0;
    DoubleType rhoNIp = 0.0;
    DoubleType rhoNm1Ip = 0.0;
    // evaluate ip values
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const DoubleType r = v_shape_function_(ip,ic);
      rhoNp1Ip += r*v_rhoNp1(ic);
      rhoNIp += r*v_rhoN(ic);
      rhoNm1Ip += r*v_rhoNm1(ic);
    }

    // start the assembly (collect ip scalings)
    const DoubleType ipFactor = v_det_j(ip)*v_ip_weight_(ip);

    // row ir
    for ( int ir = 0; ir < AlgTraits::nodesPerElement_; ++ir) {
      
      const int irNdim = ir*AlgTraits::nDim_;
      const DoubleType wIr = v_shape_function_(ip,ir);
      
      // column ic
      for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
        
        const int icNdim = ic*AlgTraits::nDim_;
        
        // component i (avoid w_uN*Ip and w_GjpIp fields)
        for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
          lhs(irNdim+i,icNdim+i) += wIr*gamma1_*rhoNp1Ip*v_shape_function_(ip,ic)/dt_*ipFactor;
          rhs(irNdim+i) -= wIr*((gamma1_*rhoNp1Ip*v_uNp1(ic,i) 
                                 + gamma2_*rhoNIp*v_uN(ic,i)
                                 + gamma3_*rhoNm1Ip*v_uNm1(ic,i))/dt_ 
                                + v_Gp(ic,i))*v_shape_function_(ip,ic)*ipFactor; 
        }
      }
    }
  }
}
  
INSTANTIATE_FEM_KERNEL(MomentumMassFemKernel);

}  // nalu
}  // sierra
