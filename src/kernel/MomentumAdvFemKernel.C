/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/MomentumAdvFemKernel.h"
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

template<typename AlgTraits>
MomentumAdvFemKernel<AlgTraits>::MomentumAdvFemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  VectorFieldType* velocity,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    shiftedGradOp_(solnOpts.get_shifted_grad_op(velocity->name())),
    includePstab_(1.0),
    skewFac_(1.0),
    om_skewFac_(1.0-skewFac_),
    projTimeScale_(1.0)
{
  // Save of required fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();

  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  density_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  coordinates_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  
  // extract master element
  MasterElement *meFEM = sierra::nalu::MasterElementRepo::get_fem_master_element(AlgTraits::topo_);
  
  // copy ip weights into our 1-d view
  for ( int k = 0; k < AlgTraits::numGp_; ++k )
    v_ip_weight_[k] = meFEM->weights_[k];
  
  // master element, shape function is shifted consistently
  if ( shiftedGradOp_ )
    get_fem_shape_fn_data<AlgTraits>([&](double* ptr){meFEM->shifted_shape_fcn(ptr);}, v_shape_function_);
  else
    get_fem_shape_fn_data<AlgTraits>([&](double* ptr){meFEM->shape_fcn(ptr);}, v_shape_function_);

  dataPreReqs.add_fem_volume_me(meFEM);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*velocityNp1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*density_, 1);
  if ( shiftedGradOp_ )
    dataPreReqs.add_master_element_call(FEM_SHIFTED_GRAD_OP, CURRENT_COORDINATES);
  else
    dataPreReqs.add_master_element_call(FEM_GRAD_OP, CURRENT_COORDINATES);

  // with pstab
  pressure_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure");
  vrtmL_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "vrtm_lagged");
  GjpL_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx_lagged");
  dataPreReqs.add_gathered_nodal_field(*vrtmL_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*pressure_, 1);
  dataPreReqs.add_gathered_nodal_field(*GjpL_, AlgTraits::nDim_);
}
  
template<typename AlgTraits>
MomentumAdvFemKernel<AlgTraits>::~MomentumAdvFemKernel()
{
  // does nothing
}

template<typename AlgTraits>
void
MomentumAdvFemKernel<AlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
  const double dt = timeIntegrator.get_time_step();
  const double gamma1 = timeIntegrator.get_gamma1();
  projTimeScale_ = dt/gamma1;
}

template<typename AlgTraits>
void
MomentumAdvFemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType**>& lhs,
  SharedMemView<DoubleType*>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  NALU_ALIGNED DoubleType w_rhoUIp[AlgTraits::nDim_];
  NALU_ALIGNED DoubleType w_dpdxIp[AlgTraits::nDim_];
  NALU_ALIGNED DoubleType w_GjpIp[AlgTraits::nDim_];
  
  SharedMemView<DoubleType**>& v_uNp1 = scratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<DoubleType*>& v_rhoNp1 = scratchViews.get_scratch_view_1D(*density_); 
  SharedMemView<DoubleType***>& v_dndx = scratchViews.get_me_views(CURRENT_COORDINATES).dndx_fem;
  SharedMemView<DoubleType*>& v_det_j = scratchViews.get_me_views(CURRENT_COORDINATES).det_j_fem;

  // pstab
  SharedMemView<DoubleType**>& v_vrtm = scratchViews.get_scratch_view_2D(*vrtmL_);
  SharedMemView<DoubleType*>& v_pressure = scratchViews.get_scratch_view_1D(*pressure_);
  SharedMemView<DoubleType**>& v_Gjp = scratchViews.get_scratch_view_2D(*GjpL_);

  for ( int ip = 0; ip < AlgTraits::numGp_; ++ip ) {
    
    // zero ip values
    for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
      w_rhoUIp[j] = 0.0;
      w_GjpIp[j] = 0.0;
      w_dpdxIp[j] = 0.0;
    }
    
    // compute ip variables
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const DoubleType r = v_shape_function_(ip,ic);
      const DoubleType rhoIc = v_rhoNp1(ic);
      const DoubleType pIc = v_pressure(ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j) {
        w_rhoUIp[j] += r*rhoIc*v_vrtm(ic,j);
        w_GjpIp[j] += r*v_Gjp(ic,j);
        w_dpdxIp[j] += v_dndx(ip,ic,j)*pIc;
      }
    }
    
    // start the assembly
    const DoubleType ipFactor = v_det_j(ip)*v_ip_weight_(ip);
    
    // row ir
    for ( int ir = 0; ir < AlgTraits::nodesPerElement_; ++ir ) {
      
      // save off 
      const int irNdim = ir*AlgTraits::nDim_;
      const DoubleType wIr = v_shape_function_(ip,ir);
      
      // compute ip mass flux - possibly including pressure stabilization
      DoubleType massFluxIp = 0.0;
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        const DoubleType pStabj = projTimeScale_*(w_dpdxIp[j] - w_GjpIp[j]);
        massFluxIp -= (w_rhoUIp[j]-includePstab_*pStabj)*v_dndx(ip,ir,j);
      }

      // column ic
      for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
        
        // save off
        const int icNdim = ic*AlgTraits::nDim_;
        const DoubleType wIc = v_shape_function_(ip,ic);
        
        // component i
        for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
          
          // term 1: IBP: -1/2*dw/dxj*[rho*uj]*ui*ipFactor
          const DoubleType lhsfacT1 = skewFac_*massFluxIp*wIc*ipFactor;
          lhs(irNdim+i,icNdim+i) += lhsfacT1;
          rhs(irNdim+i) -= lhsfacT1*v_uNp1(ic,i);
          
          DoubleType lhsSumT2 = 0.0;
          for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
            // term 2: non-IBP: 1/2*w*[rho*uj]*dui/dxj*det_j*weight
            lhsSumT2 += om_skewFac_*wIr*w_rhoUIp[j]*v_dndx(ip,ic,j)*ipFactor;
            // term 3: non-IBP: 1/2*w*[rho*uj]*duj/dxi*det_j*weight
            const DoubleType lhsfacT3 = om_skewFac_*wIr*w_rhoUIp[j]*v_dndx(ip,ic,i)*ipFactor;
            lhs(irNdim+i,icNdim+j) += lhsfacT3*0.0; // turn off duj/dxi sensitivity
            rhs(irNdim+i) -= lhsfacT3*v_uNp1(ic,j);
          }
          
          // deal with accumulated lhs and flux for term 2
          lhs(irNdim+i,icNdim+i) += lhsSumT2;
          rhs(irNdim+i) -= lhsSumT2*v_uNp1(ic,i);
        }
      } 
    }
  }
}

INSTANTIATE_FEM_KERNEL(MomentumAdvFemKernel);

}  // nalu
}  // sierra
