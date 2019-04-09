/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/MomentumDiffFemKernel.h"
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
MomentumDiffFemKernel<AlgTraits>::MomentumDiffFemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  VectorFieldType* velocity,
  ScalarFieldType* viscosity,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    viscosity_(viscosity),
    includeDivU_(solnOpts.includeDivU_),
    shiftedGradOp_(solnOpts.get_shifted_grad_op(velocity->name()))
{
  // Save of required fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();

  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
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
  dataPreReqs.add_gathered_nodal_field(*viscosity_, 1);
  if ( shiftedGradOp_ )
    dataPreReqs.add_master_element_call(FEM_SHIFTED_GRAD_OP, CURRENT_COORDINATES);
  else
    dataPreReqs.add_master_element_call(FEM_GRAD_OP, CURRENT_COORDINATES);
}

template<typename AlgTraits>
MomentumDiffFemKernel<AlgTraits>::~MomentumDiffFemKernel()
{
  // does nothing
}

template<typename AlgTraits>
void
MomentumDiffFemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType**>& lhs,
  SharedMemView<DoubleType*>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  SharedMemView<DoubleType**>& v_uNp1 = scratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<DoubleType*>& v_viscosity = scratchViews.get_scratch_view_1D(*viscosity_);
  
  SharedMemView<DoubleType***>& v_dndx = scratchViews.get_me_views(CURRENT_COORDINATES).dndx_fem;
  SharedMemView<DoubleType*>& v_det_j = scratchViews.get_me_views(CURRENT_COORDINATES).det_j_fem;

  // -dw/dxj*Fij*detJ*weight; Fij = -2*mu*S^*ij; S^*ij = Sij - 1/3*Skk*delij); Fij = -mu(dui/dxj+duj/dxi) + 2/3*divU*delij

  for ( int ip = 0; ip < AlgTraits::numGp_; ++ip ) {

    // zero ip values; sneak in divU
    DoubleType muIp = 0.0;
    DoubleType divU = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      muIp += v_viscosity(ic)*v_shape_function_(ip,ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        divU += v_dndx(ip,ic,j)*v_uNp1(ic,j);
      }
    }
    
    const DoubleType divUstress = 2.0/3.0*muIp*divU*includeDivU_;

    // start the assembly
    const DoubleType ipFactor = v_det_j(ip)*v_ip_weight_(ip);

    // row ir
    for ( int ir = 0; ir < AlgTraits::nodesPerElement_; ++ir ) {

      const int irNdim = ir*AlgTraits::nDim_;

      // column ic
      for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
        
        const int icNdim = ic*AlgTraits::nDim_;
        
        // component i
        for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
          
          // divU contribution (explicit)
          rhs(irNdim+i) -= -divUstress*v_dndx(ip,ir,i)*ipFactor;
          
          // viscous stress
          DoubleType lhsSum = 0.0;
          for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
            
            // mu*dw/dxj*dui/dxj*det_j*weight; fixed i over j loop; see below..
            lhsSum += muIp*v_dndx(ip,ir,j)*v_dndx(ip,ic,j)*ipFactor;
            
            // mu*dw/dxj*duj/dxi*det_j*weight
            const DoubleType lhsfac = muIp*v_dndx(ip,ir,j)*v_dndx(ip,ic,i)*ipFactor;
            lhs(irNdim+i,icNdim+j) += lhsfac;
            rhs(irNdim+i) -= lhsfac*v_uNp1(ic,j);
          }
          
          // deal with accumulated lhs and flux for mu*dw/dxj*dui/dxj*det_j*weight
          lhs(irNdim+i,icNdim+i) += lhsSum;
          rhs(irNdim+i) -= lhsSum*v_uNp1(ic,i);
        }
      }
    }
  }
}    

INSTANTIATE_FEM_KERNEL(MomentumDiffFemKernel);

}  // nalu
}  // sierra
