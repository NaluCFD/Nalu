/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "nso/ScalarNSOFemKernel.h"
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
ScalarNSOFemKernel<AlgTraits>::ScalarNSOFemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ScalarFieldType* scalarQ,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    shiftedGradOp_(solnOpts.get_shifted_grad_op(scalarQ->name()))
{
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  ScalarFieldType *density = metaData.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "density");

  scalarQNp1_ = &(scalarQ->field_of_state(stk::mesh::StateNP1));
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));

  if (solnOpts.does_mesh_move())
    velocityRTM_ = metaData.get_field<VectorFieldType>(
      stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = metaData.get_field<VectorFieldType>(
      stk::topology::NODE_RANK, "velocity");

  // extract field not required in execute()
  VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  // extract master element
  MasterElement *meFEM = sierra::nalu::MasterElementRepo::get_fem_master_element(AlgTraits::topo_);
  
  // copy ip weights into our 1-d view
  for ( int k = 0; k < AlgTraits::numGp_; ++k )
    v_ip_weight_[k] = meFEM->weights_[k];
  
  // master element, shape function
  get_fem_shape_fn_data<AlgTraits>([&](double* ptr){meFEM->shape_fcn(ptr);}, v_shape_function_);

  dataPreReqs.add_fem_volume_me(meFEM);

  // fields
  dataPreReqs.add_gathered_nodal_field(*scalarQNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_,1);
  dataPreReqs.add_gathered_nodal_field(*velocityRTM_, AlgTraits::nDim_);
  dataPreReqs.add_coordinates_field(*coordinates, AlgTraits::nDim_, CURRENT_COORDINATES);
  
  // master element data
  if ( shiftedGradOp_ )
    dataPreReqs.add_master_element_call(FEM_SHIFTED_GRAD_OP, CURRENT_COORDINATES);
  else
    dataPreReqs.add_master_element_call(FEM_GRAD_OP, CURRENT_COORDINATES);

  dataPreReqs.add_master_element_call(FEM_GIJ, CURRENT_COORDINATES);
}

template<typename AlgTraits>
void
ScalarNSOFemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType**>& lhs,
  SharedMemView<DoubleType *>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  NALU_ALIGNED DoubleType w_dqdxIp    [AlgTraits::nDim_];
  NALU_ALIGNED DoubleType w_rhoVrtmIp [AlgTraits::nDim_];

  SharedMemView<DoubleType*>& v_qNp1 = scratchViews.get_scratch_view_1D(*scalarQNp1_);
  SharedMemView<DoubleType**>& v_velocityRTM = scratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<DoubleType*>& v_rhoNp1 = scratchViews.get_scratch_view_1D(*densityNp1_);

  SharedMemView<DoubleType***>& v_dndx = scratchViews.get_me_views(CURRENT_COORDINATES).dndx_fem;
  SharedMemView<DoubleType*>& v_det_j = scratchViews.get_me_views(CURRENT_COORDINATES).det_j_fem;

  SharedMemView<DoubleType***>& v_gijUpper = scratchViews.get_me_views(CURRENT_COORDINATES).gijUpper;
  SharedMemView<DoubleType***>& v_gijLower = scratchViews.get_me_views(CURRENT_COORDINATES).gijLower;

  for ( int ip = 0; ip < AlgTraits::numGp_; ++ip ) {

    // zero out; scalar
    DoubleType qNp1Ip = 0.0;
    DoubleType dFdxAdv = 0.0;
    DoubleType dFdxCont = 0.0;

    // zero out; vector
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      w_dqdxIp[i] = 0.0;
      w_rhoVrtmIp[i] = 0.0;
    }

    // determine scs values of interest
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      // save off shape function
      const DoubleType r = v_shape_function_(ip,ic);

      qNp1Ip += r*v_qNp1(ic);

      // compute ip derivatives and flux derivative
      const DoubleType qIC = v_qNp1(ic);
      const DoubleType rhoIC = v_rhoNp1(ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        const DoubleType dnj = v_dndx(ip,ic,j);
        const DoubleType vrtmj = v_velocityRTM(ic,j);
        w_dqdxIp[j] += qIC*dnj;
        w_rhoVrtmIp[j] += r*rhoIC*vrtmj;
        dFdxAdv += rhoIC*vrtmj*qIC*dnj;
        dFdxCont += rhoIC*vrtmj*dnj;
      }
    }

    // compute residual for NSO; linearized first
    DoubleType residual = dFdxAdv - qNp1Ip*dFdxCont;
    for ( int j = 0; j < AlgTraits::nDim_; ++j )
      residual -= w_rhoVrtmIp[j]*w_dqdxIp[j];

    // denominator for nu as well as terms for "upwind" nu
    DoubleType gUpperMagGradQ = 0.0;
    DoubleType rhoVrtmiGLowerRhoVrtmj = 0.0;
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      const DoubleType dqdxIpi = w_dqdxIp[i];
      const DoubleType rhoVrtmi = w_rhoVrtmIp[i];
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        gUpperMagGradQ += dqdxIpi*v_gijUpper(ip,i,j)*w_dqdxIp[j];
        rhoVrtmiGLowerRhoVrtmj += rhoVrtmi*v_gijLower(ip,i,j)*w_rhoVrtmIp[j];
      }
    }

    // construct nu from residual
    const DoubleType nuResidual = stk::math::sqrt((residual*residual)/(gUpperMagGradQ+small_));

    // construct nu from first-order-like approach; SNL-internal write-up (eq 209)
    // for now, only include advection as full set of terms is too diffuse
    const DoubleType nuFirstOrder = stk::math::sqrt(rhoVrtmiGLowerRhoVrtmj);

    // limit based on first order; Cupw_ is a fudge factor similar to Guermond's approach
    const DoubleType nu = stk::math::min(Cupw_*nuFirstOrder, nuResidual);

    // start the assembly
    const DoubleType ipFactor = v_det_j(ip)*v_ip_weight_(ip);
    
    // row ir
    for ( int ir = 0; ir < AlgTraits::nodesPerElement_; ++ir) {
      
      // column ic
      DoubleType gijFac = 0.0;
      for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
        
        // save of some variables
        const DoubleType qIC = v_qNp1(ic);
        
        // NSO diffusion-like term; + nu * dw(ir)/dxi * g^ij * dq(ic)/dxj * weight(ip) * detJ(ip)
        DoubleType lhsSum = 0.0;
        for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
          const DoubleType dndxi = v_dndx(ip,ir,i);
          for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
            const DoubleType fac = dndxi* v_gijUpper(ip,i,j)*v_dndx(ip,ic,j);
            gijFac += fac*qIC;
            lhsSum += fac;
          }
        }
        
        lhs(ir,ic) += nu*lhsSum*ipFactor;
      }
      
      // residual
      rhs(ir) -= nu*gijFac*ipFactor;
    }
  }
}

INSTANTIATE_FEM_KERNEL(ScalarNSOFemKernel);

}  // nalu
}  // sierra
