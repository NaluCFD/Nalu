/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/ScalarFluxPenaltyFemKernel.h"
#include "master_element/MasterElement.h"
#include "SolutionOptions.h"
#include "TimeIntegrator.h"

// template and scratch space
#include "BuildTemplates.h"
#include "ScratchViews.h"

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra {
namespace nalu {

template<typename BcAlgTraits>
ScalarFluxPenaltyFemKernel<BcAlgTraits>::ScalarFluxPenaltyFemKernel(
  const stk::mesh::MetaData &metaData,
  const SolutionOptions &solnOpts,
  ScalarFieldType *scalarQ,
  ScalarFieldType *bcScalarQ,
  ScalarFieldType *diffFluxCoeff,
  ElemDataRequests &faceDataPreReqs,
  ElemDataRequests &elemDataPreReqs)
  : Kernel(),
    scalarQ_(scalarQ),
    bcScalarQ_(bcScalarQ),
    diffFluxCoeff_(diffFluxCoeff),
    penaltyFac_(2.0),
    shiftedGradOp_(solnOpts.get_shifted_grad_op(scalarQ->name())),
    meFEM_(sierra::nalu::MasterElementRepo::get_fem_master_element(BcAlgTraits::elemTopo_))
{
  coordinates_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  
  // extract master elements
  MasterElement* meFC = sierra::nalu::MasterElementRepo::get_fem_master_element(BcAlgTraits::faceTopo_);
  
  // copy ip weights into our 1-d view
  for ( int k = 0; k < BcAlgTraits::numFaceIp_; ++k )
    vf_ip_weight_[k] = meFC->weights_[k];
  
  // add master elements
  faceDataPreReqs.add_fem_face_me(meFC);
  elemDataPreReqs.add_fem_volume_me(meFEM_);

  // fields and data; face and then element
  faceDataPreReqs.add_gathered_nodal_field(*scalarQ_, 1);
  faceDataPreReqs.add_gathered_nodal_field(*bcScalarQ_, 1);
  faceDataPreReqs.add_gathered_nodal_field(*diffFluxCoeff_, 1);
  faceDataPreReqs.add_coordinates_field(*coordinates_, BcAlgTraits::nDim_, CURRENT_COORDINATES);

  elemDataPreReqs.add_coordinates_field(*coordinates_, BcAlgTraits::nDim_, CURRENT_COORDINATES);
  elemDataPreReqs.add_gathered_nodal_field(*scalarQ_, 1);
  
  // manage master element requirements
  elemDataPreReqs.add_master_element_call(FEM_FACE_GRAD_OP, CURRENT_COORDINATES);
  faceDataPreReqs.add_master_element_call(FEM_FACE_NORMAL, CURRENT_COORDINATES);
  faceDataPreReqs.add_master_element_call(FEM_FACE_DET_J, CURRENT_COORDINATES);

  get_face_shape_fn_data<BcAlgTraits>([&](double* ptr){meFC->shape_fcn(ptr);}, vf_shape_function_);

  if ( shiftedGradOp_ )
    NaluEnv::self().naluOutputP0() << "ScalarFluxPenaltyFemKernel is not ready for shiftedGradOp" << std::endl;
}

template<typename BcAlgTraits>
ScalarFluxPenaltyFemKernel<BcAlgTraits>::~ScalarFluxPenaltyFemKernel()
{}

template<typename BcAlgTraits>
void
ScalarFluxPenaltyFemKernel<BcAlgTraits>::execute(
  SharedMemView<DoubleType**> &lhs,
  SharedMemView<DoubleType *> &rhs,
  ScratchViews<DoubleType> &faceScratchViews,
  ScratchViews<DoubleType> &elemScratchViews,
  int elemFaceOrdinal)
{
  NALU_ALIGNED DoubleType w_dqdxBip[BcAlgTraits::nDim_];
  
  const int *face_node_ordinals = meFEM_->side_node_ordinals(elemFaceOrdinal);
 
  // face
  SharedMemView<DoubleType*>& vf_scalarQ = faceScratchViews.get_scratch_view_1D(*scalarQ_);
  SharedMemView<DoubleType*>& vf_bcScalarQ = faceScratchViews.get_scratch_view_1D(*bcScalarQ_);
  SharedMemView<DoubleType*>& vf_diffFluxCoeff = faceScratchViews.get_scratch_view_1D(*diffFluxCoeff_);
 
  // element
  SharedMemView<DoubleType*>& v_scalarQ = elemScratchViews.get_scratch_view_1D(*scalarQ_);

  // master element calls
  SharedMemView<DoubleType***>& v_dndx_fc_elem = elemScratchViews.get_me_views(CURRENT_COORDINATES).dndx_fc_elem;
  SharedMemView<DoubleType*>& vf_det_j = faceScratchViews.get_me_views(CURRENT_COORDINATES).det_j_fc;
  SharedMemView<DoubleType**>& vf_normal = faceScratchViews.get_me_views(CURRENT_COORDINATES).normal_fc_fem;

  for ( int ip = 0; ip < BcAlgTraits::numFaceIp_; ++ip) {
    
    // zero out vector quantities
    for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
      w_dqdxBip[j] = 0.0;
    }

    // form L^-1
    DoubleType inverseLengthScale = 0.0;
    for ( int ic = 0; ic < BcAlgTraits::nodesPerFace_; ++ic ) {
      const int faceNodeNumber = face_node_ordinals[ic];
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        inverseLengthScale += v_dndx_fc_elem(ip,faceNodeNumber,j)*vf_normal(ip,j);
      }
    }        
    
    // interpolate to bip
    DoubleType qBip = 0.0;
    DoubleType qbcBip = 0.0;
    DoubleType diffFluxCoeffBip = 0.0;
    for ( int ic = 0; ic < BcAlgTraits::nodesPerFace_; ++ic ) {
      const DoubleType r = vf_shape_function_(ip,ic);
      qBip += r*vf_scalarQ(ic);
      qbcBip += r*vf_bcScalarQ(ic);
      diffFluxCoeffBip += r*vf_diffFluxCoeff(ic);
    }
    
    // form dqdxBip
    for ( int ic = 0; ic < BcAlgTraits::nodesPerElement_; ++ic ) {
      const DoubleType qIc = v_scalarQ(ic);
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        w_dqdxBip[j] += v_dndx_fc_elem(ip,ic,j)*qIc;
      }
    }

    // form flux; -diffFluxCoeffBip*dqdxj*nj + penaltyFac*diffFluxCoeffBip*invL*(qBip - qbcBip)
    DoubleType flux = penaltyFac_*diffFluxCoeffBip*inverseLengthScale*(qBip - qbcBip);
    for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
      const DoubleType nj = vf_normal(ip,j);
      flux -= diffFluxCoeffBip*w_dqdxBip[j]*nj;
    }
    
    // start the assembly
    const DoubleType ipFactor = vf_det_j(ip)*vf_ip_weight_(ip);
    
    // row ir
    for ( int ir = 0; ir < BcAlgTraits::nodesPerFace_; ++ir) {
      
      const int lIr = face_node_ordinals[ir];
      const DoubleType wIr = vf_shape_function_(ip,ir);
      
      // residual
      rhs(lIr) -= wIr*flux*ipFactor;

      // face-based penalty
      for ( int ic = 0; ic < BcAlgTraits::nodesPerFace_; ++ic ) {
        const int faceNodeNumber = face_node_ordinals[ic];
        lhs(lIr,faceNodeNumber) += wIr*vf_shape_function_(ip,ic)*penaltyFac_*diffFluxCoeffBip*inverseLengthScale*ipFactor;
      }
      
      // element-based gradient
      for ( int ic = 0; ic < BcAlgTraits::nodesPerElement_; ++ic ) {
        DoubleType lhsFac = 0.0;
        for ( int j = 0; j < BcAlgTraits::nDim_; ++j )
          lhsFac += -v_dndx_fc_elem(ip,ic,j)*vf_normal(ip,j);
        lhs(lIr,ic) += wIr*diffFluxCoeffBip*lhsFac*ipFactor;
      }    
    }
  }
}

INSTANTIATE_KERNEL_FEM_FACE_ELEMENT(ScalarFluxPenaltyFemKernel);

}  // nalu
}  // sierra
