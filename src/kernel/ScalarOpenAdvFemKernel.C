/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/ScalarOpenAdvFemKernel.h"
#include "EquationSystem.h"
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
ScalarOpenAdvFemKernel<BcAlgTraits>::ScalarOpenAdvFemKernel(
  const stk::mesh::MetaData &metaData,
  const SolutionOptions &solnOpts,
  ScalarFieldType *scalarQ,
  ScalarFieldType *bcScalarQ,
  ElemDataRequests &faceDataPreReqs,
  ElemDataRequests &elemDataPreReqs)
  : Kernel(),
    scalarQ_(scalarQ),
    bcScalarQ_(bcScalarQ),
    shiftedGradOp_(solnOpts.get_shifted_grad_op(scalarQ->name())),
    projTimeScale_(1.0),
    penaltyFac_(4.0),
    meFEM_(sierra::nalu::MasterElementRepo::get_fem_master_element(BcAlgTraits::elemTopo_))
{
  // save off fields
  vrtmL_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "vrtm_lagged");
  GjpL_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx_lagged");
  pressure_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure");
  pressureBc_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure_bc");
  density_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  dynamicPressure_ = metaData.get_field<GenericFieldType>(metaData.side_rank(), "dynamic_pressure");
    
  // extract field not required in execute()
  VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

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
  faceDataPreReqs.add_gathered_nodal_field(*vrtmL_, BcAlgTraits::nDim_);
  faceDataPreReqs.add_gathered_nodal_field(*GjpL_, BcAlgTraits::nDim_);
  faceDataPreReqs.add_gathered_nodal_field(*pressure_, 1);
  faceDataPreReqs.add_gathered_nodal_field(*pressureBc_, 1);
  faceDataPreReqs.add_gathered_nodal_field(*density_, 1);
  faceDataPreReqs.add_face_field(*dynamicPressure_, BcAlgTraits::numFaceIp_);
  faceDataPreReqs.add_coordinates_field(*coordinates, BcAlgTraits::nDim_, CURRENT_COORDINATES);

  elemDataPreReqs.add_gathered_nodal_field(*pressure_, 1);
  elemDataPreReqs.add_coordinates_field(*coordinates, BcAlgTraits::nDim_, CURRENT_COORDINATES);
  
  // manage master element requirements
  elemDataPreReqs.add_master_element_call(FEM_FACE_GRAD_OP, CURRENT_COORDINATES);
  faceDataPreReqs.add_master_element_call(FEM_FACE_NORMAL, CURRENT_COORDINATES);
  faceDataPreReqs.add_master_element_call(FEM_FACE_DET_J, CURRENT_COORDINATES);

  get_face_shape_fn_data<BcAlgTraits>([&](double* ptr){meFC->shape_fcn(ptr);}, vf_shape_function_);

  if ( shiftedGradOp_ )
    NaluEnv::self().naluOutputP0() << "ScalarOpenAdvFemKernel is not ready for shiftedGradOp" << std::endl;
}

template<typename BcAlgTraits>
ScalarOpenAdvFemKernel<BcAlgTraits>::~ScalarOpenAdvFemKernel()
{}

template<typename BcAlgTraits>
void
ScalarOpenAdvFemKernel<BcAlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
  const double dt = timeIntegrator.get_time_step();
  const double gamma1 = timeIntegrator.get_gamma1();
  projTimeScale_ = dt/gamma1;
}

template<typename BcAlgTraits>
void
ScalarOpenAdvFemKernel<BcAlgTraits>::execute(
  SharedMemView<DoubleType**> &lhs,
  SharedMemView<DoubleType *> &rhs,
  ScratchViews<DoubleType> &faceScratchViews,
  ScratchViews<DoubleType> &elemScratchViews,
  int elemFaceOrdinal)
{
  NALU_ALIGNED DoubleType w_rho_vrtmBip[BcAlgTraits::nDim_];
  NALU_ALIGNED DoubleType w_GpdxBip[BcAlgTraits::nDim_];
  NALU_ALIGNED DoubleType w_dpdxBip[BcAlgTraits::nDim_];
 
  const int *face_node_ordinals = meFEM_->side_node_ordinals(elemFaceOrdinal);
 
  // face
  SharedMemView<DoubleType*>& vf_scalarQ = faceScratchViews.get_scratch_view_1D(*scalarQ_);
  SharedMemView<DoubleType*>& vf_bcScalarQ = faceScratchViews.get_scratch_view_1D(*bcScalarQ_);
  SharedMemView<DoubleType**>& vf_vrtm = faceScratchViews.get_scratch_view_2D(*vrtmL_);
  SharedMemView<DoubleType**>& vf_Gpdx = faceScratchViews.get_scratch_view_2D(*GjpL_);
  SharedMemView<DoubleType*>& vf_pressure = faceScratchViews.get_scratch_view_1D(*pressure_);
  SharedMemView<DoubleType*>& vf_pressureBc = faceScratchViews.get_scratch_view_1D(*pressureBc_);
  SharedMemView<DoubleType*>& vf_density = faceScratchViews.get_scratch_view_1D(*density_);
  SharedMemView<DoubleType*>& vf_dynamicP = faceScratchViews.get_scratch_view_1D(*dynamicPressure_);
   
  // element
  SharedMemView<DoubleType*>& v_pressure = elemScratchViews.get_scratch_view_1D(*pressure_);
  
  // master element calls
  SharedMemView<DoubleType***>& v_dndx_fc_elem = elemScratchViews.get_me_views(CURRENT_COORDINATES).dndx_fc_elem;
  SharedMemView<DoubleType*>& vf_det_j = faceScratchViews.get_me_views(CURRENT_COORDINATES).det_j_fc;
  SharedMemView<DoubleType**>& vf_normal = faceScratchViews.get_me_views(CURRENT_COORDINATES).normal_fc_fem;
    
  for (int ip = 0; ip < BcAlgTraits::numFaceIp_; ++ip) {

    const DoubleType ipFactor = vf_det_j(ip)*vf_ip_weight_(ip);
        
    // zero out vector quantities
    for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
      w_rho_vrtmBip[j] = 0.0;
      w_GpdxBip[j] = 0.0;
      w_dpdxBip[j] = 0.0;
    }
        
    // form L^-1
    DoubleType inverseLengthScale = 0.0;
    for ( int ic = 0; ic < BcAlgTraits::nodesPerFace_; ++ic ) {
      const int faceNodeNumber = face_node_ordinals[ic];
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        inverseLengthScale += v_dndx_fc_elem(ip,faceNodeNumber,j)*vf_normal(ip,j);
      }
    }
    
    // compute quantities at bip
    DoubleType qBip = 0.0;
    DoubleType qBipEntrain = 0.0;
    DoubleType pBip = 0.0;
    DoubleType pbcBip = -vf_dynamicP(ip);
    for ( int ic = 0; ic < BcAlgTraits::nodesPerFace_; ++ic ) {
      const DoubleType r = vf_shape_function_(ip,ic);
      qBip += r*vf_scalarQ(ic);
      qBipEntrain += r*vf_bcScalarQ(ic);
      pBip += r*vf_pressure(ic);
      pbcBip += r*vf_pressureBc(ic);
      const DoubleType rhoIc = vf_density(ic);
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        w_rho_vrtmBip[j] += r*rhoIc*vf_vrtm(ic,j);
        w_GpdxBip[j] += r*vf_Gpdx(ic,j);
      }
    }
           
    // form dpdxBip
    for ( int ic = 0; ic < BcAlgTraits::nodesPerElement_; ++ic ) {
      const DoubleType pIc = v_pressure(ic);
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        w_dpdxBip[j] += v_dndx_fc_elem(ip,ic,j)*pIc;
      }
    }
    
    // form mass flux:
    // [rho*uj - projT*(dpdxj - GjpL)]*nj + penaltyFac*projTimeScale*invL*(pBip - pbcBip)
    DoubleType massFlux = penaltyFac_*projTimeScale_*inverseLengthScale*(pBip - pbcBip);
    for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
      massFlux += (w_rho_vrtmBip[j] - projTimeScale_*(w_dpdxBip[j] - w_GpdxBip[j]))*vf_normal(ip,j);
    }
  
    //================================
    // advection first
    //================================
    
    // account for both cases, i.e., leaving or entering to avoid hard loop over SIMD length
    const DoubleType mdotLeaving = stk::math::if_then_else(massFlux > 0.0, 1.0, 0.0);
    const DoubleType om_mdotLeaving = 1.0 - mdotLeaving;    

    // total advection (leaving and entraining)
    const DoubleType advLeaving = massFlux*qBip*ipFactor;
    const DoubleType advEntraining = massFlux*qBipEntrain*ipFactor;
    
    // central factor leaving the domain
    const DoubleType facLeaving = massFlux*ipFactor*mdotLeaving;
    
    // row ir
    for ( int ir = 0; ir < BcAlgTraits::nodesPerFace_; ++ir) {
      
      const DoubleType wIr = vf_shape_function_(ip,ir);      
      const int indexR = face_node_ordinals[ir];
        
      //===================
      // flow is leaving
      //===================
      rhs(indexR) -= wIr*advLeaving*mdotLeaving;
        
      // implicit
      for ( int ic = 0; ic < BcAlgTraits::nodesPerFace_; ++ic ) {
        const int nn = face_node_ordinals[ic];
        lhs(indexR,nn) += wIr*vf_shape_function_(ip,ic)*facLeaving;
      }
        
      //===================
      // flow is extraining
      //===================
      rhs(indexR) -= wIr*advEntraining*om_mdotLeaving;
    }
  }
}

INSTANTIATE_KERNEL_FEM_FACE_ELEMENT(ScalarOpenAdvFemKernel);

}  // nalu
}  // sierra
