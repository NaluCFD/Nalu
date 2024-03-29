/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/MomentumOpenAdvDiffElemKernel.h"
#include "EquationSystem.h"
#include "master_element/MasterElement.h"
#include "PecletFunction.h"
#include "SolutionOptions.h"

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
MomentumOpenAdvDiffElemKernel<BcAlgTraits>::MomentumOpenAdvDiffElemKernel(
  const stk::mesh::MetaData &metaData,
  const SolutionOptions &solnOpts,
  EquationSystem *eqSystem,
  VectorFieldType *velocity,
  GenericFieldType *Gjui,
  ScalarFieldType *viscosity,
  ElemDataRequests &faceDataPreReqs,
  ElemDataRequests &elemDataPreReqs)
  : Kernel(),
    viscosity_(viscosity),
    Gjui_(Gjui),
    alphaUpw_(solnOpts.get_alpha_upw_factor("velocity")),
    om_alphaUpw_(1.0-alphaUpw_),
    hoUpwind_(solnOpts.get_upw_factor("velocity")),
    includeDivU_(solnOpts.includeDivU_),
    meshVelocityCorrection_(solnOpts.does_mesh_move() ? 1.0 : 0.0),
    shiftedGradOp_(solnOpts.get_shifted_grad_op(velocity->name())),
    faceIpNodeMap_(sierra::nalu::MasterElementRepo::get_surface_master_element(BcAlgTraits::faceTopo_)->ipNodeMap()),
    meSCS_(sierra::nalu::MasterElementRepo::get_surface_master_element(BcAlgTraits::elemTopo_)),
    pecletFunction_(eqSystem->create_peclet_function<DoubleType>(velocity->name()))
{
  // save off fields
  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  if ( solnOpts.does_mesh_move() ) {
    velocityRTM_ = metaData.get_field<double>(stk::topology::NODE_RANK, "velocity_rtm");
    meshVelocity_ = metaData.get_field<double>(stk::topology::NODE_RANK, "mesh_velocity");
  }
  else {
    velocityRTM_ = velocityNp1_;
    meshVelocity_ = velocityNp1_;
  }
  coordinates_ = metaData.get_field<double>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  density_ = metaData.get_field<double>(stk::topology::NODE_RANK, "density");
  exposedAreaVec_ = metaData.get_field<double>(metaData.side_rank(), "exposed_area_vector");
  openMassFlowRate_ = metaData.get_field<double>(metaData.side_rank(), "open_mass_flow_rate");
  velocityBc_ = metaData.get_field<double>(stk::topology::NODE_RANK, "open_velocity_bc");
  
  // extract master elements
  MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(BcAlgTraits::faceTopo_);
  
  // add master elements
  faceDataPreReqs.add_cvfem_face_me(meFC);
  elemDataPreReqs.add_cvfem_surface_me(meSCS_);

  // fields and data; face and then element
  faceDataPreReqs.add_gathered_nodal_field(*density_, 1);
  faceDataPreReqs.add_gathered_nodal_field(*viscosity_, 1);
  faceDataPreReqs.add_gathered_nodal_field(*velocityNp1_, BcAlgTraits::nDim_);
  faceDataPreReqs.add_gathered_nodal_field(*meshVelocity_, BcAlgTraits::nDim_);
  faceDataPreReqs.add_gathered_nodal_field(*velocityBc_, BcAlgTraits::nDim_);
  faceDataPreReqs.add_gathered_nodal_field(*Gjui_, BcAlgTraits::nDim_, BcAlgTraits::nDim_);
  faceDataPreReqs.add_coordinates_field(*coordinates_, BcAlgTraits::nDim_, CURRENT_COORDINATES);
  faceDataPreReqs.add_face_field(*exposedAreaVec_, BcAlgTraits::numFaceIp_, BcAlgTraits::nDim_);
  faceDataPreReqs.add_face_field(*openMassFlowRate_, BcAlgTraits::numFaceIp_);
  
  elemDataPreReqs.add_coordinates_field(*coordinates_, BcAlgTraits::nDim_, CURRENT_COORDINATES);
  elemDataPreReqs.add_gathered_nodal_field(*velocityNp1_, BcAlgTraits::nDim_);
  elemDataPreReqs.add_gathered_nodal_field(*velocityRTM_, BcAlgTraits::nDim_);
  elemDataPreReqs.add_gathered_nodal_field(*viscosity_, 1);
  elemDataPreReqs.add_gathered_nodal_field(*density_, 1);
 
  if ( shiftedGradOp_ )
    elemDataPreReqs.add_master_element_call(SCS_SHIFTED_FACE_GRAD_OP, CURRENT_COORDINATES);
  else
    elemDataPreReqs.add_master_element_call(SCS_FACE_GRAD_OP, CURRENT_COORDINATES);

  // never shift properties
  get_face_shape_fn_data<BcAlgTraits>([&](double* ptr){meFC->shape_fcn(ptr);}, vf_shape_function_);

  // advection operator may use shifted
  const bool skewSymmetric = solnOpts.get_skew_symmetric(velocity->name());
  get_face_shape_fn_data<BcAlgTraits>([&](double* ptr){skewSymmetric ? meFC->shifted_shape_fcn(ptr) : meFC->shape_fcn(ptr);}, 
                                      vf_adv_shape_function_);
}

template<typename BcAlgTraits>
MomentumOpenAdvDiffElemKernel<BcAlgTraits>::~MomentumOpenAdvDiffElemKernel()
{
  delete pecletFunction_;
}


template<typename BcAlgTraits>
void
MomentumOpenAdvDiffElemKernel<BcAlgTraits>::execute(
  SharedMemView<DoubleType**> &lhs,
  SharedMemView<DoubleType *> &rhs,
  ScratchViews<DoubleType> &faceScratchViews,
  ScratchViews<DoubleType> &elemScratchViews,
  int elemFaceOrdinal)
{
  NALU_ALIGNED DoubleType w_uBip[BcAlgTraits::nDim_];
  NALU_ALIGNED DoubleType w_rho_vBip[BcAlgTraits::nDim_];
  NALU_ALIGNED DoubleType w_uBipExtrap[BcAlgTraits::nDim_];
  NALU_ALIGNED DoubleType w_uspecBip[BcAlgTraits::nDim_];
  NALU_ALIGNED DoubleType w_coordBip[BcAlgTraits::nDim_];
  NALU_ALIGNED DoubleType w_nx[BcAlgTraits::nDim_];

  const int *face_node_ordinals = meSCS_->side_node_ordinals(elemFaceOrdinal);
 
  // face
  SharedMemView<DoubleType*>& vf_density = faceScratchViews.get_scratch_view_1D(*density_);
  SharedMemView<DoubleType*>& vf_viscosity = faceScratchViews.get_scratch_view_1D(*viscosity_);
  SharedMemView<DoubleType**>& vf_velocityNp1 = faceScratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<DoubleType**>& vf_meshVelocity = faceScratchViews.get_scratch_view_2D(*meshVelocity_);
  SharedMemView<DoubleType**>& vf_bcVelocity = faceScratchViews.get_scratch_view_2D(*velocityBc_);
  SharedMemView<DoubleType***>& vf_Gjui = faceScratchViews.get_scratch_view_3D(*Gjui_);
  SharedMemView<DoubleType**>& vf_coordinates = faceScratchViews.get_scratch_view_2D(*coordinates_);
  SharedMemView<DoubleType**>& vf_exposedAreaVec = faceScratchViews.get_scratch_view_2D(*exposedAreaVec_);
  SharedMemView<DoubleType*>& vf_openMassFlowRate = faceScratchViews.get_scratch_view_1D(*openMassFlowRate_);
  
  // element
  SharedMemView<DoubleType**>& v_coordinates = elemScratchViews.get_scratch_view_2D(*coordinates_);
  SharedMemView<DoubleType**>& v_velocityNp1 = elemScratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<DoubleType**>& v_vrtm = elemScratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<DoubleType*>& v_viscosity = elemScratchViews.get_scratch_view_1D(*viscosity_);
  SharedMemView<DoubleType*>& v_density = elemScratchViews.get_scratch_view_1D(*density_);
  SharedMemView<DoubleType***>& v_dndx_fc_elem = shiftedGradOp_
    ? elemScratchViews.get_me_views(CURRENT_COORDINATES).dndx_shifted_fc_elem
    : elemScratchViews.get_me_views(CURRENT_COORDINATES).dndx_fc_elem;

  for (int ip=0; ip < BcAlgTraits::numFaceIp_; ++ip) {
    
    const int opposingNode = meSCS_->opposingNodes(elemFaceOrdinal,ip); // "Left"
    const int nearestNode = meSCS_->ipNodeMap(elemFaceOrdinal)[ip]; // "Right"
    const int localFaceNode = faceIpNodeMap_[ip];
    
    // zero out vector quantities
    DoubleType asq = 0.0;
    for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
      w_uBip[j] = 0.0;
      w_rho_vBip[j] = 0.0;
      w_uspecBip[j] = 0.0;
      w_coordBip[j] = 0.0;
      const DoubleType axj = vf_exposedAreaVec(ip,j);
      asq += axj*axj;
    }
    const DoubleType amag = stk::math::sqrt(asq);
    
    // interpolate to bip
    DoubleType rhoBip = 0.0;
    DoubleType viscBip = 0.0;
    for ( int ic = 0; ic < BcAlgTraits::nodesPerFace_; ++ic ) {
      const DoubleType r = vf_shape_function_(ip,ic);
      const DoubleType rAdv = vf_adv_shape_function_(ip,ic);
      rhoBip += r*vf_density(ic);
      viscBip += r*vf_viscosity(ic);
      DoubleType rhoIc = vf_density(ic);
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        w_uspecBip[j] += rAdv*vf_bcVelocity(ic,j);
        w_uBip[j] += rAdv*vf_velocityNp1(ic,j);
        w_coordBip[j] += rAdv*vf_coordinates(ic,j);
        w_rho_vBip[j] += r*rhoIc*vf_meshVelocity(ic,j);
      }
    }
        
    // Peclet factor; along the edge is fine  
    DoubleType udotx = 0.0;
    for ( int i = 0; i < BcAlgTraits::nDim_; ++i ) {
      const DoubleType dxi = v_coordinates(nearestNode,i) - v_coordinates(opposingNode,i);
      udotx += 0.5*dxi*(v_vrtm(opposingNode,i) + v_vrtm(nearestNode,i));
      w_nx[i] = vf_exposedAreaVec(ip,i)/amag; 
      // extrapolation
      DoubleType duR = 0.0;
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        DoubleType dxj = w_coordBip[j] - v_coordinates(nearestNode,j);
        duR += dxj*vf_Gjui(localFaceNode,i,j)*hoUpwind_;
      }
      w_uBipExtrap[i] = v_velocityNp1(nearestNode,i) + duR;
    }
    
    // Peclet factor; along the edge is fine
    const DoubleType diffIp = 0.5*(v_viscosity(opposingNode)/v_density(opposingNode)
                                   + v_viscosity(nearestNode)/v_density(nearestNode));
    const DoubleType pecFuncArg = stk::math::abs(udotx)/(diffIp+small_);
    const DoubleType pecfac = pecletFunction_->execute(pecFuncArg);
    const DoubleType om_pecfac = 1.0-pecfac;
   
    //================================
    // advection first
    //================================
    const DoubleType tmdot = vf_openMassFlowRate(ip);

    // account for both cases, i.e., leaving or entering to avoid hard loop over SIMD length
    const DoubleType mdotLeaving = stk::math::if_then_else(tmdot > 0, 1.0, 0.0);
    const DoubleType om_mdotLeaving = 1.0 - mdotLeaving;    

    // entrainment magnitude (must correct for possible mesh motion at the open bc)
    DoubleType mvc = 0.0;
    for ( int j = 0; j < BcAlgTraits::nDim_; ++j )
      mvc += w_rho_vBip[j]*vf_exposedAreaVec(ip,j);
    const DoubleType uEntrain = tmdot/(rhoBip*amag) + mvc/(rhoBip*amag)*meshVelocityCorrection_;

    // user specified entrainment
    DoubleType uspecbipnx = 0.0;
    for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
      const DoubleType nj = w_nx[j];
      uspecbipnx += w_uspecBip[j]*nj;
    }

    for ( int i = 0; i < BcAlgTraits::nDim_; ++i ) {
      
      const int indexR = nearestNode*BcAlgTraits::nDim_ + i;
      
      //================
      // flow is leaving
      //================

      // central
      const DoubleType uiIp = w_uBip[i];

      // upwind
      const DoubleType uiUpwind = alphaUpw_*w_uBipExtrap[i] + om_alphaUpw_*uiIp;

      // total advection; pressure contribution in time expression
      const DoubleType afluxLeaving = tmdot*(pecfac*uiUpwind+om_pecfac*uiIp);

      rhs(indexR) -= afluxLeaving*mdotLeaving;

      // upwind lhs
      lhs(indexR,nearestNode*BcAlgTraits::nDim_+i) += tmdot*pecfac*alphaUpw_*mdotLeaving;

      // central part
      const DoubleType fac = tmdot*(pecfac*om_alphaUpw_+om_pecfac)*mdotLeaving;
      for ( int ic = 0; ic < BcAlgTraits::nodesPerFace_; ++ic ) {
        const int nn = face_node_ordinals[ic];
        lhs(indexR,nn*BcAlgTraits::nDim_+i) += vf_adv_shape_function_(ip,ic)*fac;
      }

      //===================
      // flow is extraining
      //===================
      
      const DoubleType nxi = w_nx[i];
      
      // total advection; with normal and tangeant entrainment
      const DoubleType afluxEntraining = tmdot*uEntrain*nxi
        + tmdot*(w_uspecBip[i] - uspecbipnx*nxi);
      
      rhs(indexR) -= afluxEntraining*om_mdotLeaving;
    }
    
    //================================
    // diffusion second
    //================================
    for ( int ic = 0; ic < BcAlgTraits::nodesPerElement_; ++ic ) {
      
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        
        const DoubleType axj = vf_exposedAreaVec(ip,j);
        const DoubleType dndxj = v_dndx_fc_elem(ip,ic,j);
        const DoubleType uxj = v_velocityNp1(ic,j);
        
        const DoubleType divUstress = 2.0/3.0*viscBip*dndxj*uxj*axj*includeDivU_;
        
        for ( int i = 0; i < BcAlgTraits::nDim_; ++i ) {
          
          // matrix entries
          int indexR = nearestNode*BcAlgTraits::nDim_ + i;
          
          const DoubleType dndxi = v_dndx_fc_elem(ip,ic,i);
          const DoubleType uxi = v_velocityNp1(ic,i);
          const DoubleType nxi = w_nx[i];
          const DoubleType om_nxinxi = 1.0-nxi*nxi;
          
          // -mu*dui/dxj*Aj(1.0-nini); sneak in divU (explicit)
          DoubleType lhsfac = -viscBip*dndxj*axj*om_nxinxi;
          lhs(indexR,ic*BcAlgTraits::nDim_+i) += lhsfac;
          rhs(indexR) -= lhsfac*uxi + divUstress*om_nxinxi;
          
          // -mu*duj/dxi*Aj(1.0-nini)
          lhsfac = -viscBip*dndxi*axj*om_nxinxi;
          lhs(indexR,ic*BcAlgTraits::nDim_+j) += lhsfac;
          rhs(indexR) -= lhsfac*uxj;
          
          // now we need the -nx*ny*Fy - nx*nz*Fz part
          for ( int l = 0; l < BcAlgTraits::nDim_; ++l ) {
            
            if ( i != l ) {
              const DoubleType nxinxl = nxi*w_nx[l];
              const DoubleType uxl = v_velocityNp1(ic,l);
              const DoubleType dndxl = v_dndx_fc_elem(ip,ic,l);
              
              // +ni*nl*mu*dul/dxj*Aj; sneak in divU (explicit)
              lhsfac = viscBip*dndxj*axj*nxinxl;
              lhs(indexR,ic*BcAlgTraits::nDim_+l) += lhsfac;
              rhs(indexR) -= lhsfac*uxl + divUstress*nxinxl;
              
              // +ni*nl*mu*duj/dxl*Aj
              lhsfac = viscBip*dndxl*axj*nxinxl;
              lhs(indexR,ic*BcAlgTraits::nDim_+j) += lhsfac;
              rhs(indexR) -= lhsfac*uxj;
            }
          }
        }
      }
    }
  }
}

INSTANTIATE_KERNEL_FACE_ELEMENT(MomentumOpenAdvDiffElemKernel);

}  // nalu
}  // sierra
