/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/MomentumOpenAdvDiffElemKernel.h"
#include "BcAlgTraits.h"
#include "master_element/MasterElement.h"
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
  VectorFieldType *velocity,
  GenericFieldType *Gjui,
  ScalarFieldType *viscosity,
  ElemDataRequests &faceDataPreReqs,
  ElemDataRequests &elemDataPreReqs)
  : Kernel(),
    viscosity_(viscosity),
    Gjui_(Gjui),
    alphaUpw_(solnOpts..get_alpha_upw_factor("velocity")),
    om_alphaUpw_(1.0-alphaUpw_),
    hoUpwind_(solnOpts.get_upw_factor("velocity")),
    nfEntrain_(solnOpts.nearestFaceEntrain_),
    om_nfEntrain_(1.0-nfEntrain_),
    includeDivU_(solnOpts.includeDivU_),
    shiftedGradOp_(solnOpts.get_shifted_grad_op(velocity->name())),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_surface_master_element(BcAlgTraits::elemTopo_)->ipNodeMap()),
    faceIpNodeMap_(sierra::nalu::MasterElementRepo::get_face_master_element(BcAlgTraits::faceTopo_)->ipNodeMap()),
{
  // save off fields
  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  if ( solnOpts.does_mesh_move() )
    velocityRTM_ = msetaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = velocity;
  coordinates_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  density_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  exposedAreaVec_ = metaData.get_field<GenericFieldType>(meta_data.side_rank(), "exposed_area_vector");
  openMassFlowRate_ = metaData.get_field<GenericFieldType>(meta_data.side_rank(), "open_mass_flow_rate");
  velocityBc_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "open_velocity_bc");
  
  // create the peclet blending function
  pecletFunction_ = eqSystem->create_peclet_function(velocity_->name());

  // extract master elements
  MasterElement* meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(BcAlgTraits::faceTopo_);
  MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(BcAlgTraits::elemTopo_);
  
  // add master elements
  faceDataPreReqs.add_cvfem_face_me(meFC);
  elemDataPreReqs.add_cvfem_surface_me(meSCS);

  // fields and data; face and then element
  faceDataPreReqs.add_gathered_nodal_field(*viscosity_, 1);
  faceDataPreReqs.add_gathered_nodal_field(*velocityNp1_, BcAlgTraits::nDim_);
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
  get_scs_shape_fn_data<BcAlgTraits>([&](double* ptr){meSCS->shape_fcn(ptr);}, v_shape_function_);
}

template<typename BcAlgTraits>
MomentumOpenAdvDiffElemKernel<BcAlgTraits>::~MomentumOpenAdvDiffElemKernel()
{}

template<typename BcAlgTraits>
void
MomentumOpenAdvDiffElemKernel<BcAlgTraits>::execute(
  SharedMemView<DoubleType**> &lhs,
  SharedMemView<DoubleType *> &rhs,
  ScratchViews<DoubleType> &faceScratchViews,
  ScratchViews<DoubleType> &elemScratchViews)
{
  DoubleType w_uBip[BcAlgTraits::nDim_];
  DoubleType w_uScs[BcAlgTraits::nDim_];
  DoubleType w_uBipExtrap[BcAlgTraits::nDim_];
  DoubleType w_uspecBip[BcAlgTraits::nDim_];
  DoubleType w_coordBip[BcAlgTraits::nDim_];
  DoubleType w_nx[BcAlgTraits::nDim_];

  // face
  SharedMemView<DoubleType*>& vf_viscosity = faceScratchViews.get_scratch_view_1D(*viscosity_);
  SharedMemView<DoubleType**>& vf_velocityNp1 = faceScratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<DoubleType**>& vf_bcVelocity = faceScratchViews.get_scratch_view_2D(*velocityBc_);
  SharedMemView<DoubleType***>& vf_Gjui = faceScratchViews.get_scratch_view_3D(*Gjui_);
  SharedMemView<DoubleType**>& vf_coordinates = faceScratchViews.get_scratch_view_2D(*coordinates_);
  SharedMemView<DoubleType**>& vf_exposedAreaVec = faceScratchViews.get_scratch_view_2D(*exposedAreaVec_);
  SharedMemView<DoubleType*>& vf_openMassFlowRate = faceScratchViews.get_scratch_view_1D(*openMassFlowRate_);
  
  // element
  SharedMemView<DoubleType**>& v_coordinates = elemScratchViews.get_scratch_view_2D(*coordinates_);
  SharedMemView<DoubleType**>& v_velocityNp1 = elemScratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<DoubleType**>& v_vrtm = elemScratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<DoubleType**>& v_viscosity = elemScratchViews.get_scratch_view_1D(*viscosity_);
  SharedMemView<DoubleType**>& v_density = elemScratchViews.get_scratch_view_1D(*density_);
  SharedMemView<DoubleType***>& v_dndx = shiftedGradOp_
    ? elemScratchViews.get_me_views(CURRENT_COORDINATES).dndx_shifted_fc_scs
    : elemScratchViews.get_me_views(CURRENT_COORDINATES).dndx_fc_scs;

  // fixme: need master elements
  MasterElement* meFC = NULL;
  MasterElement* meSCS = NULL;

  for (int ip=0; ip < BcAlgTraits::numFaceIp_; ++ip) {
    
    // fixme: need face_ordinal, hack it in...
    const int face_ordinal = -1;

    // fixme: need the master element
    const int opposingNode = meSCS->opposingNodes(face_ordinal,ip);
    const int nearestNode = ipNodeMap_[ip];
    const int opposingScsIp = meSCS->opposingFace(face_ordinal,ip);
    const int localFaceNode = faceIpNodeMap[ip];
    
    // fixme: left and right nodes; right is on the face; left is the opposing node
    stk::mesh::Entity nodeL = stk::mesh::Entity(); // elem_node_rels[opposingNode];
    stk::mesh::Entity nodeR = stk::mesh::Entity(); // elem_node_rels[nearestNode];
    
    // zero out vector quantities
    DoubleType asq = 0.0;
    for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
      w_uBip[j] = 0.0;
      w_uScs[j] = 0.0;
      w_uspecBip[j] = 0.0;
      w_coordBip[j] = 0.0;
      const double axj = vf_exposedAreaVec(ip,j);
      asq += axj*axj;
    }
    const DoubleType amag = std::sqrt(asq);
    
    // interpolate to bip
    DoubleType viscBip = 0.0;
    for ( int ic = 0; ic < NcAlgTraits::nodesPerFace_; ++ic ) {
      const DoubleType r = vf_shape_function_(ip,ic);
      viscBip += r*vf_viscosity(ic);
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        w_uspecBip[j] += r*vf_bcVelocity(ic,j);
        w_uBip[j] += r*vf_velocityNp1(ic,j);
        w_coordBip[j] += r*vf_coordinates(ic,j);
      }
    }
    
    // data at interior opposing face
    for ( int ic = 0; ic < BcAlgtraits::nodesPerElement_; ++ic ) {
      const DoubleType r = v_shape_function(ip,ic);
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        p_uScs[j] += r*v_velocityNp1(ic,j);
      }
    }
    
    // fixme: is this okay? Peclet factor; along the edge is fine
    const DoubleType densL   = *stk::mesh::field_data(densityNp1, nodeL);
    const DoubleType densR   = *stk::mesh::field_data(densityNp1, nodeR);
    const DoubleType viscL   = *stk::mesh::field_data(*viscosity_, nodeL);
    const DoubleType viscR   = *stk::mesh::field_data(*viscosity_, nodeR);
    const DoubleType *uNp1R  =  stk::mesh::field_data(velocityNp1, nodeR);
    const DoubleType *vrtmL  =  stk::mesh::field_data(*velocityRTM_, nodeL);
    const DoubleType *vrtmR  =  stk::mesh::field_data(*velocityRTM_, nodeR);
    
    const DoubleType *coordL =  stk::mesh::field_data(*coordinates_, nodeL);
    const DoubleType *coordR =  stk::mesh::field_data(*coordinates_, nodeR);
    
    DoubleType udotx = 0.0;
    for ( int i = 0; i < BcAlgTraits::nDim_; ++i ) {
      const DoubleType dxi = coordR[i]  - coordL[i];
      udotx += 0.5*dxi*(vrtmL[i] + vrtmR[i]);
      w_nx[i] = vf_exposedAreaVec(ip,i)/amag;
      // extrapolation
      DoubleType duR = 0.0;
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        DoubleType dxj = w_coordBip[j] - coordR[j];
        duR += dxj*vf_Gjui(localFaceNode,i,j)*hoUpwind;
      }
      w_uBipExtrap[i] = uNp1R[i] + duR;
    }
    
    const DoubleType diffIp = 0.5*(viscL/densL + viscR/densR);
    const DoubleType pecFuncArg = stk::math::abs(udotx)/(diffIp+small_);
    // FIXME: modify pecletFunction to be double type
    DoubleType pecfac = 0.0;
    for(int simdIndex=0; simdIndex<stk::simd::ndoubles; ++simdIndex) {
      stk::simd::set_data(pecfac, simdIndex, pecletFunction_->execute(stk::simd::get_data(pecFuncArg, simdIndex)));
    }
    const DoubleType om_pecfac = 1.0-pecfac;
   
    //================================
    // advection first
    //================================
    const DoubleType tmdot = vf_openMassFlowRate(ip);
    
    // advection; leaving the domain
    // fixme... oh my, do I need to wrap this in a stk::math::if_then_else(tmdot > 0, doThis, doThat);
    if ( tmdot > 0.0 ) {
      
      for ( int i = 0; i < BcAlgTraits::nDim_; ++i ) {
        
        const int indexR = nearestNode*BcAlgTraits::nDim_ + i;
        
        // central
        const DoubleType uiIp = w_uBip[i];
        
        // upwind
        const DoubleType uiUpwind = alphaUpw*w_uBipExtrap[i] + om_alphaUpw*uiIp;
        
        // total advection; pressure contribution in time expression
        const DoubleType aflux = tmdot*(pecfac*uiUpwind+om_pecfac*uiIp);
        
        rhs(indexR) -= aflux;
        
        // upwind lhs
        lhs(indexR, nearestNode*BcAlgTraits::nDim_+i) += tmdot*pecfac*alphaUpw;
        
        // central part
        const DoubleType fac = tmdot*(pecfac*om_alphaUpw+om_pecfac);
        for ( int ic = 0; ic < NcAlgTraits::nodesPerFace_; ++ic ) {
          const DoubleType r = vf_shape_function(ip,ic);
          // fixme: where to extract face_node_ordinals
          const int nn = face_node_ordinals[ic];
          lhs(indexR, nn*BcAlgTraits::nDim_+i) += r*fac;
        }
      }
    }
    else {
      // extrainment
      DoubleType ubipnx = 0.0;
      DoubleType ubipExtrapnx = 0.0;
      DoubleType uscsnx = 0.0;
      DoubleType uspecbipnx = 0.0;
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        const DoubleType nj = w_nx[j];
        ubipnx += w_uBip[j]*nj;
        ubipExtrapnx += w_uBipExtrap[j]*nj;
        uscsnx += w_uScs[j]*nj;
        uspecbipnx += w_uspecBip[j]*nj;
      }
      
      for ( int i = 0; i < BcAlgTraits::nDim_; ++i ) {
        const int indexR = nearestNode*BcAlgTraits::nDim_ + i;
        const DoubleType nxi = w_nx[i];
        
        // total advection; with tangeant entrain
        const DoubleType aflux = tmdot*(pecfac*ubipExtrapnx+om_pecfac*
                                        (nfEntrain*ubipnx + om_nfEntrain*uscsnx))*nxi
          + tmdot*(w_uspecBip[i] - uspecbipnx*nxi);
        
        rhs(indexR) -= aflux;
        
        // upwind and central
        for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
          const DoubleType nxinxj = nxi*w_nx[j];
          
          // upwind
          lhs(indexR,nearestNode*BcAlgTraits::nDim_+j) += tmdot*pecfac*alphaUpw*nxinxj;
          
          // central part; exposed face
          DoubleType fac = tmdot*(pecfac*om_alphaUpw+om_pecfac*nfEntrain)*nxinxj;
          for ( int ic = 0; ic < NcAlgTraits::nodesPerFace_; ++ic ) {
            const DoubleType r = vf_shape_function(ip,ic);
            // fixme: where to extract face_node_ordinals
            
            const int nn = face_node_ordinals[ic];
            lhs(indexR,nn*BcAlgTraits::nDim_+j) += r*fac;
          }
          
          // central part; scs face
          fac = tmdot*om_pecfac*om_nfEntrain*nxinxj;
          for ( int ic = 0; ic < BcAlgtraits::nodesPerElement_; ++ic ) {
            const DoubleType r = p_shape_function[offSetSF_elem+ic];
            lhs(indexR, ic*BcAlgTraits::nDim_+j) += r*fac;
          }
          
        }
      }
    }
    
    //================================
    // diffusion second
    //================================
    for ( int ic = 0; ic < BcAlgtraits::nodesPerElement_; ++ic ) {
      
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        
        const DoubleType axj = vf_exposedAreaVec(ip,j);
        const DoubleType dndxj = v_dndx(ip,ic,j);
        const DoubleType uxj = v_velocityNp1(ic,j);
        
        const DoubleType divUstress = 2.0/3.0*viscBip*dndxj*uxj*axj*includeDivU_;
        
        for ( int i = 0; i < BcAlgTraits::nDim_; ++i ) {
          
          // matrix entries
          int indexR = nearestNode*BcAlgTraits::nDim_ + i;
          
          const DoubleType dndxi = v_dndx(ip,ic,i);
          const DoubleType uxi = v_velocityNp1(ic,i);
          const DoubleType nxi = w_nx[i];
          const DoubleType om_nxinxi = 1.0-nxi*nxi;
          
          // -mu*dui/dxj*Aj(1.0-nini); sneak in divU (explicit)
          DoubleType lhsfac = -viscBip*dndxj*axj*om_nxinxi;
          lhs(indexR,ic*BcAlgTraits::nDim_+i) += lhsfac;
          rhs(indexR) -= lhsfac*uxi + divUstress*om_nxinxi;
          
          // -mu*duj/dxi*Aj(1.0-nini)
          lhsfac = -viscBip*dndxi*axj*om_nxinxi;
          lhs[indexR,ic*BcAlgTraits::nDim_+j) += lhsfac;
          rhs(indexR) -= lhsfac*uxj;
          
          // now we need the -nx*ny*Fy - nx*nz*Fz part
          for ( int l = 0; l < BcAlgTraits::nDim_; ++l ) {
            
            if ( i != l ) {
              const DoubleType nxinxl = nxi*w_nx[l];
              const DoubleType uxl = v_velocityNp1(ic,l);
              const DoubleType dndxl = v_dndx(ip,ic,l);
              
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

template class MomentumOpenAdvDiffElemKernel <BcAlgTraitsHex8Quad4>;

}  // nalu
}  // sierra
