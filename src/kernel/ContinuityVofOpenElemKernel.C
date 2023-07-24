/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/ContinuityVofOpenElemKernel.h"
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
ContinuityVofOpenElemKernel<BcAlgTraits>::ContinuityVofOpenElemKernel(
  const stk::mesh::MetaData &metaData,
  const SolutionOptions &solnOpts,
  ElemDataRequests &faceDataPreReqs,
  ElemDataRequests &elemDataPreReqs)
  : Kernel(),
    shiftedGradOp_(solnOpts.get_shifted_grad_op("pressure")),
    reducedSensitivities_(solnOpts.cvfemReducedSensPoisson_),
    n_(solnOpts.localVofN_),
    m_(solnOpts.localVofM_),
    c_(solnOpts.localVofC_),
    meSCS_(sierra::nalu::MasterElementRepo::get_surface_master_element(BcAlgTraits::elemTopo_))
{
  if ( solnOpts.does_mesh_move())
    velocityRTM_ = metaData.get_field<double>(stk::topology::NODE_RANK, "velocity_rtm");
  else
    velocityRTM_ = metaData.get_field<double>(stk::topology::NODE_RANK, "velocity");
  Gpdx_ = metaData.get_field<double>(stk::topology::NODE_RANK, "dpdx");
  coordinates_ = metaData.get_field<double>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  pressure_ = metaData.get_field<double>(stk::topology::NODE_RANK, "pressure");
  pressureBc_ = metaData.get_field<double>(stk::topology::NODE_RANK, "pressure_bc");
  density_ = metaData.get_field<double>(stk::topology::NODE_RANK, "density");
  interfaceCurvature_ = metaData.get_field<double>(stk::topology::NODE_RANK, "interface_curvature");
  surfaceTension_ = metaData.get_field<double>(stk::topology::NODE_RANK, "surface_tension");
  vof_ = metaData.get_field<double>(stk::topology::NODE_RANK, "volume_of_fluid");
  exposedAreaVec_ = metaData.get_field<double>(metaData.side_rank(), "exposed_area_vector");
  dynamicPressure_ = metaData.get_field<double>(metaData.side_rank(), "dynamic_pressure");
  for (int i = 0; i < BcAlgTraits::nDim_; ++i)
    gravity_(i) = solnOpts.gravity_[i];

  if (solnOpts.buoyancyPressureStab_)
    buoyancyWeight_ = 1.0;
  else
    buoyancyWeight_ = 0.0;

  // extract master elements
  MasterElement* meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(BcAlgTraits::faceTopo_);
  
  // add master elements
  faceDataPreReqs.add_cvfem_face_me(meFC);
  elemDataPreReqs.add_cvfem_surface_me(meSCS_);

  // fields and data; face and then element
  faceDataPreReqs.add_gathered_nodal_field(*pressure_, 1);
  faceDataPreReqs.add_gathered_nodal_field(*pressureBc_, 1);
  faceDataPreReqs.add_gathered_nodal_field(*density_, 1);
  faceDataPreReqs.add_gathered_nodal_field(*interfaceCurvature_, 1);
  faceDataPreReqs.add_gathered_nodal_field(*surfaceTension_, 1);
  faceDataPreReqs.add_gathered_nodal_field(*velocityRTM_, BcAlgTraits::nDim_);
  faceDataPreReqs.add_gathered_nodal_field(*Gpdx_, BcAlgTraits::nDim_);  
  faceDataPreReqs.add_gathered_nodal_field(*vof_, 1);
  faceDataPreReqs.add_face_field(*exposedAreaVec_, BcAlgTraits::numFaceIp_, BcAlgTraits::nDim_);
  faceDataPreReqs.add_face_field(*dynamicPressure_, BcAlgTraits::numFaceIp_);
  faceDataPreReqs.add_coordinates_field(*coordinates_, BcAlgTraits::nDim_, CURRENT_COORDINATES);
  
  elemDataPreReqs.add_coordinates_field(*coordinates_, BcAlgTraits::nDim_, CURRENT_COORDINATES);
  elemDataPreReqs.add_gathered_nodal_field(*pressure_, 1);
  elemDataPreReqs.add_gathered_nodal_field(*vof_, 1);

  // manage dndx
  if ( !shiftedGradOp_ || !reducedSensitivities_ )
    elemDataPreReqs.add_master_element_call(SCS_FACE_GRAD_OP, CURRENT_COORDINATES);
  if ( shiftedGradOp_ || reducedSensitivities_ )
    elemDataPreReqs.add_master_element_call(SCS_SHIFTED_FACE_GRAD_OP, CURRENT_COORDINATES);

  if ( solnOpts.cvfemShiftMdot_ )
    get_face_shape_fn_data<BcAlgTraits>([&](double* ptr){meFC->shifted_shape_fcn(ptr);}, vf_shape_function_);
  else
    get_face_shape_fn_data<BcAlgTraits>([&](double* ptr){meFC->shape_fcn(ptr);}, vf_shape_function_);

  NaluEnv::self().naluOutputP0() << "VOF ContinuityVofOpenElemKernel is active" << std::endl;
}

template<typename BcAlgTraits>
ContinuityVofOpenElemKernel<BcAlgTraits>::~ContinuityVofOpenElemKernel()
{}

template<typename BcAlgTraits>
void
ContinuityVofOpenElemKernel<BcAlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
  const double dt = timeIntegrator.get_time_step();
  const double gamma1 = timeIntegrator.get_gamma1();
  projTimeScale_ = dt/gamma1;
}

template<typename BcAlgTraits>
void
ContinuityVofOpenElemKernel<BcAlgTraits>::execute(
  SharedMemView<DoubleType**> &lhs,
  SharedMemView<DoubleType *> &rhs,
  ScratchViews<DoubleType> &faceScratchViews,
  ScratchViews<DoubleType> &elemScratchViews,
  int elemFaceOrdinal)
{
  NALU_ALIGNED DoubleType w_uBip[BcAlgTraits::nDim_];
  NALU_ALIGNED DoubleType w_GpdxBip[BcAlgTraits::nDim_];
  NALU_ALIGNED DoubleType w_dpdxBip[BcAlgTraits::nDim_];
 
  const int *face_node_ordinals = meSCS_->side_node_ordinals(elemFaceOrdinal);
 
  // face
  SharedMemView<DoubleType*>& vf_pressure = faceScratchViews.get_scratch_view_1D(*pressure_);
  SharedMemView<DoubleType*>& vf_pressureBc = faceScratchViews.get_scratch_view_1D(*pressureBc_);
  SharedMemView<DoubleType**>& vf_Gpdx = faceScratchViews.get_scratch_view_2D(*Gpdx_);
  SharedMemView<DoubleType*>& vf_density = faceScratchViews.get_scratch_view_1D(*density_);
  SharedMemView<DoubleType*>& vf_kappa = faceScratchViews.get_scratch_view_1D(*interfaceCurvature_);
  SharedMemView<DoubleType*>& vf_sigma = faceScratchViews.get_scratch_view_1D(*surfaceTension_);
  SharedMemView<DoubleType**>& vf_vrtm = faceScratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<DoubleType**>& vf_exposedAreaVec = faceScratchViews.get_scratch_view_2D(*exposedAreaVec_);
  SharedMemView<DoubleType*>& vf_dynamicP = faceScratchViews.get_scratch_view_1D(*dynamicPressure_);
  SharedMemView<DoubleType*>& vf_vof = faceScratchViews.get_scratch_view_1D(*vof_);
 
  // element
  SharedMemView<DoubleType*>& v_pressure = elemScratchViews.get_scratch_view_1D(*pressure_);
  SharedMemView<DoubleType*>& v_vof = elemScratchViews.get_scratch_view_1D(*vof_);

  // dndx for both rhs and lhs
  SharedMemView<DoubleType***>& v_dndx_fc_elem = shiftedGradOp_ 
    ? elemScratchViews.get_me_views(CURRENT_COORDINATES).dndx_shifted_fc_elem
    : elemScratchViews.get_me_views(CURRENT_COORDINATES).dndx_fc_elem;
  SharedMemView<DoubleType***>& v_dndx_fc_elem_lhs = (shiftedGradOp_ || reducedSensitivities_)
    ? elemScratchViews.get_me_views(CURRENT_COORDINATES).dndx_shifted_fc_elem
    : elemScratchViews.get_me_views(CURRENT_COORDINATES).dndx_fc_elem;

  for (int ip=0; ip < BcAlgTraits::numFaceIp_; ++ip) {
    
    const int nearestNode = meSCS_->ipNodeMap(elemFaceOrdinal)[ip]; // "Right"
    
    // zero out vector quantities; form aMag
    DoubleType aMag = 0.0;
    for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
      w_uBip[j] = 0.0;
      w_GpdxBip[j] = 0.0;
      w_dpdxBip[j] = 0.0;
      const DoubleType axj = vf_exposedAreaVec(ip,j);
      aMag += axj*axj;
    }
    aMag = stk::math::sqrt(aMag);
    
    // form L^-1
    DoubleType inverseLengthScale = 0.0;
    for ( int ic = 0; ic < BcAlgTraits::nodesPerFace_; ++ic ) {
      const int faceNodeNumber = face_node_ordinals[ic];
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        inverseLengthScale += v_dndx_fc_elem(ip,faceNodeNumber,j)*vf_exposedAreaVec(ip,j);
      }
    }        
    inverseLengthScale /= aMag;

    // interpolate to bip
    DoubleType pBip = 0.0;
    DoubleType rhoBip = 0.0;
    DoubleType pbcBip = -vf_dynamicP(ip);
    DoubleType sigmaKappaBip = 0.0;
    DoubleType vofBip = 0.0;
    for ( int ic = 0; ic < BcAlgTraits::nodesPerFace_; ++ic ) {
      const DoubleType r = vf_shape_function_(ip,ic);
      pBip += r*vf_pressure(ic);
      pbcBip += r*vf_pressureBc(ic);
      rhoBip += r*vf_density(ic);
      sigmaKappaBip += r*vf_sigma(ic)*vf_kappa(ic);
      vofBip += r*vf_vof(ic);
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        w_uBip[j] += r*vf_vrtm(ic,j);
        w_GpdxBip[j] += r*vf_Gpdx(ic,j);
      }
    }
    
    // form dpdxBip
    DoubleType dvofdaBip = 0.0;
    for ( int ic = 0; ic < BcAlgTraits::nodesPerElement_; ++ic ) {
      const DoubleType pIc = v_pressure(ic);
      const DoubleType vofIc = v_vof(ic);
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        const DoubleType dxj = v_dndx_fc_elem(ip,ic,j);
        w_dpdxBip[j] += dxj*pIc;
        dvofdaBip += dxj*vofIc*vf_exposedAreaVec(ip,j);
      }
    }

    // correct for localized approach
    sigmaKappaBip *= c_*stk::math::pow(vofBip,n_)*stk::math::pow(1.0-vofBip,m_);

    // form vdot; uj*Aj - projTS*(dpdxj/rhoBip - Gjph)*Aj + penaltyFac*projTS/rhoBip*invL*(pBip - pbcBip)*aMag + BF
    DoubleType vdot = penaltyFac_*projTimeScale_/rhoBip*inverseLengthScale*(pBip - pbcBip)*aMag
      + projTimeScale_*sigmaKappaBip*dvofdaBip/rhoBip;
    for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
      const DoubleType axj = vf_exposedAreaVec(ip,j);
      vdot += (w_uBip[j] - projTimeScale_*((w_dpdxBip[j]-buoyancyWeight_*rhoBip*gravity_(j))/rhoBip - w_GpdxBip[j]))*axj;
    }
    
    // face-based penalty; divide by projTimeScale
    for ( int ic = 0; ic < BcAlgTraits::nodesPerFace_; ++ic ) {
      const int faceNodeNumber = face_node_ordinals[ic];
      const DoubleType r = vf_shape_function_(ip,ic);
      lhs(nearestNode,faceNodeNumber) += r*penaltyFac_/rhoBip*inverseLengthScale*aMag;
    }
    
    // element-based gradient; divide by projTimeScale
    for ( int ic = 0; ic < BcAlgTraits::nodesPerElement_; ++ic ) {
      DoubleType lhsFac = 0.0;
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j )
        lhsFac += -v_dndx_fc_elem_lhs(ip,ic,j)*vf_exposedAreaVec(ip,j);
      lhs(nearestNode,ic) += lhsFac/rhoBip;
    }
    
    // residual
    rhs(nearestNode) -= vdot/projTimeScale_;
  }
}

INSTANTIATE_KERNEL_FACE_ELEMENT(ContinuityVofOpenElemKernel);

}  // nalu
}  // sierra
