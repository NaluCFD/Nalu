/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/ContinuityAdvFemKernel.h"
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
ContinuityAdvFemKernel<AlgTraits>::ContinuityAdvFemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    meshMotion_(solnOpts.does_mesh_move()),
    shiftMdot_(solnOpts.cvfemShiftMdot_),
    shiftedGradOp_(solnOpts.get_shifted_grad_op("pressure")),
    reducedSensitivities_(solnOpts.cvfemReducedSensPoisson_),
    projTimeScale_(1.0)
{
  // Save of required fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  std::string velocityRTM_name = meshMotion_ ? "velocity_rtm" : "velocity";
  velocityRTM_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, velocityRTM_name);
  Gjp_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx");
  pressure_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure");
  ScalarFieldType *density = metaData.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "density");
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  coordinates_ = metaData.get_field<VectorFieldType>( stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  // extract master element
  MasterElement *meFEM = sierra::nalu::MasterElementRepo::get_fem_master_element(AlgTraits::topo_);
  
  // copy ip weights into our 1-d view
  for ( int k = 0; k < AlgTraits::numGp_; ++k )
    v_ip_weight_[k] = meFEM->weights_[k];

  dataPreReqs.add_fem_volume_me(meFEM);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*velocityRTM_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*pressure_, 1);
  dataPreReqs.add_gathered_nodal_field(*Gjp_, AlgTraits::nDim_);
  
  // manage dndx
  if ( !shiftedGradOp_ || !reducedSensitivities_ )
    dataPreReqs.add_master_element_call(FEM_GRAD_OP, CURRENT_COORDINATES);
  if ( shiftedGradOp_ || reducedSensitivities_ )
    dataPreReqs.add_master_element_call(FEM_SHIFTED_GRAD_OP, CURRENT_COORDINATES);
  
  // finally, mdot shifting
  if ( shiftMdot_ )
    get_fem_shape_fn_data<AlgTraits>([&](double* ptr){meFEM->shifted_shape_fcn(ptr);}, v_shape_function_);
  else
    get_fem_shape_fn_data<AlgTraits>([&](double* ptr){meFEM->shape_fcn(ptr);}, v_shape_function_);

  // sanity checks for un-supported options
  if ( reducedSensitivities_ )
    throw std::runtime_error("ContinuityAdvFemKernel is not ready for reducedSensitivities_");
  if ( shiftMdot_ && !shiftedGradOp_ )
    throw std::runtime_error("ContinuityAdvFemKernel is not ready for shiftMdot_ && !shiftedGradOp_");
  if ( !shiftMdot_ && shiftedGradOp_ )
    throw std::runtime_error("ContinuityAdvFemKernel is not ready for !shiftMdot_ && shiftedGradOp_");
}

template<typename AlgTraits>
ContinuityAdvFemKernel<AlgTraits>::~ContinuityAdvFemKernel()
{
  // nothing to do
}

template<typename AlgTraits>
void
ContinuityAdvFemKernel<AlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
  const double dt = timeIntegrator.get_time_step();
  const double gamma1 = timeIntegrator.get_gamma1();
  projTimeScale_ = dt/gamma1;
}

template<typename AlgTraits>
void
ContinuityAdvFemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType**>& lhs,
  SharedMemView<DoubleType*>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  // Work arrays (fixed size)
  NALU_ALIGNED DoubleType w_rhoUip [AlgTraits::nDim_];
  NALU_ALIGNED DoubleType w_GpIp[AlgTraits::nDim_];
  NALU_ALIGNED DoubleType w_dpdxIp [AlgTraits::nDim_];

  SharedMemView<DoubleType*>& v_densityNp1 = scratchViews.get_scratch_view_1D(*densityNp1_);
  SharedMemView<DoubleType*>& v_pressure = scratchViews.get_scratch_view_1D(*pressure_);

  SharedMemView<DoubleType**>& v_vrtm = scratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<DoubleType**>& v_Gp = scratchViews.get_scratch_view_2D(*Gjp_);

  SharedMemView<DoubleType*>& v_det_j = scratchViews.get_me_views(CURRENT_COORDINATES).det_j_fem;
  
  SharedMemView<DoubleType***>& v_dndx = shiftedGradOp_ ?
    scratchViews.get_me_views(CURRENT_COORDINATES).dndx_fem : scratchViews.get_me_views(CURRENT_COORDINATES).dndx_fem;
  SharedMemView<DoubleType***>& v_dndx_lhs = (shiftedGradOp_ || reducedSensitivities_)?
    scratchViews.get_me_views(CURRENT_COORDINATES).dndx_fem : scratchViews.get_me_views(CURRENT_COORDINATES).dndx_fem;
  
  for ( int ip = 0; ip < AlgTraits::numGp_; ++ip ) {

    // zero out and compute w_rhoUip
    for (int j = 0; j < AlgTraits::nDim_; ++j) {
      w_rhoUip[j] = 0.0;
      w_GpIp[j] = 0.0;
      w_dpdxIp[j] = 0.0;
    }

    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const DoubleType r = v_shape_function_(ip,ic);
      const DoubleType rhoIc = v_densityNp1(ic);
      const DoubleType pIc = v_pressure(ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        w_rhoUip[j] += r*rhoIc*v_vrtm(ic,j);
        w_GpIp[j] += r*v_Gp(ic,j);
        w_dpdxIp[j] += v_dndx(ip,ic,j)*pIc;
      }
    }
    
    // start the assembly
    const DoubleType ipFactor = v_det_j(ip)*v_ip_weight_(ip);
    
    // row ir
    for ( int ir = 0; ir < AlgTraits::nodesPerElement_; ++ir) {

      DoubleType rhsSum = 0.0;
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        rhsSum -= (w_rhoUip[j] -projTimeScale_*(w_dpdxIp[j] - w_GpIp[j]))*v_dndx(ip,ir,j);
      }

      // column ic
      for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
        
        DoubleType lhsSum = 0.0;
        for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
          const DoubleType facLhs = v_dndx_lhs(ip,ir,j)*v_dndx_lhs(ip,ic,j);
          lhsSum += facLhs;
        }
        lhs(ir,ic) += lhsSum*ipFactor;
      }
      rhs(ir) -= rhsSum*ipFactor/projTimeScale_;
    }
  }
}

INSTANTIATE_FEM_KERNEL(ContinuityAdvFemKernel);

}  // nalu
}  // sierra
