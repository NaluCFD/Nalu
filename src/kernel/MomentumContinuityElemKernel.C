/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/MomentumContinuityElemKernel.h"
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

template<class AlgTraits>
MomentumContinuityElemKernel<AlgTraits>::MomentumContinuityElemKernel(
  const stk::mesh::BulkData &bulkData,
  const SolutionOptions &solnOpts,
  VectorFieldType *velocity,
  ScalarFieldType *density,
  ScalarFieldType *viscosity,
  const bool lumpedMass,
  ElemDataRequests &dataPreReqs)
  : Kernel(),
    viscosity_(viscosity),
    includeDivU_(solnOpts.includeDivU_),
    lrscv_(sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_)->adjacentNodes()),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap()),
    dofSize_(AlgTraits::nDim_+1)
{
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  densityN_ = &(density->field_of_state(stk::mesh::StateN));
  if (density->number_of_states() == 2)
    densityNm1_ = densityN_;
  else
    densityNm1_ = &(density->field_of_state(stk::mesh::StateNM1));

  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  velocityN_ = &(velocity->field_of_state(stk::mesh::StateN));
  if (velocity->number_of_states() == 2)
    velocityNm1_ = velocityN_;
  else
    velocityNm1_ = &(velocity->field_of_state(stk::mesh::StateNM1));
  
  Gjp_ = metaData.get_field<double>(
    stk::topology::NODE_RANK, "dpdx");
  GjpOld_ = metaData.get_field<double>(
    stk::topology::NODE_RANK, "dpdx_old");
  
  coordinates_ = metaData.get_field<double>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  pressure_ = metaData.get_field<double>(
    stk::topology::NODE_RANK, "pressure");

  // fields and data; mdot not gathered as element data
  dataPreReqs.add_gathered_nodal_field(*densityNm1_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityN_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*velocityNp1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityN_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*velocityNm1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*Gjp_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*GjpOld_, AlgTraits::nDim_);
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*pressure_, 1);
  dataPreReqs.add_gathered_nodal_field(*viscosity_, 1);

  // master element registrations; surface
  MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_);
  dataPreReqs.add_cvfem_surface_me(meSCS);
  dataPreReqs.add_master_element_call(SCS_AREAV, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(SCS_GRAD_OP, CURRENT_COORDINATES);

  // master element registrations; volume
  MasterElement* meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);
  dataPreReqs.add_cvfem_volume_me(meSCV);
  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);

  // compute shape functions
  get_scs_shape_fn_data<AlgTraits>([&](double* ptr){meSCS->shape_fcn(ptr);}, v_shape_function_scs_);
  if ( lumpedMass )
    get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shifted_shape_fcn(ptr);}, v_shape_function_scv_);
  else
    get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shape_fcn(ptr);}, v_shape_function_scv_);
  
  // error checks - we are not supporting skew symmetric or shifted grad-ops
  const bool skewSymmetric = solnOpts.get_skew_symmetric(velocity->name());
  const bool shiftedGradOpU = solnOpts.get_shifted_grad_op(velocity->name());
  const bool shiftedGradOpP = solnOpts.get_shifted_grad_op(pressure_->name());
  if ( skewSymmetric || shiftedGradOpU || shiftedGradOpP )
    throw std::runtime_error("MomentumContinuityElemKernel::error: skewSymmetric || shiftedGradOpU || shiftedGradOpP");
}
  
template<class AlgTraits>
MomentumContinuityElemKernel<AlgTraits>::~MomentumContinuityElemKernel()
{}

template<typename AlgTraits>
void
MomentumContinuityElemKernel<AlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
  dt_ = timeIntegrator.get_time_step();
  gamma1_ = timeIntegrator.get_gamma1();
  gamma2_ = timeIntegrator.get_gamma2();
  gamma3_ = timeIntegrator.get_gamma3();
  projTimeScale_ = dt_/gamma1_;
}
  
template<class AlgTraits>
void
MomentumContinuityElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>& /*lhs*/,
  SharedMemView<DoubleType *>& /*rhs*/,
  ScratchViews<DoubleType>& /*scratchViews*/)
{
  throw std::runtime_error("MomentumContinuityElemKernel::execute() Error: Not implemented!");
}

INSTANTIATE_KERNEL(MomentumContinuityElemKernel);

}  // nalu
}  // sierra
