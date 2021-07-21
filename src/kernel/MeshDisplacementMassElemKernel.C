/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/MeshDisplacementMassElemKernel.h"
#include "AlgTraits.h"
#include "master_element/MasterElement.h"
#include "TimeIntegrator.h"
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
MeshDisplacementMassElemKernel<AlgTraits>::MeshDisplacementMassElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  VectorFieldType* meshDisplacement,
  ElemDataRequests& dataPreReqs,
  const bool lumpedMass)
  : Kernel(),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap()),
    lumpedMass_(lumpedMass)
{
  // save off fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();

  meshDisplacement_ = &(meshDisplacement->field_of_state(stk::mesh::StateN));
  meshDisplacementNp1_ = &(meshDisplacement->field_of_state(stk::mesh::StateNP1));

  const int numStates = meshDisplacement->number_of_states(); 
  meshDisplacementNm1_ = (numStates == 2) ? meshDisplacement_ : &(meshDisplacement->field_of_state(stk::mesh::StateNM1));

  coordinates_ = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  density_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");

  MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);

  // compute shape function
  if ( lumpedMass_ )
    get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shifted_shape_fcn(ptr);}, v_shape_function_);
  else
    get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shape_fcn(ptr);}, v_shape_function_);

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*meshDisplacementNp1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*meshDisplacementNm1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*meshDisplacement_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*density_, 1);

  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
}

template<typename AlgTraits>
MeshDisplacementMassElemKernel<AlgTraits>::~MeshDisplacementMassElemKernel()
{}

template<typename AlgTraits>
void
MeshDisplacementMassElemKernel<AlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
  const double dt_ = timeIntegrator.get_time_step();
  const double dt_nm1_ = timeIntegrator.get_time_step(NALU_STATE_NM1);
  const double dt_avg_ = 0.5*dt_+dt_nm1_;
  gamma_1_ = 1.0/(dt_*dt_avg_);
  gamma_2_ = -1.0/(dt_*dt_avg_)-1.0/(dt_nm1_*dt_avg_);
  gamma_3_ = 1.0/(dt_nm1_*dt_avg_);

}

template<typename AlgTraits>
void
MeshDisplacementMassElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>& lhs,
  SharedMemView<DoubleType *>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  NALU_ALIGNED DoubleType w_uIp[AlgTraits::nDim_];
  NALU_ALIGNED DoubleType meshDispIp[AlgTraits::nDim_];
  NALU_ALIGNED DoubleType meshDispNp1Ip[AlgTraits::nDim_];
  NALU_ALIGNED DoubleType meshDispNm1Ip[AlgTraits::nDim_];

  SharedMemView<DoubleType**>& meshDisp = scratchViews.get_scratch_view_2D(*meshDisplacement_);
  SharedMemView<DoubleType**>& meshDispNp1 = scratchViews.get_scratch_view_2D(*meshDisplacementNp1_);
  SharedMemView<DoubleType**>& meshDispNm1 = scratchViews.get_scratch_view_2D(*meshDisplacementNm1_);
  SharedMemView<DoubleType*>& rho = scratchViews.get_scratch_view_1D(*density_);

  SharedMemView<DoubleType*>& v_scv_volume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;


  for ( int ip = 0; ip < AlgTraits::numScvIp_; ++ip ) {

    // nearest node to ip
    const int nearestNode = ipNodeMap_[ip];
    const int nNodeNdim = nearestNode*AlgTraits::nDim_;

    // zero out vector
    for (int j = 0; j < AlgTraits::nDim_; ++j) {
      meshDispIp[j] = 0.0;
      meshDispNp1Ip[j] = 0.0;
      meshDispNm1Ip[j] = 0.0;
      w_uIp[j] = 0.0;
    }

    // zero out scalar
    DoubleType densityIp = 0.0; 

    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      // save off shape function
      const DoubleType r = v_shape_function_(ip,ic);
      for (int j = 0; j < AlgTraits::nDim_; ++j) {
        meshDispIp[j] += meshDisp(ic,j)*r;
        meshDispNp1Ip[j] += meshDispNp1(ic,j)*r;
        meshDispNm1Ip[j] += meshDispNm1(ic,j)*r;
      }
      densityIp += rho(ic)*r; 
    }

    const DoubleType scV = v_scv_volume(ip);
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const DoubleType lhsMass = densityIp*v_shape_function_(ip,ic)*gamma_1_;

      const int icNdim = ic*AlgTraits::nDim_;
      for (int j = 0; j < AlgTraits::nDim_; ++j) {
        lhs(nNodeNdim+j,icNdim+j) += (lhsMass)*scV;
      }
    }

    for (int j = 0; j < AlgTraits::nDim_; ++j) {
      // Estimate 2nd derivative
      const DoubleType mass = densityIp*(gamma_1_*meshDispNp1Ip[j] + gamma_2_*meshDispIp[j] + gamma_3_*meshDispNm1Ip[j]);
      rhs(nNodeNdim+j) -= mass*scV;
    }

  }
}

INSTANTIATE_KERNEL(MeshDisplacementMassElemKernel);

}  // nalu
}  // sierra
