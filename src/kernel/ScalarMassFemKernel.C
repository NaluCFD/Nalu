/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/ScalarMassFemKernel.h"
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
ScalarMassFemKernel<AlgTraits>::ScalarMassFemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ScalarFieldType* scalarQ,
  ScalarFieldType* density,
  ElemDataRequests& dataPreReqs)
  : Kernel()
{
  // Save of required fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  coordinates_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  scalarQNp1_ = &(scalarQ->field_of_state(stk::mesh::StateNP1));
  scalarQN_ = &(scalarQ->field_of_state(stk::mesh::StateN));
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  densityN_ = &(density->field_of_state(stk::mesh::StateN));
  if (scalarQ->number_of_states() == 2) {
    scalarQNm1_ = scalarQN_;
    densityNm1_ = densityN_;
  }
  else {
    scalarQNm1_ = &(scalarQ->field_of_state(stk::mesh::StateNM1));
    densityNm1_ = &(density->field_of_state(stk::mesh::StateNM1));
  }
  
  // extract master element
  MasterElement *meFEM = sierra::nalu::MasterElementRepo::get_fem_master_element(AlgTraits::topo_);
  
  // copy ip weights into our 1-d view
  for ( int k = 0; k < AlgTraits::numGp_; ++k )
    v_ip_weight_[k] = meFEM->weights_[k];
  
  // master element, shape function is shifted consistently
  if ( solnOpts.get_shifted_grad_op(scalarQ->name()) )
    get_fem_shape_fn_data<AlgTraits>([&](double* ptr){meFEM->shifted_shape_fcn(ptr);}, v_shape_function_);
  else
    get_fem_shape_fn_data<AlgTraits>([&](double* ptr){meFEM->shape_fcn(ptr);}, v_shape_function_);
  
  // add FEM master element
  dataPreReqs.add_fem_volume_me(meFEM);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*scalarQNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*scalarQN_, 1);
  dataPreReqs.add_gathered_nodal_field(*scalarQNm1_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityN_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityNm1_, 1);
  
  dataPreReqs.add_master_element_call(FEM_DET_J, CURRENT_COORDINATES);
}

template<typename AlgTraits>
ScalarMassFemKernel<AlgTraits>::~ScalarMassFemKernel()
{
  // does nothing
}

template<typename AlgTraits>
void
ScalarMassFemKernel<AlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
  dt_ = timeIntegrator.get_time_step();
  gamma1_ = timeIntegrator.get_gamma1();
  gamma2_ = timeIntegrator.get_gamma2();
  gamma3_ = timeIntegrator.get_gamma3(); // gamma3 may be zero
}

template<typename AlgTraits>
void
ScalarMassFemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType**>& lhs,
  SharedMemView<DoubleType*>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  SharedMemView<DoubleType*>& v_scalarQNp1 = scratchViews.get_scratch_view_1D(*scalarQNp1_);
  SharedMemView<DoubleType*>& v_scalarQN = scratchViews.get_scratch_view_1D(*scalarQN_);
  SharedMemView<DoubleType*>& v_scalarQNm1 = scratchViews.get_scratch_view_1D(*scalarQNm1_);
  SharedMemView<DoubleType*>& v_densityNp1 = scratchViews.get_scratch_view_1D(*densityNp1_);
  SharedMemView<DoubleType*>& v_densityN = scratchViews.get_scratch_view_1D(*densityN_);
  SharedMemView<DoubleType*>& v_densityNm1 = scratchViews.get_scratch_view_1D(*densityNm1_);
  
  SharedMemView<DoubleType*>& v_det_j = scratchViews.get_me_views(CURRENT_COORDINATES).det_j_fem;

  for ( int ip = 0; ip < AlgTraits::numGp_; ++ip ) {
    // compute ip property
    DoubleType qNp1Ip = 0.0;
    DoubleType qNIp = 0.0;
    DoubleType qNm1Ip = 0.0;
    DoubleType rhoNp1Ip = 0.0;
    DoubleType rhoNIp= 0.0;
    DoubleType rhoNm1Ip = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const DoubleType r = v_shape_function_(ip,ic);
      qNp1Ip += r*v_scalarQNp1(ic);
      qNIp += r*v_scalarQN(ic);
      qNm1Ip += r*v_scalarQNm1(ic);
      rhoNp1Ip += r*v_densityNp1(ic);
      rhoNIp += r*v_densityN(ic);
      rhoNm1Ip += r*v_densityNm1(ic);
    }

    // start the assembly (collect ip scalings)
    const DoubleType ipFactor = v_det_j(ip)*v_ip_weight_(ip);

    // row ir
    for ( int ir = 0; ir < AlgTraits::nodesPerElement_; ++ir) {
      
      const DoubleType wIr = v_shape_function_(ip,ir);
      
      // column ic
      for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
        lhs(ir,ic) += gamma1_*rhoNp1Ip*wIr*v_shape_function_(ip,ic)/dt_*ipFactor;
      }
      rhs(ir) -= wIr*(gamma1_*qNp1Ip*rhoNp1Ip + gamma2_*qNIp*rhoNIp + gamma3_*qNm1Ip*rhoNm1Ip)/dt_*ipFactor;
    }
  }
}
  
INSTANTIATE_FEM_KERNEL(ScalarMassFemKernel);

}  // nalu
}  // sierra
