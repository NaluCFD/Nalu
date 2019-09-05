/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/ContinuityInflowFemKernel.h"
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
ContinuityInflowFemKernel<BcAlgTraits>::ContinuityInflowFemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions &solnOpts,
  const bool &useShifted,
  ElemDataRequests &dataPreReqs)
  : Kernel(),
    useShifted_(useShifted),
    projTimeScale_(1.0)
{
  // save off fields
  const stk::mesh::MetaData &metaData = bulkData.mesh_meta_data();
  velocityBC_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "cont_velocity_bc");
  densityBC_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  
  // extract field not required in execute()
  VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  
  // extract master element
  MasterElement *meFC = sierra::nalu::MasterElementRepo::get_fem_master_element(BcAlgTraits::topo_);

  // copy ip weights into our 1-d view
  for ( int k = 0; k < BcAlgTraits::numFaceIp_; ++k )
    v_ip_weight_[k] = meFC->weights_[k];

  // add master elements
  dataPreReqs.add_fem_volume_me(meFC);
 
  // required fields
  dataPreReqs.add_gathered_nodal_field(*velocityBC_, BcAlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*densityBC_, 1);
  dataPreReqs.add_coordinates_field(*coordinates, BcAlgTraits::nDim_, CURRENT_COORDINATES);
  
  // manage master element requirements
  dataPreReqs.add_master_element_call(FEM_DET_J, CURRENT_COORDINATES);  
  dataPreReqs.add_master_element_call(FEM_NORMAL, CURRENT_COORDINATES);  

  if ( useShifted )
    get_face_shape_fn_data<BcAlgTraits>([&](double* ptr){meFC->shifted_shape_fcn(ptr);}, v_shape_function_);
  else
    get_face_shape_fn_data<BcAlgTraits>([&](double* ptr){meFC->shape_fcn(ptr);}, v_shape_function_);
 }
  
template<typename BcAlgTraits>
ContinuityInflowFemKernel<BcAlgTraits>::~ContinuityInflowFemKernel()
{}

template<typename BcAlgTraits>
void
ContinuityInflowFemKernel<BcAlgTraits>::setup(const TimeIntegrator &timeIntegrator)
{
  const double dt = timeIntegrator.get_time_step();
  const double gamma1 = timeIntegrator.get_gamma1();
  projTimeScale_ = dt/gamma1;
}

template<typename BcAlgTraits>
void
ContinuityInflowFemKernel<BcAlgTraits>::execute(
  SharedMemView<DoubleType **>&/*lhs*/,
  SharedMemView<DoubleType *>&rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  NALU_ALIGNED DoubleType w_rho_uBip[BcAlgTraits::nDim_];
  
  SharedMemView<DoubleType**>& v_velocityBC = scratchViews.get_scratch_view_2D(*velocityBC_);
  SharedMemView<DoubleType*>& v_density = scratchViews.get_scratch_view_1D(*densityBC_);

  SharedMemView<DoubleType*>& v_det_j = scratchViews.get_me_views(CURRENT_COORDINATES).det_j_fem;
  SharedMemView<DoubleType**>& v_normal = scratchViews.get_me_views(CURRENT_COORDINATES).normal_fem;

  for (int ip = 0; ip < BcAlgTraits::numFaceIp_; ++ip) {
     
    // zero out vector quantities
    for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
      w_rho_uBip[j] = 0.0;
    }
    
    for ( int ic = 0; ic < BcAlgTraits::nodesPerFace_; ++ic ) {
      const DoubleType r = v_shape_function_(ip,ic);
      const DoubleType rhoIC = v_density(ic);
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        w_rho_uBip[j] += r*rhoIC*v_velocityBC(ic,j);
      }
    }
    
    DoubleType massFlux = 0.0;
    for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
      massFlux += w_rho_uBip[j]*v_normal(ip,j);
    }
    
    // start the assembly (collect ip scalings)
    const DoubleType ipFactor = v_det_j(ip)*v_ip_weight_(ip);
 
    // row ir
    for ( int ir = 0; ir < BcAlgTraits::nodesPerElement_; ++ir) {  
      const DoubleType wIr = v_shape_function_(ip,ir);
      rhs(ir) -= wIr*massFlux/projTimeScale_*ipFactor;
    } 
  }
}

INSTANTIATE_FEM_KERNEL_FACE(ContinuityInflowFemKernel);

}  // nalu
}  // sierra
