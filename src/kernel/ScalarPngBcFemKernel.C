/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/ScalarPngBcFemKernel.h"
#include "master_element/MasterElement.h"
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

template<class BcAlgTraits>
ScalarPngBcFemKernel<BcAlgTraits>::ScalarPngBcFemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  const std::string fieldName,
  ElemDataRequests& dataPreReqs)
  : Kernel()
{
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  scalarQ_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, fieldName);
  VectorFieldType *coordinates = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
 
  // extract master element
  MasterElement *meFEM = sierra::nalu::MasterElementRepo::get_fem_master_element(BcAlgTraits::topo_);

  // copy ip weights into our 1-d view
  for ( int k = 0; k < BcAlgTraits::numFaceIp_; ++k )
    v_ip_weight_[k] = meFEM->weights_[k];
  
  // compute and save shape function
  if ( solnOpts.get_shifted_grad_op(fieldName) )
    get_face_shape_fn_data<BcAlgTraits>([&](double* ptr){meFEM->shifted_shape_fcn(ptr);}, vf_shape_function_);
  else
    get_face_shape_fn_data<BcAlgTraits>([&](double* ptr){meFEM->shape_fcn(ptr);}, vf_shape_function_);
  
  // add master elements
  dataPreReqs.add_fem_volume_me(meFEM);
  
  // fields and data; mdot not gathered as element data
  dataPreReqs.add_coordinates_field(*coordinates, BcAlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*scalarQ_, 1);
  dataPreReqs.add_master_element_call(FEM_DET_J, CURRENT_COORDINATES);  
  dataPreReqs.add_master_element_call(FEM_NORMAL, CURRENT_COORDINATES);  
}

template<class BcAlgTraits>
ScalarPngBcFemKernel<BcAlgTraits>::~ScalarPngBcFemKernel()
{}

template<class BcAlgTraits>
void
ScalarPngBcFemKernel<BcAlgTraits>::execute(
  SharedMemView<DoubleType **>& lhs,
  SharedMemView<DoubleType *>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  SharedMemView<DoubleType*>& v_scalarQ = scratchViews.get_scratch_view_1D(*scalarQ_);
  SharedMemView<DoubleType*>& v_det_j = scratchViews.get_me_views(CURRENT_COORDINATES).det_j_fem;
  SharedMemView<DoubleType**>& v_normal = scratchViews.get_me_views(CURRENT_COORDINATES).normal_fem;
  
  for ( int ip = 0; ip < BcAlgTraits::numFaceIp_; ++ip ) {
    
    // interpolate to bip
    DoubleType qBip = 0.0;
    for ( int ic = 0; ic < BcAlgTraits::nodesPerFace_; ++ic ) {
      const DoubleType r = vf_shape_function_(ip,ic);
      qBip += r*v_scalarQ(ic);
    }
    
    // start the assembly (collect ip scalings)
    const DoubleType ipFactor = v_det_j(ip)*v_ip_weight_(ip);
    
    // row ir
    for ( int ir = 0; ir < BcAlgTraits::nodesPerElement_; ++ir) {
      
      const DoubleType wIr = vf_shape_function_(ip,ir);
      
      // component j
      for ( int j = 0; j < BcAlgTraits::nDim_; ++j ) {
        rhs(ir*BcAlgTraits::nDim_+j) += wIr*qBip*ipFactor*v_normal(ip,j);
      }
    }
  }
}

INSTANTIATE_FEM_KERNEL_FACE(ScalarPngBcFemKernel);

}  // nalu
}  // sierra
