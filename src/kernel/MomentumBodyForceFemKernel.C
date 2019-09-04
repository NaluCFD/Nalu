/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/MomentumBodyForceFemKernel.h"
#include "AlgTraits.h"
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

template<typename AlgTraits>
MomentumBodyForceFemKernel<AlgTraits>::MomentumBodyForceFemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ElemDataRequests& dataPreReqs)
  : Kernel()
{
  // Save of required fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  coordinates_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  
  // extract master element
  MasterElement *meFEM = sierra::nalu::MasterElementRepo::get_fem_master_element(AlgTraits::topo_);
  
  // copy ip weights into our 1-d view
  for ( int k = 0; k < AlgTraits::numGp_; ++k )
    v_ip_weight_[k] = meFEM->weights_[k];
  
  // extract the body force
  const std::map<std::string, std::vector<double> >::const_iterator iparams
    = solnOpts.elemSrcTermParamMap_.find("momentum");
  if ( iparams != solnOpts.elemSrcTermParamMap_.end() ) {
    std::vector<double> bf = (*iparams).second;
    // copy constant source into our 1-d view
    for ( int k = 0; k < AlgTraits::nDim_; ++k ) {
      v_body_force_[k] = bf[k];
    }
  }
  else {
    throw std::runtime_error("MomentumBodyForceFemKernel is missing body force specification");
  }
  
  // master element, shape function is shifted consistently
  if ( solnOpts.get_shifted_grad_op("velocity") )
    get_fem_shape_fn_data<AlgTraits>([&](double* ptr){meFEM->shifted_shape_fcn(ptr);}, v_shape_function_);
  else
    get_fem_shape_fn_data<AlgTraits>([&](double* ptr){meFEM->shape_fcn(ptr);}, v_shape_function_);
  
  // add FEM master element
  dataPreReqs.add_fem_volume_me(meFEM);
  dataPreReqs.add_master_element_call(FEM_DET_J, CURRENT_COORDINATES);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);  
}

template<typename AlgTraits>
MomentumBodyForceFemKernel<AlgTraits>::~MomentumBodyForceFemKernel()
{
  // does nothing
}

template<typename AlgTraits>
void
MomentumBodyForceFemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType**>& lhs,
  SharedMemView<DoubleType*>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  SharedMemView<DoubleType*>& v_det_j = scratchViews.get_me_views(CURRENT_COORDINATES).det_j_fem;

  for ( int ip = 0; ip < AlgTraits::numGp_; ++ip ) {
  
    // start the assembly (collect ip scalings)
    const DoubleType ipFactor = v_det_j(ip)*v_ip_weight_(ip);

    // row ir
    for ( int ir = 0; ir < AlgTraits::nodesPerElement_; ++ir) {
      
      const int irNdim = ir*AlgTraits::nDim_;
      const DoubleType wIr = v_shape_function_(ip,ir);
                      
      // component i
      for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
        rhs(irNdim+i) += wIr*v_body_force_[i]*ipFactor; 
      }
    }
  }
}
  
INSTANTIATE_FEM_KERNEL(MomentumBodyForceFemKernel);

}  // nalu
}  // sierra
