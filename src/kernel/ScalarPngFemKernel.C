/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/ScalarPngFemKernel.h"
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
ScalarPngFemKernel<AlgTraits>::ScalarPngFemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  std::string independentDofName,
  std::string dofName,
  ElemDataRequests& dataPreReqs)
  : Kernel()
{
  // Save of required fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  scalarQ_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, independentDofName);
  Gjq_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, dofName);

  coordinates_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  // extract master element
  MasterElement *meFEM = sierra::nalu::MasterElementRepo::get_fem_master_element(AlgTraits::topo_);
  
  // copy ip weights into our 1-d view
  for ( int k = 0; k < AlgTraits::numGp_; ++k )
    v_ip_weight_[k] = meFEM->weights_[k];

  // master element, shape function is shifted consistently
  if ( solnOpts.get_shifted_grad_op(dofName) )
    get_fem_shape_fn_data<AlgTraits>([&](double* ptr){meFEM->shifted_shape_fcn(ptr);}, v_shape_function_);
  else
    get_fem_shape_fn_data<AlgTraits>([&](double* ptr){meFEM->shape_fcn(ptr);}, v_shape_function_);
  
  // add FEM master element
  dataPreReqs.add_fem_volume_me(meFEM);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*scalarQ_, 1);
  dataPreReqs.add_gathered_nodal_field(*Gjq_, AlgTraits::nDim_);
  
  if ( solnOpts.get_shifted_grad_op(dofName) )
    dataPreReqs.add_master_element_call(FEM_SHIFTED_GRAD_OP, CURRENT_COORDINATES);
  else
    dataPreReqs.add_master_element_call(FEM_GRAD_OP, CURRENT_COORDINATES);
 
}

template<typename AlgTraits>
ScalarPngFemKernel<AlgTraits>::~ScalarPngFemKernel()
{
  // does nothing
}

template<typename AlgTraits>
void
ScalarPngFemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType**>& lhs,
  SharedMemView<DoubleType*>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  NALU_ALIGNED DoubleType w_GjqIp[AlgTraits::nDim_];

  SharedMemView<DoubleType*>& v_scalarQ = scratchViews.get_scratch_view_1D(*scalarQ_);
  SharedMemView<DoubleType**>& v_Gjq = scratchViews.get_scratch_view_2D(*Gjq_);
  
  SharedMemView<DoubleType***>& v_dndx = scratchViews.get_me_views(CURRENT_COORDINATES).dndx_fem;
  SharedMemView<DoubleType*>& v_det_j = scratchViews.get_me_views(CURRENT_COORDINATES).det_j_fem;

  for ( int ip = 0; ip < AlgTraits::numGp_; ++ip ) {
    
    // compute ip values
    DoubleType qIp = 0.0;
    for ( int j = 0; j < AlgTraits::nDim_; ++j )
      w_GjqIp[j] = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const DoubleType r = v_shape_function_(ip,ic);
      qIp += r*v_scalarQ(ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        w_GjqIp[j] += r*v_Gjq(ic,j);
      }
    }

    // start the assembly (collect ip scalings)
    const DoubleType ipFactor = v_det_j(ip)*v_ip_weight_(ip);

    // row ir
    for ( int ir = 0; ir < AlgTraits::nodesPerElement_; ++ir) {
      
      const DoubleType wIr = v_shape_function_(ip,ir);
      
      // component j
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        // column ic
        for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
          lhs(ir*AlgTraits::nDim_+j,ic*AlgTraits::nDim_+j) += wIr*v_shape_function_(ip,ic)*ipFactor;
        }
        rhs(ir*AlgTraits::nDim_+j) -= (wIr*w_GjqIp[j] + v_dndx(ip,ir,j)*qIp)*ipFactor;
      }
    }
  }
}
  
INSTANTIATE_FEM_KERNEL(ScalarPngFemKernel);

}  // nalu
}  // sierra
