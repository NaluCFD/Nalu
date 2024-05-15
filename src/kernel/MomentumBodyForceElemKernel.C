/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/MomentumBodyForceElemKernel.h"
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
MomentumBodyForceElemKernel<AlgTraits>::MomentumBodyForceElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap())
{
  // Save of required fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  coordinates_ = metaData.get_field<double>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  // extract master element
  MasterElement* meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);
  
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
    throw std::runtime_error("MomentumBodyForceElemKernel is missing body force specification");
  }
    
  // add master element
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
}

template<typename AlgTraits>
MomentumBodyForceElemKernel<AlgTraits>::~MomentumBodyForceElemKernel()
{
  // does nothing
}

template<typename AlgTraits>
void
MomentumBodyForceElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType**>& /*lhs*/,
  SharedMemView<DoubleType*>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  SharedMemView<DoubleType*>& v_scv_volume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;

  for (int ip=0; ip < AlgTraits::numScvIp_; ++ip) {

    const int nearestNode = ipNodeMap_[ip];

    // Compute RHS
    const DoubleType scV = v_scv_volume(ip);
    const int nnNdim = nearestNode*AlgTraits::nDim_;
    for (int i=0; i < AlgTraits::nDim_; ++i) {
      rhs(nnNdim+i) += scV*v_body_force_[i]; 
    }
  }
}
  
INSTANTIATE_KERNEL(MomentumBodyForceElemKernel);

}  // nalu
}  // sierra
