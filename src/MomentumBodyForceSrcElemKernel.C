/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "MomentumBodyForceSrcElemKernel.h"
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

template<class AlgTraits>
MomentumBodyForceSrcElemKernel<AlgTraits>::MomentumBodyForceSrcElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    ipNodeMap_(sierra::nalu::get_volume_master_element(AlgTraits::topo_)->ipNodeMap())
{
  auto isrc = solnOpts.srcTermParamMap_.find("body_force");
  std::vector<double> theParams = (isrc != solnOpts.srcTermParamMap_.end()) ? isrc->second : std::vector<double>();
  if (theParams.empty()) {
    std::string msg = "Body force vector not set.  Must specify source_term_parameters for body force in input file";
    throw std::runtime_error(msg);
  }

  ThrowRequireMsg(theParams.size() == AlgTraits::nDim_, "Body force vector has an incorrect number of components");
  for (int j = 0; j < AlgTraits::nDim_; ++j) {
    v_body_force_(j) = theParams[j];
  }

  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  coordinates_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  MasterElement* meSCV = sierra::nalu::get_volume_master_element(AlgTraits::topo_);

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
}

template<class AlgTraits>
MomentumBodyForceSrcElemKernel<AlgTraits>::~MomentumBodyForceSrcElemKernel()
{}

template<typename AlgTraits>
void
MomentumBodyForceSrcElemKernel<AlgTraits>::execute(
  SharedMemView<double**>& /* lhs */,
  SharedMemView<double*>& rhs,
  stk::mesh::Entity /*element*/,
  ScratchViews& scratchViews)
{
  SharedMemView<double*>& v_scv_volume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;
  for (int ip=0; ip < AlgTraits::numScvIp_; ++ip) {
    const int nearestNode = ipNodeMap_[ip];

    const double scV = v_scv_volume(ip);
    const int nnNdim = nearestNode * AlgTraits::nDim_;
    for (int j=0; j < AlgTraits::nDim_; j++) {
      rhs(nnNdim + j) += scV * v_body_force_(j);
    }
  }
}
INSTANTIATE_KERNEL(MomentumBodyForceSrcElemKernel);

}  // nalu
}  // sierra
