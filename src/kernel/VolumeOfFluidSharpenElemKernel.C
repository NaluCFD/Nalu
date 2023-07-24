/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/VolumeOfFluidSharpenElemKernel.h"
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
VolumeOfFluidSharpenElemKernel<AlgTraits>::VolumeOfFluidSharpenElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ScalarFieldType* vof,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    cAlpha_(solnOpts.vofCalpha_),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap())
{
  // save off fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  
  vofNp1_ = &(vof->field_of_state(stk::mesh::StateNP1));
  interfaceNormal_ = metaData.get_field<double>(stk::topology::NODE_RANK, "interface_normal");
  std::string velocity_name = solnOpts.does_mesh_move() ? "velocity_rtm" : "velocity";
  velocityRTM_ = metaData.get_field<double>(stk::topology::NODE_RANK, velocity_name);
  coordinates_ = metaData.get_field<double>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  
  MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*velocityRTM_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*interfaceNormal_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*vofNp1_, 1);
  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(SCV_GRAD_OP, CURRENT_COORDINATES);
}

template<typename AlgTraits>
VolumeOfFluidSharpenElemKernel<AlgTraits>::~VolumeOfFluidSharpenElemKernel()
{}

template<typename AlgTraits>
void
VolumeOfFluidSharpenElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>& lhs,
  SharedMemView<DoubleType *>&rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  NALU_ALIGNED DoubleType w_magU[AlgTraits::nodesPerElement_];

  SharedMemView<DoubleType*>& v_vofNp1 = scratchViews.get_scratch_view_1D(*vofNp1_);
  SharedMemView<DoubleType**>& v_vrtm = scratchViews.get_scratch_view_2D(*velocityRTM_);
  SharedMemView<DoubleType**>& v_interfaceNormal = scratchViews.get_scratch_view_2D(*interfaceNormal_);
  SharedMemView<DoubleType*>& v_scVolume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;
  SharedMemView<DoubleType***>& v_dndx = scratchViews.get_me_views(CURRENT_COORDINATES).dndx_scv;

  // determine nodal velocity magnitude
  for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
    DoubleType magU = 0.0;
    for ( int k = 0; k < AlgTraits::nDim_; ++k ) {
      magU += v_vrtm(ic,k)*v_vrtm(ic,k);
    }
    w_magU[ic] = stk::math::sqrt(magU);
  }

  for ( int ip = 0; ip < AlgTraits::numScvIp_; ++ip ) {

    // nearest node to ip
    const int nearestNode = ipNodeMap_[ip];
    const DoubleType scV = v_scVolume(ip);
    
    DoubleType sharpen = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      DoubleType sum = 0.0;
      DoubleType lhsOne = 0.0;
      DoubleType lhsTwo = 0.0;
      const DoubleType vofIc = v_vofNp1(ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        const DoubleType ucj = cAlpha_*w_magU[ic]*v_interfaceNormal(ic,j);
        lhsOne += ucj*v_dndx(ip,ic,j);
        lhsTwo -= 2.0*ucj*vofIc*v_dndx(ip,ic,j);
        sum += ucj*v_dndx(ip,ic,j);
      }
      sharpen += sum*(vofIc*(1.0-vofIc));
      lhs(nearestNode,ic) += (lhsOne*1.0+lhsTwo*0.0)*scV;
    }
    
    // assemble rhs
    rhs(nearestNode) -= sharpen*scV;
  }
}

INSTANTIATE_KERNEL(VolumeOfFluidSharpenElemKernel);

}  // nalu
}  // sierra
