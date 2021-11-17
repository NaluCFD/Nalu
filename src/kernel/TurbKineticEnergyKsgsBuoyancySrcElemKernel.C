/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/TurbKineticEnergyKsgsBuoyancySrcElemKernel.h"
#include "AlgTraits.h"
#include "Enums.h"
#include "SolutionOptions.h"
#include "master_element/MasterElement.h"

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
TurbKineticEnergyKsgsBuoyancySrcElemKernel<AlgTraits>::TurbKineticEnergyKsgsBuoyancySrcElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    CbTwo_(solnOpts.get_turb_model_constant(TM_CbTwo)),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap())
{
  // save off fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  coordinates_ = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  ScalarFieldType *tke = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_ke");
  tkeNp1_ = &(tke->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType *density = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNp1_ = &(density->field_of_state(stk::mesh::StateNP1));
  dualNodalVolume_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);
  get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shape_fcn(ptr);}, v_shape_function_);

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // required fields
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*tkeNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*dualNodalVolume_, 1);
  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(SCV_GRAD_OP, CURRENT_COORDINATES);

  // extract gravity from solution options
  std::array<double, 3> soGravity = solnOpts.gravity_;
  
  // error checks
  if ( AlgTraits::nDim_ != 3 || soGravity.size() != 3 ) 
    throw std::runtime_error("TurbKineticEnergyKsgsBuoyancySrcElemKernel is intended to be 3D only (check gravity)");
  
  for ( size_t k = 0; k < soGravity.size(); ++k )
    gravity_[k] = soGravity[k];
}
  
template<typename AlgTraits>
TurbKineticEnergyKsgsBuoyancySrcElemKernel<AlgTraits>::~TurbKineticEnergyKsgsBuoyancySrcElemKernel()
{}

template<typename AlgTraits>
void
TurbKineticEnergyKsgsBuoyancySrcElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>&lhs,
  SharedMemView<DoubleType *>&rhs,
  ScratchViews<DoubleType>& scratchViews)
{ 
  SharedMemView<DoubleType*>& v_tkeNp1 = scratchViews.get_scratch_view_1D(
    *tkeNp1_);
  SharedMemView<DoubleType*>& v_densityNp1 = scratchViews.get_scratch_view_1D(
    *densityNp1_);
  SharedMemView<DoubleType*>& v_dualNodalVolume = scratchViews.get_scratch_view_1D(
    *dualNodalVolume_);
  SharedMemView<DoubleType*>& v_scv_volume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;
  SharedMemView<DoubleType***>& v_dndx = scratchViews.get_me_views(CURRENT_COORDINATES).dndx_scv;

  for (int ip=0; ip < AlgTraits::numScvIp_; ++ip) {
    
    const int nearestNode = ipNodeMap_[ip];
        
    // zero
    DoubleType tkeIp = 0.0;
    DoubleType dualNodalVolIp = 0.0;
    for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
      w_drhodx_[j] = 0.0;
    }

    // compute integration point values, including spatial derivative
    for (int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic) {      
      const DoubleType r = v_shape_function_(ip, ic);      
      tkeIp += r*v_tkeNp1(ic);
      dualNodalVolIp += r*v_dualNodalVolume(ic);
      const DoubleType rhoIc = v_densityNp1[ic];
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        w_drhodx_[j] += v_dndx(ip,ic,j)*rhoIc;
      }
    }
    
    // dot product; gradRho(dot)g
    DoubleType dotProduct = 0.0;
    for ( int i = 0; i < AlgTraits::nDim_; ++i )
      dotProduct += w_drhodx_[i]*gravity_[i];

    // filter size
    const DoubleType filter = stk::math::cbrt(dualNodalVolIp);

    // cross product magnitude; | gradRho x g |
    const DoubleType crossProductMag = cross_product_magnitude();
    
    // assemble rhs; no chance for stable LHS
    rhs[nearestNode] += CbTwo_*filter*stk::math::sqrt(tkeIp)*(crossProductMag - dotProduct)*v_scv_volume(ip);
  }
} 

//--------------------------------------------------------------------------
//-------- cross_product_magnitude -----------------------------------------
//--------------------------------------------------------------------------
template<class AlgTraits>
DoubleType
TurbKineticEnergyKsgsBuoyancySrcElemKernel<AlgTraits>::cross_product_magnitude()
{
  const DoubleType crossX =   w_drhodx_[1]*gravity_[2] - w_drhodx_[2]*gravity_[1];
  const DoubleType crossY = -(w_drhodx_[0]*gravity_[2] - w_drhodx_[2]*gravity_[0]);
  const DoubleType crossZ =   w_drhodx_[0]*gravity_[1] - w_drhodx_[1]*gravity_[0];
  return stk::math::sqrt(crossX*crossX + crossY*crossY + crossZ*crossZ);
}  
  
INSTANTIATE_KERNEL(TurbKineticEnergyKsgsBuoyancySrcElemKernel);

}  // nalu
}  // sierra
