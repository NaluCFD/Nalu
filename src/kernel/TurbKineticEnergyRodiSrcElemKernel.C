/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/TurbKineticEnergyRodiSrcElemKernel.h"
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
TurbKineticEnergyRodiSrcElemKernel<AlgTraits>::TurbKineticEnergyRodiSrcElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    beta_(solnOpts.thermalExpansionCoeff_),
    turbPr_(solnOpts.get_turb_prandtl("enthalpy")),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap())
{
  // save off fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  coordinates_ = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  dhdx_ = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "dhdx");
  specificHeat_ = metaData.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "specific_heat");
  tvisc_ = metaData.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "turbulent_viscosity");

  const std::vector<double>& solnOptsGravity = solnOpts.get_gravity_vector(AlgTraits::nDim_);
  for (int i = 0; i < AlgTraits::nDim_; i++)
    gravity_(i) = solnOptsGravity[i];
  
  MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // required fields
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*dhdx_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*specificHeat_, 1);
  dataPreReqs.add_gathered_nodal_field(*tvisc_, 1);
  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
}

template<typename AlgTraits>
TurbKineticEnergyRodiSrcElemKernel<AlgTraits>::~TurbKineticEnergyRodiSrcElemKernel()
{}

template<typename AlgTraits>
void
TurbKineticEnergyRodiSrcElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>&lhs,
  SharedMemView<DoubleType *>&rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  SharedMemView<DoubleType**>& v_dhdx = scratchViews.get_scratch_view_2D(*dhdx_);
  SharedMemView<DoubleType*>& v_specificHeat = scratchViews.get_scratch_view_1D(
    *specificHeat_);
  SharedMemView<DoubleType*>& v_tvisc = scratchViews.get_scratch_view_1D(
    *tvisc_);
  SharedMemView<DoubleType*>& v_scv_volume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;

  for (int ip=0; ip < AlgTraits::numScvIp_; ++ip) {
    const int nearestNode = ipNodeMap_[ip];

    DoubleType gidhdxi = 0.0;
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      gidhdxi += gravity_[i]*v_dhdx(nearestNode,i);
    }

    // rhs assembly, all lumped; no implicit contribution
    const DoubleType scvol = v_scv_volume(ip);
    rhs(nearestNode) += 
      beta_*v_tvisc(nearestNode)/turbPr_*gidhdxi/v_specificHeat(nearestNode)*scvol;
  }
}

INSTANTIATE_KERNEL(TurbKineticEnergyRodiSrcElemKernel);

}  // nalu
}  // sierra
