/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "TurbKineticEnergyKsgsSrcElemKernel.h"
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
TurbKineticEnergyKsgsSrcElemKernel<AlgTraits>::TurbKineticEnergyKsgsSrcElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    cEps_(solnOpts.get_turb_model_constant(TM_cEps)),
    tkeProdLimitRatio_(solnOpts.get_turb_model_constant(TM_tkeProdLimitRatio)),
    ipNodeMap_(sierra::nalu::get_volume_master_element(AlgTraits::topo_)->ipNodeMap())
{
  // save off fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  tkeNp1_ = metaData.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "turbulent_ke");
  densityNp1_ = metaData.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "density");
  tvisc_ = metaData.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "turbulent_viscosity");
  dualNodalVolume_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  Gju_ = metaData.get_field<GenericFieldType>(
    stk::topology::NODE_RANK, "dudx");

  MasterElement *meSCV = sierra::nalu::get_volume_master_element(AlgTraits::topo_);

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // required fields
  dataPreReqs.add_gathered_nodal_field(*tkeNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*tvisc_, 1);
  dataPreReqs.add_gathered_nodal_field(*dualNodalVolume_, 1);
  dataPreReqs.add_gathered_nodal_field(*Gju_, AlgTraits::nDim_, AlgTraits::nDim_);
  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
}

template<typename AlgTraits>
TurbKineticEnergyKsgsSrcElemKernel<AlgTraits>::~TurbKineticEnergyKsgsSrcElemKernel()
{}

template<typename AlgTraits>
void
TurbKineticEnergyKsgsSrcElemKernel<AlgTraits>::execute(
  SharedMemView<double **>&lhs,
  SharedMemView<double *>&rhs,
  ScratchViews& scratchViews)
{
  SharedMemView<double*>& v_tkeNp1 = scratchViews.get_scratch_view_1D(
    *tkeNp1_);
  SharedMemView<double*>& v_densityNp1 = scratchViews.get_scratch_view_1D(
    *densityNp1_);
  SharedMemView<double*>& v_tvisc = scratchViews.get_scratch_view_1D(
    *tvisc_);
  SharedMemView<double*>& v_dualNodalVolume = scratchViews.get_scratch_view_1D(
    *dualNodalVolume_);
  SharedMemView<double***>& v_Gju = scratchViews.get_scratch_view_3D(*Gju_);
  SharedMemView<double*>& v_scv_volume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;

  for (int ip=0; ip < AlgTraits::numScvIp_; ++ip) {
    const int nearestNode = ipNodeMap_[ip];

    // filter
    double filter = std::pow(v_dualNodalVolume(nearestNode), 1.0/AlgTraits::nDim_);

    double Pk = 0.0;
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        Pk += v_Gju(nearestNode,i,j)*(v_Gju(nearestNode,i,j) + v_Gju(nearestNode,j,i));
      }
    }
    Pk *= v_tvisc(nearestNode);

    // save off multiply used values
    const double tke = v_tkeNp1(nearestNode);
    const double rho = v_densityNp1(nearestNode);

    // dissipation and production
    double Dk = cEps_*rho*std::pow(tke, 1.5)/filter;
    if ( Pk > tkeProdLimitRatio_*Dk )
      Pk = tkeProdLimitRatio_*Dk;
    
    // lhs assembly, all lumped
    const double scvol = v_scv_volume(ip);
    rhs(nearestNode) += (Pk - Dk)*scvol;
    lhs(nearestNode,nearestNode) += 1.5*cEps_*rho*std::sqrt(tke)/filter*scvol;
  }
}

INSTANTIATE_KERNEL(TurbKineticEnergyKsgsSrcElemKernel);

}  // nalu
}  // sierra
