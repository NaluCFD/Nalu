
/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/TurbDissipationKEpsilonSrcElemKernel.h"
#include "FieldTypeDef.h"
#include "SolutionOptions.h"

#include "BuildTemplates.h"
#include "ScratchViews.h"

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

// c++
#include<limits>

namespace sierra {
namespace nalu {

template <typename AlgTraits>
TurbDissipationKEpsilonSrcElemKernel<AlgTraits>::
  TurbDissipationKEpsilonSrcElemKernel(
    const stk::mesh::BulkData& bulkData,
    const SolutionOptions& solnOpts,
    ElemDataRequests& dataPreReqs,
    const bool lumpedMass)
  : Kernel(),
    lumpedMass_(lumpedMass),
    shiftedGradOp_(solnOpts.get_shifted_grad_op("velocity")),
    cEpsOne_(solnOpts.get_turb_model_constant(TM_cEpsOne)),
    cEpsTwo_(solnOpts.get_turb_model_constant(TM_cEpsTwo)),
    tkeProdLimitRatio_(solnOpts.get_turb_model_constant(TM_tkeProdLimitRatio)),
    includeDivU_(solnOpts.includeDivU_),
    twoThirds_(2.0/3.0),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap())
{
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  coordinates_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  VectorFieldType* velocity = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  velocityNp1_ = &(velocity->field_of_state(stk::mesh::StateNP1));
  ScalarFieldType* tke = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_ke");
  tkeNp1_ = &tke->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType* eps = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_dissipation");
  epsNp1_ = &eps->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType* density = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  densityNp1_ = &density->field_of_state(stk::mesh::StateNP1);
  tvisc_ = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity");

  MasterElement* meSCV =
    sierra::nalu::MasterElementRepo::get_volume_master_element(
      AlgTraits::topo_);

  // compute shape function
  if (lumpedMass_)
    get_scv_shape_fn_data<AlgTraits>(
      [&](double* ptr) { meSCV->shifted_shape_fcn(ptr); }, v_shape_function_);
  else
    get_scv_shape_fn_data<AlgTraits>(
      [&](double* ptr) { meSCV->shape_fcn(ptr); }, v_shape_function_);

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*velocityNp1_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*tkeNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*epsNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*densityNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*tvisc_, 1);
  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
  if (shiftedGradOp_)
    dataPreReqs.add_master_element_call(
      SCV_SHIFTED_GRAD_OP, CURRENT_COORDINATES);
  else
    dataPreReqs.add_master_element_call(SCV_GRAD_OP, CURRENT_COORDINATES);
}

template <typename AlgTraits>
TurbDissipationKEpsilonSrcElemKernel<
  AlgTraits>::~TurbDissipationKEpsilonSrcElemKernel()
{
}

template <typename AlgTraits>
void
TurbDissipationKEpsilonSrcElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType**>& lhs,
  SharedMemView<DoubleType*>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  NALU_ALIGNED DoubleType w_dudx[AlgTraits::nDim_][AlgTraits::nDim_];

  SharedMemView<DoubleType**>& v_velocityNp1 = scratchViews.get_scratch_view_2D(*velocityNp1_);
  SharedMemView<DoubleType*>& v_tkeNp1 = scratchViews.get_scratch_view_1D(*tkeNp1_);
  SharedMemView<DoubleType*>& v_epsNp1 = scratchViews.get_scratch_view_1D(*epsNp1_);
  SharedMemView<DoubleType*>& v_densityNp1 = scratchViews.get_scratch_view_1D(*densityNp1_);
  SharedMemView<DoubleType*>& v_tvisc = scratchViews.get_scratch_view_1D(*tvisc_);
  SharedMemView<DoubleType***>& v_dndx =
    shiftedGradOp_
      ? scratchViews.get_me_views(CURRENT_COORDINATES).dndx_scv_shifted
      : scratchViews.get_me_views(CURRENT_COORDINATES).dndx_scv;
  SharedMemView<DoubleType*>& v_scv_volume =
    scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;

  for (int ip = 0; ip < AlgTraits::numScvIp_; ++ip) {

    // nearest node to ip
    const int nearestNode = ipNodeMap_[ip];

    // save off scvol
    const DoubleType scV = v_scv_volume(ip);

    DoubleType epsIp = 0.0;
    DoubleType tkeIp = 0.0;
    DoubleType rhoIp = 0.0;
    DoubleType tviscIp = 0.0;
    for (int i = 0; i < AlgTraits::nDim_; ++i) {
      for (int j = 0; j < AlgTraits::nDim_; ++j) {
        w_dudx[i][j] = 0.0;
      }
    }

    for (int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic) {

      const DoubleType r = v_shape_function_(ip, ic);

      epsIp += r*v_epsNp1(ic);
      tkeIp += r* v_tkeNp1(ic);
      rhoIp += r*v_densityNp1(ic);
      tviscIp += r*v_tvisc(ic);

      for (int i = 0; i < AlgTraits::nDim_; ++i) {
        const DoubleType ui = v_velocityNp1(ic, i);
        for (int j = 0; j < AlgTraits::nDim_; ++j) {
          w_dudx[i][j] += v_dndx(ip, ic, j) * ui;
        }
      }
    }

    DoubleType divU = 0.0;
    for ( int i = 0; i < AlgTraits::nDim_; ++i ) {
      divU += w_dudx[i][i];
    }
    divU *= includeDivU_;
    
    DoubleType Pk = 0.0;
    for (int i = 0; i < AlgTraits::nDim_; ++i) {
      for (int j = 0; j < AlgTraits::nDim_; ++j) {
        Pk += w_dudx[i][j]*(w_dudx[i][j] + w_dudx[j][i]);
      }
    }
    Pk = tviscIp*(Pk-twoThirds_*divU*divU) - twoThirds_*rhoIp*tkeIp*divU;

    // clip
    Pk = stk::math::max(0.0,Pk);
    const DoubleType tkeIpC 
      = stk::math::max(std::numeric_limits<double>::min(),tkeIp);

    // dissipation and production (limited)
    DoubleType Dk = rhoIp*epsIp;
    Pk = stk::math::min(Pk, tkeProdLimitRatio_*Dk);

    // assemble RHS and LHS
    rhs(nearestNode) += epsIp/tkeIpC*(cEpsOne_*Pk - cEpsTwo_*Dk)*scV;
    for (int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic) {
      lhs(nearestNode, ic) += v_shape_function_(ip,ic)*cEpsTwo_*rhoIp*epsIp/tkeIpC*scV;
    }
  }
}

INSTANTIATE_KERNEL(TurbDissipationKEpsilonSrcElemKernel);

} // namespace nalu
} // namespace sierra
