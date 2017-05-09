/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "ScalarAdvDiffElemKernel.h"
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
ScalarAdvDiffElemKernel<AlgTraits>::ScalarAdvDiffElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ScalarFieldType* scalarQ,
  ScalarFieldType* diffFluxCoeff,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    scalarQ_(scalarQ),
    diffFluxCoeff_(diffFluxCoeff),
    lrscv_(sierra::nalu::get_surface_master_element(AlgTraits::topo_)->adjacentNodes())
{
  // Save of required fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  coordinates_ = metaData.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, solnOpts.get_coordinates_name());
  massFlowRate_ = metaData.get_field<GenericFieldType>(
    stk::topology::ELEMENT_RANK, "mass_flow_rate_scs");

  MasterElement *meSCS = sierra::nalu::get_surface_master_element(AlgTraits::topo_);
  meSCS->shape_fcn(&v_shape_function_(0,0));

  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_);
  dataPreReqs.add_gathered_nodal_field(*scalarQ_, 1);
  dataPreReqs.add_gathered_nodal_field(*diffFluxCoeff_, 1);
  dataPreReqs.add_master_element_call(SCS_AREAV);
  dataPreReqs.add_master_element_call(SCS_GRAD_OP);
}

template<typename AlgTraits>
ScalarAdvDiffElemKernel<AlgTraits>::~ScalarAdvDiffElemKernel()
{}

template<typename AlgTraits>
void
ScalarAdvDiffElemKernel<AlgTraits>::execute(
  SharedMemView<double**>& lhs,
  SharedMemView<double*>& rhs,
  stk::mesh::Entity element,
  ScratchViews& scratchViews)
{
  SharedMemView<double*>& v_scalarQ = scratchViews.get_scratch_view_1D(*scalarQ_);
  SharedMemView<double*>& v_diffFluxCoeff = scratchViews.get_scratch_view_1D(*diffFluxCoeff_);

  SharedMemView<double**>& v_scs_areav = scratchViews.get_me_views().scs_areav;
  SharedMemView<double***>& v_dndx = scratchViews.get_me_views().dndx;

  // ip data for this element
  const double *mdot = stk::mesh::field_data(*massFlowRate_, element);

  // start the assembly
  for ( int ip = 0; ip < AlgTraits::numScsIp_; ++ip ) {

    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];

    // save off mdot
    const double tmdot = mdot[ip];

    // compute ip property and
    double diffFluxCoeffIp = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const double r = v_shape_function_(ip,ic);
      diffFluxCoeffIp += r*v_diffFluxCoeff(ic);
    }

    // advection and diffusion
    double qAdv = 0.0;
    double qDiff = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {

      // advection
      const double r = v_shape_function_(ip,ic);
      const double lhsfacAdv = r*tmdot;
      qAdv += lhsfacAdv*v_scalarQ(ic);

      // diffusion
      double lhsfacDiff = 0.0;
      for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
        lhsfacDiff += -diffFluxCoeffIp*v_dndx(ip,ic,j)*v_scs_areav(ip,j);
      }
      qDiff += lhsfacDiff*v_scalarQ(ic);

      // lhs; il then ir
      lhs(il,ic) += lhsfacAdv + lhsfacDiff;
      lhs(ir,ic) -= lhsfacAdv + lhsfacDiff;
    }

    // rhs; il then ir
    rhs(il) -= qAdv + qDiff;
    rhs(ir) += qAdv + qDiff;
  }
}

INSTANTIATE_KERNEL(ScalarAdvDiffElemKernel);

}  // nalu
}  // sierra
