/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "user_functions/SteadyTaylorVortexContinuitySrcElemKernel.h"
#include "AlgTraits.h"
#include "master_element/MasterElement.h"
#include "SolutionOptions.h"
#include "TimeIntegrator.h"

// template and scratch space
#include "BuildTemplates.h"
#include "ScratchViews.h"

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

#include <cmath>

namespace sierra {
namespace nalu {

template<typename AlgTraits>
SteadyTaylorVortexContinuitySrcElemKernel<AlgTraits>::SteadyTaylorVortexContinuitySrcElemKernel(
  const stk::mesh::BulkData& bulkData,
  SolutionOptions& solnOpts,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    ipNodeMap_(sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_)->ipNodeMap()),
    rhoP_(1.0),
    rhoS_(1.0),
    unot_(1.0),
    vnot_(1.0),
    znot_(1.0),
    pnot_(1.0),
    a_(20.0),
    amf_(10.0),
    Sc_(0.9),
    pi_(acos(-1.0)),
    projTimeScale_(1.0)
{
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  coordinates_ = metaData.get_field<double>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(AlgTraits::topo_);
  get_scv_shape_fn_data<AlgTraits>([&](double* ptr){meSCV->shape_fcn(ptr);}, v_shape_function_);

  // add master elements
  dataPreReqs.add_cvfem_volume_me(meSCV);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(SCV_VOLUME, CURRENT_COORDINATES);
}

template<typename AlgTraits>
void
SteadyTaylorVortexContinuitySrcElemKernel<AlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
  projTimeScale_ = timeIntegrator.get_time_step() / timeIntegrator.get_gamma1();
}

template<typename AlgTraits>
void
SteadyTaylorVortexContinuitySrcElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType**>& /* lhs */,
  SharedMemView<DoubleType *>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  // Forcing nDim = 3 instead of using AlgTraits::nDim_ here to avoid compiler
  // warnings when this template is instantiated for 2-D topologies. 
  DoubleType w_scvCoords[3];

  SharedMemView<DoubleType**>& v_coordinates = scratchViews.get_scratch_view_2D(*coordinates_);
  SharedMemView<DoubleType*>& v_scv_volume = scratchViews.get_me_views(CURRENT_COORDINATES).scv_volume;

  // interpolate to ips and evaluate source
  for ( int ip = 0; ip < AlgTraits::numScvIp_; ++ip ) {

    // nearest node to ip
    const int nearestNode = ipNodeMap_[ip];

    // zero out
    for ( int j =0; j < AlgTraits::nDim_; ++j )
      w_scvCoords[j] = 0.0;

    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const DoubleType r = v_shape_function_(ip,ic);
      for ( int j = 0; j < AlgTraits::nDim_; ++j )
        w_scvCoords[j] += r*v_coordinates(ic,j);
    }

    // save off scv coords and src
    const DoubleType x = w_scvCoords[0];
    const DoubleType y = w_scvCoords[1];

    const DoubleType src = 0.10e1 * stk::math::pow(znot_ * stk::math::cos(amf_ * pi_ * x) * stk::math::cos(amf_ * pi_ * y) / rhoP_ + (0.1e1 - znot_ * stk::math::cos(amf_ * pi_ * x) * stk::math::cos(amf_ * pi_ * y)) / rhoS_, -0.2e1) * unot_ * stk::math::cos(a_ * pi_ * x) * stk::math::sin(a_ * pi_ * y) * (-znot_ * stk::math::sin(amf_ * pi_ * x) * amf_ * pi_ * stk::math::cos(amf_ * pi_ * y) / rhoP_ + znot_ * stk::math::sin(amf_ * pi_ * x) * amf_ * pi_ * stk::math::cos(amf_ * pi_ * y) / rhoS_) + 0.10e1 / (znot_ * stk::math::cos(amf_ * pi_ * x) * stk::math::cos(amf_ * pi_ * y) / rhoP_ + (0.1e1 - znot_ * stk::math::cos(amf_ * pi_ * x) * stk::math::cos(amf_ * pi_ * y)) / rhoS_) * unot_ * stk::math::sin(a_ * pi_ * x) * a_ * pi_ * stk::math::sin(a_ * pi_ * y) - 0.10e1 * stk::math::pow(znot_ * stk::math::cos(amf_ * pi_ * x) * stk::math::cos(amf_ * pi_ * y) / rhoP_ + (0.1e1 - znot_ * stk::math::cos(amf_ * pi_ * x) * stk::math::cos(amf_ * pi_ * y)) / rhoS_, -0.2e1) * vnot_ * stk::math::sin(a_ * pi_ * x) * stk::math::cos(a_ * pi_ * y) * (-znot_ * stk::math::cos(amf_ * pi_ * x) * stk::math::sin(amf_ * pi_ * y) * amf_ * pi_ / rhoP_ + znot_ * stk::math::cos(amf_ * pi_ * x) * stk::math::sin(amf_ * pi_ * y) * amf_ * pi_ / rhoS_) - 0.10e1 / (znot_ * stk::math::cos(amf_ * pi_ * x) * stk::math::cos(amf_ * pi_ * y) / rhoP_ + (0.1e1 - znot_ * stk::math::cos(amf_ * pi_ * x) * stk::math::cos(amf_ * pi_ * y)) / rhoS_) * vnot_ * stk::math::sin(a_ * pi_ * x) * stk::math::sin(a_ * pi_ * y) * a_ * pi_;
  
    rhs(nearestNode) += src*v_scv_volume(ip)/projTimeScale_;
  }
}

INSTANTIATE_KERNEL(SteadyTaylorVortexContinuitySrcElemKernel);

}  // nalu
}  // sierra
