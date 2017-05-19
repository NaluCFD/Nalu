/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "ScalarDiffFemKernel.h"
#include "AlgTraits.h"
#include "master_element/Hex8FEM.h"
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
ScalarDiffFemKernel<AlgTraits>::ScalarDiffFemKernel(
  const stk::mesh::BulkData& bulkData,
  ScalarFieldType* scalarQ,
  ScalarFieldType* diffFluxCoeff,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    bulkData_(&bulkData),
    scalarQ_(scalarQ),
    diffFluxCoeff_(diffFluxCoeff),
    meFEM_(new Hex8FEM()),
    ipWeight_(&meFEM_->weights_[0])
{
  ThrowRequireMsg(AlgTraits::topo_ == stk::topology::HEX_8,
                  "FEM_DIFF only available for hexes currently");

  // Save of required fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  coordinates_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates");

  // master element
  meFEM_->shape_fcn(&v_shape_function_(0,0));
  dataPreReqs.add_fem_volume_me(meFEM_);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*scalarQ_, 1);
  dataPreReqs.add_gathered_nodal_field(*diffFluxCoeff_, 1);
  dataPreReqs.add_master_element_call(FEM_GRAD_OP, CURRENT_COORDINATES);
}

template<typename AlgTraits>
ScalarDiffFemKernel<AlgTraits>::~ScalarDiffFemKernel()
{
  delete meFEM_;
}

template<typename AlgTraits>
void
ScalarDiffFemKernel<AlgTraits>::execute(
  SharedMemView<double**>& lhs,
  SharedMemView<double*>& rhs,
  stk::mesh::Entity /* element */,
  ScratchViews& scratchViews)
{
  SharedMemView<double*>& v_scalarQ = scratchViews.get_scratch_view_1D(*scalarQ_);
  SharedMemView<double*>& v_diffFluxCoeff = scratchViews.get_scratch_view_1D(*diffFluxCoeff_);
  
  SharedMemView<double***>& v_dndx = scratchViews.get_me_views(CURRENT_COORDINATES).dndx_fem;
  SharedMemView<double*>& v_det_j = scratchViews.get_me_views(CURRENT_COORDINATES).det_j_fem;

  for ( int ip = 0; ip < AlgTraits::numGp_; ++ip ) {

    // compute ip property
    double diffFluxCoeffIp = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const double r = v_shape_function_(ip,ic);
      diffFluxCoeffIp += r*v_diffFluxCoeff(ic);
    }

    // start the assembly
    const double ipFactor = v_det_j(ip)*ipWeight_[ip];

    // row ir
    for ( int ir = 0; ir < AlgTraits::nodesPerElement_; ++ir) {

      // column ic
      double rhsSum = 0.0;
      for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {

        double lhsSum = 0.0;
        double scalarQ = v_scalarQ(ic);
        for ( int j = 0; j < AlgTraits::nDim_; ++j ) {
          const double fac = v_dndx(ip,ir,j)*diffFluxCoeffIp*v_dndx(ip,ic,j);
          lhsSum += fac;
          rhsSum += fac*scalarQ;
        }
        lhs(ir,ic) += lhsSum*ipFactor;
      }
      rhs(ir) -= rhsSum*ipFactor;
    }
  }
}

INSTANTIATE_KERNEL(ScalarDiffFemKernel);

}  // nalu
}  // sierra
