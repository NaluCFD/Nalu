/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/HeatCondMassFemKernel.h"
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
HeatCondMassFemKernel<AlgTraits>::HeatCondMassFemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ScalarFieldType* temperature,
  ScalarFieldType* density,
  ScalarFieldType* specHeat,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    bulkData_(&bulkData),
    density_(density),
    specHeat_(specHeat)
{
  // Save of required fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  coordinates_ = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, solnOpts.get_coordinates_name());

  temperatureN_ = &(temperature->field_of_state(stk::mesh::StateN));
  temperatureNp1_ = &(temperature->field_of_state(stk::mesh::StateNP1));
  if (temperature->number_of_states() == 2)
    temperatureNm1_ = temperatureN_;
  else
    temperatureNm1_ = &(temperature->field_of_state(stk::mesh::StateNM1));

  // extract master element
  MasterElement *meFEM = sierra::nalu::MasterElementRepo::get_fem_master_element(AlgTraits::topo_);
  
  // copy ip weights into our 1-d view
  for ( int k = 0; k < AlgTraits::numGp_; ++k )
    v_ip_weight_[k] = meFEM->weights_[k];
  
  // master element, shape function is shifted consistently
  if ( solnOpts.get_shifted_grad_op(temperature->name()) )
    get_fem_shape_fn_data<AlgTraits>([&](double* ptr){meFEM->shifted_shape_fcn(ptr);}, v_shape_function_);
  else
    get_fem_shape_fn_data<AlgTraits>([&](double* ptr){meFEM->shape_fcn(ptr);}, v_shape_function_);
  
  // add FEM master element
  dataPreReqs.add_fem_volume_me(meFEM);

  // fields and data
  dataPreReqs.add_coordinates_field(*coordinates_, AlgTraits::nDim_, CURRENT_COORDINATES);
  dataPreReqs.add_gathered_nodal_field(*temperatureNp1_, 1);
  dataPreReqs.add_gathered_nodal_field(*temperatureN_, 1);
  dataPreReqs.add_gathered_nodal_field(*temperatureNm1_, 1);
  dataPreReqs.add_gathered_nodal_field(*density_, 1);
  dataPreReqs.add_gathered_nodal_field(*specHeat_, 1);
  
  //dataPreReqs.add_master_element_call(FEM_GRAD_OP, CURRENT_COORDINATES);
  dataPreReqs.add_master_element_call(FEM_DET_J, CURRENT_COORDINATES);
}

template<typename AlgTraits>
HeatCondMassFemKernel<AlgTraits>::~HeatCondMassFemKernel()
{
  // does nothing
}

template<typename AlgTraits>
void
HeatCondMassFemKernel<AlgTraits>::setup(const TimeIntegrator& timeIntegrator)
{
  dt_ = timeIntegrator.get_time_step();
  gamma1_ = timeIntegrator.get_gamma1();
  gamma2_ = timeIntegrator.get_gamma2();
  gamma3_ = timeIntegrator.get_gamma3(); // gamma3 may be zero
}

template<typename AlgTraits>
void
HeatCondMassFemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType**>& lhs,
  SharedMemView<DoubleType*>& rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  SharedMemView<DoubleType*>& v_temperatureNp1 = scratchViews.get_scratch_view_1D(*temperatureNp1_);
  SharedMemView<DoubleType*>& v_temperatureN = scratchViews.get_scratch_view_1D(*temperatureN_);
  SharedMemView<DoubleType*>& v_temperatureNm1 = scratchViews.get_scratch_view_1D(*temperatureNm1_);
  SharedMemView<DoubleType*>& v_density = scratchViews.get_scratch_view_1D(*density_);
  SharedMemView<DoubleType*>& v_specHeat = scratchViews.get_scratch_view_1D(*specHeat_);
  
  SharedMemView<DoubleType*>& v_det_j = scratchViews.get_me_views(CURRENT_COORDINATES).det_j_fem;

  for ( int ip = 0; ip < AlgTraits::numGp_; ++ip ) {

    // compute ip property
    DoubleType densityIp = 0.0;
    DoubleType specHeatIp = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
      const DoubleType r = v_shape_function_(ip,ic);
      densityIp += r*v_density(ic);
      specHeatIp += r*v_specHeat(ic);
    }

    // start the assembly (collect ip scalings)
    const DoubleType ipFactor = densityIp*specHeatIp*v_det_j(ip)*v_ip_weight_(ip);

    // row ir
    for ( int ir = 0; ir < AlgTraits::nodesPerElement_; ++ir) {
      
      const DoubleType wIr = v_shape_function_(ip,ir);
      
      // column ic
      DoubleType rhsSum = 0.0;
      for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {
        DoubleType tdot = (gamma1_*v_temperatureNp1(ic) + gamma2_*v_temperatureN(ic) + gamma3_*v_temperatureNm1(ic))/dt_;
        const DoubleType fac = wIr*v_shape_function_(ip,ic);
        rhsSum += fac*tdot; 
        lhs(ir,ic) += gamma1_*ipFactor*fac/dt_;
      }
      rhs(ir) -= rhsSum*ipFactor;
    }
  }
}

INSTANTIATE_FEM_KERNEL(HeatCondMassFemKernel);

}  // nalu
}  // sierra
