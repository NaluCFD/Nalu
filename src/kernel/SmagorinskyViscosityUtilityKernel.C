/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <kernel/SmagorinskyViscosityUtilityKernel.h>

#include <EquationSystem.h>
#include <SolverAlgorithm.h>
#include <master_element/MasterElement.h>

#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <Realm.h>
#include <kernel/Kernel.h>
#include <TimeIntegrator.h>
#include <SolutionOptions.h>

#include <master_element/TensorOps.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk topo
#include <stk_topology/topology.hpp>

#include <KokkosInterface.h>
#include <ScratchViews.h>
#include <CopyAndInterleave.h>

namespace sierra{
namespace nalu{

SmagorinskyViscosity::SmagorinskyViscosity(
  const stk::mesh::BulkData& bulk,
  const SolutionOptions& solnOpts,
  ElemDataRequests& inputData,
  ElemDataRequests& outputData)
: NodalUtilityKernel()
{
  const auto& meta = bulk.mesh_meta_data();
  dudx_ = meta.get_field<GenericFieldType>(stk::topology::NODE_RANK, "dudx");
  dualVolume_ = meta.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  density_ = meta.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  turbVisc_ = meta.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity");

  nDim_ = meta.spatial_dimension();
  inputData.add_gathered_nodal_field(*dudx_, nDim_, nDim_);
  inputData.add_gathered_nodal_field(*dualVolume_,1);
  inputData.add_gathered_nodal_field(*density_,1);
  outputData.add_gathered_nodal_field(*turbVisc_, 1);

  Cs_ = solnOpts.turbModelConstantMap_.at(TM_cmuCs);
}

void SmagorinskyViscosity::execute(
  const NodalScratchData<DoubleType>& inputData,
  NodalScratchData<DoubleType>& outputData) const
{
  const SharedMemView<const DoubleType**>& v_dudx = inputData.get_2D(*dudx_);
  const DoubleType SijMag = mag_symm_part(nDim_, v_dudx);
  const DoubleType filterSize = stk::math::if_then_else(nDim_ == 3,
    stk::math::cbrt(inputData.get_value(*dualVolume_)),
    stk::math::sqrt(inputData.get_value(*dualVolume_))
  );
  const DoubleType rho = inputData.get_value(*density_);

  outputData.get_value(*turbVisc_) = rho * pow2(Cs_ * filterSize) * SijMag;
}

} // namespace nalu
} // namespace Sierra
