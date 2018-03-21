/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <TurbViscSmagorinskyAlgorithm.h>
#include <Algorithm.h>
#include <FieldTypeDef.h>
#include <Realm.h>

#include <kernel/UtilityKernelManager.h>
#include <kernel/NodalUtilityKernelManager.h>
#include <kernel/SmagorinskyViscosityUtilityKernel.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// TurbViscSmagorinskyAlgorithm - compute tvisc for Smagorinsky model
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TurbViscSmagorinskyAlgorithm::TurbViscSmagorinskyAlgorithm(
  Realm &realm,
  stk::mesh::Part *part)
  : Algorithm(realm, part),
    dudx_(NULL),
    density_(NULL),
    tvisc_(NULL),
    dualNodalVolume_(NULL),
    cmuCs_(realm.get_turb_model_constant(TM_cmuCs))
{
  kernelManager_ = UtilityKernelManager::create(
    stk::topology::NODE,
    realm_.bulk_data(),
    *realm_.solutionOptions_,
    *realm_.timeIntegrator_,
    in,
    out
  );

  kernelManager_->build_utility_kernel<SmagorinskyViscosity>("smag_visc",
    realm_.bulk_data(), *realm_.solutionOptions_, in, out);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
TurbViscSmagorinskyAlgorithm::execute()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  stk::mesh::Selector s_all_nodes = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
      & stk::mesh::selectField(*tvisc_);

  kernelManager_->execute_kernels(meta_data.universal_part());
}

} // namespace nalu
} // namespace Sierra
