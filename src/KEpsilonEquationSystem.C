/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "KEpsilonEquationSystem.h"
#include "AlgorithmDriver.h"
#include "ComputeSSTMaxLengthScaleElemAlgorithm.h"
#include "FieldFunctions.h"
#include "master_element/MasterElement.h"
#include "NaluEnv.h"
#include "TurbDissipationEquationSystem.h"
#include "SolutionOptions.h"
#include "TurbKineticEnergyEquationSystem.h"
#include "Realm.h"

// stk_util
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/MetaData.hpp>

// stk_io
#include <stk_io/IossBridge.hpp>

// basic c++
#include <cmath>
#include <vector>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// KEpsilonEquationSystem - manage SST
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
KEpsilonEquationSystem::KEpsilonEquationSystem(
 EquationSystems& eqSystems,
 const bool outputClippingDiag)
  : EquationSystem(eqSystems, "KEpsilonWrap"),
    outputClippingDiag_(outputClippingDiag),
    tkeEqSys_(NULL),
    epsEqSys_(NULL),
    tke_(NULL),
    eps_(NULL),
    isInit_(true)
{
  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  // create momentum and pressure
  tkeEqSys_= new TurbKineticEnergyEquationSystem(eqSystems);
  epsEqSys_ = new TurbDissipationEquationSystem(eqSystems);
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
KEpsilonEquationSystem::~KEpsilonEquationSystem()
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
KEpsilonEquationSystem::initialize()
{
  // let equation systems that are owned some information
  tkeEqSys_->convergenceTolerance_ = convergenceTolerance_;
  epsEqSys_->convergenceTolerance_ = convergenceTolerance_;
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
KEpsilonEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();
  const int numStates = realm_.number_of_states();

  // re-register tke and sdr for convenience; other specifics managed by EQS
  tke_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_ke", numStates));
  stk::mesh::put_field_on_mesh(*tke_, *part, nullptr);
  eps_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_dissipation", numStates));
  stk::mesh::put_field_on_mesh(*eps_, *part, nullptr);
}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
KEpsilonEquationSystem::solve_and_update()
{
  if ( isInit_ ) {
    // compute projected nodal gradients
    tkeEqSys_->compute_projected_nodal_gradient();
    epsEqSys_->compute_projected_nodal_gradient();
    isInit_ = false;
  }

  // k-eps effective viscosity for k and omega
  tkeEqSys_->compute_effective_diff_flux_coeff();
  epsEqSys_->compute_effective_diff_flux_coeff();

  // wall values
  tkeEqSys_->compute_wall_model_parameters();
  epsEqSys_->compute_wall_model_parameters();

  // start the iteration loop
  for ( int k = 0; k < maxIterations_; ++k ) {

    NaluEnv::self().naluOutputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << name_ << std::endl;

    // tke and sdr assemble, load_complete and solve; Jacobi iteration
    tkeEqSys_->assemble_and_solve(tkeEqSys_->kTmp_);
    epsEqSys_->assemble_and_solve(epsEqSys_->eTmp_);
    
    // update each
    update_and_clip();

    // compute projected nodal gradients
    tkeEqSys_->compute_projected_nodal_gradient();
    epsEqSys_->compute_projected_nodal_gradient();
  }

}

//--------------------------------------------------------------------------
//-------- initial_work ----------------------------------------------------
//--------------------------------------------------------------------------
void
KEpsilonEquationSystem::initial_work()
{
  // do not let the user specify a negative field
  const double clipValue = 1.0e-8;
  const double viscFac = 1.0;
  const double cMu = realm_.get_turb_model_constant(TM_cMu);

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // required fields
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  ScalarFieldType *viscosity = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity");

  // required fields with state
  ScalarFieldType &epsNp1 = eps_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &tkeNp1 = tke_->field_of_state(stk::mesh::StateNP1);

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*eps_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    const double *visc = stk::mesh::field_data(*viscosity, b);
    const double *rho = stk::mesh::field_data(*density, b);
    double *tke = stk::mesh::field_data(tkeNp1, b);
    double *eps = stk::mesh::field_data(epsNp1, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      const double tkeNew = tke[k];
      const double epsNew = eps[k];
      const double viscM = visc[k]*viscFac;
      if ( (tkeNew >= 0.0) && (epsNew > 0.0) ) {
        // nothing
      }
      else if ( (tkeNew < 0.0) && (epsNew < 0.0) ) {
        // both negative;
        tke[k] = clipValue;
        eps[k] = cMu*rho[k]*clipValue*clipValue/viscM;
      }
      else if ( tkeNew < 0.0 ) {
        tke[k] = std::sqrt(viscM*epsNew/(rho[k]*cMu));
        eps[k] = epsNew;
      }
      else {
        eps[k] = cMu*rho[k]*tkeNew*tkeNew/viscM;
        tke[k] = tkeNew;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- update_and_clip() -----------------------------------------------
//--------------------------------------------------------------------------
void
KEpsilonEquationSystem::update_and_clip()
{
  const double clipTke = 1.0e-8;
  const double viscFac = 1.0;
  const double cMu = realm_.get_turb_model_constant(TM_cMu);
  
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // clip diagnosis
  size_t numClip[2] = {0,0};
  double minTke = +1.0e16;
  double minEps = +1.0e16;
  const double small = 1.0e-16;

  // required fields
  ScalarFieldType *density = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "density");
  ScalarFieldType *viscosity = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity");
  ScalarFieldType *turbViscosity = meta_data.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity");

  // required fields with state
  ScalarFieldType &epsNp1 = eps_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &tkeNp1 = tke_->field_of_state(stk::mesh::StateNP1);

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*turbViscosity);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    const double *visc = stk::mesh::field_data(*viscosity, b);
    const double *rho = stk::mesh::field_data(*density, b);
    const double *kTmp = stk::mesh::field_data(*tkeEqSys_->kTmp_, b);
    const double *eTmp = stk::mesh::field_data(*epsEqSys_->eTmp_, b);
    double *tke = stk::mesh::field_data(tkeNp1, b);
    double *eps = stk::mesh::field_data(epsNp1, b);
    double *tvisc = stk::mesh::field_data(*turbViscosity, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      const double tkeNew = tke[k] + kTmp[k];
      const double epsNew = eps[k] + eTmp[k];
      const double viscM = visc[k]*viscFac;

      if ( (tkeNew > small ) && (epsNew > small) ) {
        // if all is well
        tke[k] = tkeNew;
        eps[k] = epsNew;
      }
      else if ( (tkeNew < small) && (epsNew < small) ) {
        // both negative; set k to small, tvisc to molecular visc*fac and use Prandtl/Kolm for eps
        tke[k] = clipTke;
        eps[k] = cMu*rho[k]*clipTke*clipTke/viscM;
        minTke = std::min(minTke, tkeNew);
        minEps = std::min(minEps, epsNew);
        numClip[0]++;
        numClip[1]++;
        tvisc[k] = viscM;
      }
      else if ( tkeNew < small ) {
        // only tke is off; reset tvisc to molecular visc*fac and compute new tke appropriately
        tke[k] = std::sqrt(viscM*epsNew/(rho[k]*cMu));
        eps[k] = epsNew;
        minTke = std::min(minTke, tkeNew);
        numClip[0]++;
        tvisc[k] = visc[k];
      }
      else {
        // only eps if off; reset tvisc to molecular visc and compute new eps appropriately
        eps[k] = cMu*rho[k]*tkeNew*tkeNew/viscM;
        tke[k] = tkeNew;
        minEps = std::min(minEps, epsNew);
        numClip[1]++;
        tvisc[k] = viscM;
      }
    }
  }

  // output clipping
  if ( outputClippingDiag_ ) {
    size_t g_numClip[2] = {};
    stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
    stk::all_reduce_sum(comm, numClip, g_numClip, 2);
    
    if ( g_numClip[0] > 0 ) {
      double g_minTke = 0;
      stk::all_reduce_min(comm, &minTke, &g_minTke, 1);
      NaluEnv::self().naluOutputP0() << "TKE clipped (-) " << g_numClip[0] << " times; min: " << g_minTke << std::endl;
    }
    if ( g_numClip[1] > 0 ) {
      double g_minEps = 0;
      stk::all_reduce_min(comm, &minEps, &g_minEps, 1);
      NaluEnv::self().naluOutputP0() << "EPS clipped (-) " << g_numClip[1] << " times; min: " << g_minEps << std::endl;
    }
  }

}
  
} // namespace nalu
} // namespace Sierra
