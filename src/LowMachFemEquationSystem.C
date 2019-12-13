/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "LowMachFemEquationSystem.h"
#include "AlgorithmDriver.h"
#include "AuxFunctionAlgorithm.h"
#include "ConstantAuxFunction.h"
#include "CopyFieldAlgorithm.h"
#include "DirichletBC.h"
#include "Enums.h"
#include "EquationSystem.h"
#include "EquationSystems.h"
#include "FieldFunctions.h"
#include "LinearSolver.h"
#include "LinearSolvers.h"
#include "LinearSystem.h"
#include "master_element/MasterElement.h"
#include "NaluEnv.h"
#include "NaluParsing.h"
#include "ProjectedNodalGradientEquationSystem.h"
#include "Realm.h"
#include "Realms.h"
#include "Simulation.h"
#include "SolutionOptions.h"
#include "SolverAlgorithmDriver.h"

// classic algorithms
#include "AssembleCourantReynoldsFemAlgorithm.h"

// algorithm drivers

// template for kernels
#include "AlgTraits.h"
#include "kernel/KernelBuilder.h"
#include "kernel/KernelBuilderLog.h"

// kernels
#include "kernel/ContinuityAdvFemKernel.h"
#include "kernel/MomentumAdvFemKernel.h"
#include "kernel/MomentumDiffFemKernel.h"
#include "kernel/MomentumMassFemKernel.h"

// source kernels
#include "kernel/MomentumBodyForceFemKernel.h"

// bc kernels
#include "kernel/ContinuityOpenFemKernel.h"
#include "kernel/ContinuityInflowFemKernel.h"
#include "kernel/MomentumOpenAdvDiffFemKernel.h"

// nso - n/a

// hybrid turbulence -n/a

// overset
#include "overset/UpdateOversetFringeAlgorithmDriver.h"

// user function
#include "user_functions/TaylorGreenPressureAuxFunction.h"
#include "user_functions/TaylorGreenVelocityAuxFunction.h"
#include "user_functions/WindEnergyTaylorVortexAuxFunction.h"
#include "user_functions/OneTwoTenVelocityAuxFunction.h"
#include "user_functions/PulseVelocityAuxFunction.h"

// stk_util
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/SortAndUnique.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>

// stk_topo
#include <stk_topology/topology.hpp>

#include <utils/StkHelpers.h>

// basic c++
#include <vector>


namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// LowMachFemEquationSystem - manage the low Mach FEM equation system (uvw_p)
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
LowMachFemEquationSystem::LowMachFemEquationSystem(
  EquationSystems& eqSystems)
  : EquationSystem(eqSystems, "LowMachEOSWrap","low_mach_type"),
    density_(NULL),
    viscosity_(NULL),
    dpdxL_(NULL),
    vrtmL_(NULL),
    isInit_(true)
{
  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  // create momentum and pressure
  momentumEqSys_= new MomentumFemEquationSystem(eqSystems);
  continuityEqSys_ = new ContinuityFemEquationSystem(eqSystems);

  // inform realm
  realm_.hasFluids_ = true;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
LowMachFemEquationSystem::~LowMachFemEquationSystem()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
LowMachFemEquationSystem::initialize()
{
  // let equation systems that are owned some information
  momentumEqSys_->convergenceTolerance_ = convergenceTolerance_;
  continuityEqSys_->convergenceTolerance_ = convergenceTolerance_;
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
LowMachFemEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register lagged variables to avoid storage of an (ip) field mDot
  const int nDim = meta_data.spatial_dimension();
  dpdxL_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx_lagged"));
  stk::mesh::put_field_on_mesh(*dpdxL_, *part, nDim, nullptr);
  vrtmL_ = &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "vrtm_lagged"));
  stk::mesh::put_field_on_mesh(*vrtmL_, *part, nDim, nullptr);
 
  // add properties; denisty needs to be a restart field
  const int numStates = realm_.number_of_states();
  density_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "density", numStates));
  stk::mesh::put_field_on_mesh(*density_, *part, nullptr);
  realm_.augment_restart_variable_list("density");

  viscosity_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity"));
  stk::mesh::put_field_on_mesh(*viscosity_, *part, nullptr);

  // push to property list
  realm_.augment_property_map(DENSITY_ID, density_);
  realm_.augment_property_map(VISCOSITY_ID, viscosity_);

  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && (!realm_.restarted_simulation() || realm_.support_inconsistent_restart()) ) {
    ScalarFieldType &densityN = density_->field_of_state(stk::mesh::StateN);
    ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               &densityNp1, &densityN,
                               0, 1,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlg);
  }

  // register the fringe nodal field 
  if ( realm_.query_for_overset() && realm_.has_mesh_motion() ) {
    ScalarFieldType *fringeNode
      = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "fringe_node"));
    stk::mesh::put_field_on_mesh(*fringeNode, *part, nullptr);
  }
}

//--------------------------------------------------------------------------
//-------- register_element_fields -------------------------------------------
//--------------------------------------------------------------------------
void
LowMachFemEquationSystem::register_element_fields(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();
    
  // ip fields - n/a

  // register the intersected elemental field
  if ( realm_.query_for_overset() ) {
    const int sizeOfElemField = 1;
    GenericFieldType *intersectedElement
      = &(meta_data.declare_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "intersected_element"));
    stk::mesh::put_field_on_mesh(*intersectedElement, *part, sizeOfElemField, nullptr);
  }

  // provide mean element Peclet and Courant fields; always...
  GenericFieldType *elemReynolds
    = &(meta_data.declare_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "element_reynolds"));
  stk::mesh::put_field_on_mesh(*elemReynolds, *part, 1, nullptr);
  GenericFieldType *elemCourant
    = &(meta_data.declare_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "element_courant"));
  stk::mesh::put_field_on_mesh(*elemCourant, *part, 1, nullptr);
}

//--------------------------------------------------------------------------
//-------- register_open_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
LowMachFemEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &partTopo,
  const OpenBoundaryConditionData &openBCData)
{  
  // Motivation: register common fields between momentum and continuity

  // set face/element required sort
  realm_.solutionOptions_->set_consolidated_bc_solver_alg();

  // extract types of data
  stk::mesh::MetaData &metaData = realm_.meta_data();
  OpenUserData userData = openBCData.userData_;

  // register boundary data
  ScalarFieldType *pressureBC 
    = &(metaData.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure_bc"));
  stk::mesh::put_field_on_mesh(*pressureBC, *part, nullptr);

  // algorithm to populate specified pressure
  Pressure pSpec = userData.p_;
  std::vector<double> userSpecPbc(1);
  userSpecPbc[0] = pSpec.pressure_;
  
  // new it
  ConstantAuxFunction *theAuxFuncPbc = new ConstantAuxFunction(0, 1, userSpecPbc);
  
  // bc data alg
  AuxFunctionAlgorithm *auxAlgPbc
    = new AuxFunctionAlgorithm(realm_, part,
                               pressureBC, theAuxFuncPbc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlgPbc);

  // integration point data
  MasterElement *meFC = sierra::nalu::MasterElementRepo::get_fem_master_element(partTopo);
  const int numIntPoints = meFC->numIntPoints_;
  
  // pbip; always register (initial value of zero)
  std::vector<double> zeroVec(numIntPoints,0.0);
  GenericFieldType *pBip 
    = &(metaData.declare_field<GenericFieldType>(static_cast<stk::topology::rank_t>(metaData.side_rank()), 
                                                 "dynamic_pressure"));
  stk::mesh::put_field_on_mesh(*pBip, *part, numIntPoints, zeroVec.data());
  
  // check for total bc to create an algorithm
  if ( userData.useTotalP_ ) {
    throw std::runtime_error("LowMachFemEquationSystem::Error: use total pressure not supported: need dypP and PNG");
  }
}

//--------------------------------------------------------------------------
//-------- pre_iter_work ---------------------------------------------------
//--------------------------------------------------------------------------
void
LowMachFemEquationSystem::pre_iter_work()
{
  momentumEqSys_->pre_iter_work();
  continuityEqSys_->pre_iter_work();
}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
LowMachFemEquationSystem::solve_and_update()
{
  // wrap timing
  double timeA, timeB;
  if ( isInit_ ) {
    timeA = NaluEnv::self().nalu_time();
    // TBD: compute_dynamic_pressure();
    continuityEqSys_->compute_projected_nodal_gradient();
    copy_lagged();
    timeB = NaluEnv::self().nalu_time();
    continuityEqSys_->timerMisc_ += (timeB-timeA);
    isInit_ = false;
  }
  
  // start the iteration loop
  for ( int k = 0; k < maxIterations_; ++k ) {

    NaluEnv::self().naluOutputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << userSuppliedName_ << std::endl;

    // momentum assemble, load_complete and solve
    momentumEqSys_->assemble_and_solve(momentumEqSys_->uTmp_);

    // update all of velocity
    timeA = NaluEnv::self().nalu_time();
    field_axpby(
      realm_.meta_data(),
      realm_.bulk_data(),
      1.0, *momentumEqSys_->uTmp_,
      1.0, momentumEqSys_->velocity_->field_of_state(stk::mesh::StateNP1),
      realm_.get_activate_aura());
    timeB = NaluEnv::self().nalu_time();
    momentumEqSys_->timerAssemble_ += (timeB-timeA);

    // compute velocity relative to mesh with new velocity
    realm_.compute_vrtm();

    // continuity assemble, load_complete and solve
    continuityEqSys_->assemble_and_solve(continuityEqSys_->pTmp_);
    
    // update pressure
    timeA = NaluEnv::self().nalu_time();
    field_axpby(
      realm_.meta_data(),
      realm_.bulk_data(),
      1.0, *continuityEqSys_->pTmp_,
      1.0, *continuityEqSys_->pressure_,
      realm_.get_activate_aura());
    timeB = NaluEnv::self().nalu_time();
    continuityEqSys_->timerAssemble_ += (timeB-timeA);
  
    // copy dpdx and velocity before we update for usage in rho*ujHat
    copy_lagged();
    
    // project nodal velocity
    timeA = NaluEnv::self().nalu_time();
    project_nodal_velocity();
    timeB = NaluEnv::self().nalu_time();
    timerMisc_ += (timeB-timeA);

    // compute velocity relative to mesh with new velocity
    realm_.compute_vrtm();
  }

  // process CFL/Reynolds
  momentumEqSys_->cflReyAlgDriver_->execute();
}

//--------------------------------------------------------------------------
//-------- copy_lagged -----------------------------------------------------
//--------------------------------------------------------------------------
void
LowMachFemEquationSystem::copy_lagged()
{
  // copy velocity and projected nodal pressure gradient to a lagged set of fields
  const std::string vrtmName = realm_.does_mesh_move() ? "velocity_rtm" : "velocity";
  VectorFieldType *vrtm = realm_.meta_data().get_field<VectorFieldType>(stk::topology::NODE_RANK, vrtmName);
  field_copy(realm_.meta_data(), realm_.bulk_data(), *continuityEqSys_->dpdx_, *dpdxL_, realm_.get_activate_aura());
  field_copy(realm_.meta_data(), realm_.bulk_data(), *vrtm, *vrtmL_, realm_.get_activate_aura());
}

//--------------------------------------------------------------------------
//-------- project_nodal_velocity ------------------------------------------
//--------------------------------------------------------------------------
void
LowMachFemEquationSystem::project_nodal_velocity()
{
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // time step
  const double dt = realm_.get_time_step();
  const double gamma1 = realm_.get_gamma1();
  const double projTimeScale = dt/gamma1;
  
  const int nDim = meta_data.spatial_dimension();

  // field that we need
  VectorFieldType *velocity = momentumEqSys_->velocity_;
  VectorFieldType &velocityNp1 = velocity->field_of_state(stk::mesh::StateNP1);
  VectorFieldType *uTmp = momentumEqSys_->uTmp_;
  VectorFieldType *dpdx = continuityEqSys_->dpdx_;
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  //==========================================================
  // save off dpdx to uTmp (do it everywhere)
  //==========================================================
 
  // selector (everywhere dpdx lives) and node_buckets 
  stk::mesh::Selector s_nodes = stk::mesh::selectField(*dpdx);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * ut = stk::mesh::field_data(*uTmp, b);
    double * dp = stk::mesh::field_data(*dpdx, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j ) {
        ut[offSet+j] = dp[offSet+j];
      }
    }
  }

  //==========================================================
  // safe to update pressure gradient
  //==========================================================
  continuityEqSys_->compute_projected_nodal_gradient();

  //==========================================================
  // project u, u^n+1 = u^k+1 - dt/rho*(Gjp^N+1 - uTmp);
  //==========================================================
  
  // selector and node_buckets (only projected nodes)
  stk::mesh::Selector s_projected_nodes
    = (!stk::mesh::selectUnion(momentumEqSys_->notProjectedPart_)) &
    stk::mesh::selectField(*dpdx);
  stk::mesh::BucketVector const& p_node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_projected_nodes );
  
  // process loop
  for ( stk::mesh::BucketVector::const_iterator ib = p_node_buckets.begin() ;
        ib != p_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * uNp1 = stk::mesh::field_data(velocityNp1, b);
    double * ut = stk::mesh::field_data(*uTmp, b);
    double * dp = stk::mesh::field_data(*dpdx, b);
    double * rho = stk::mesh::field_data(densityNp1, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      // Get scaling factor
      const double fac = projTimeScale/rho[k];
      
      // projection step
      const size_t offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j ) {
        const double gdpx = dp[offSet+j] - ut[offSet+j];
        uNp1[offSet+j] -= fac*gdpx;
      }
    }
  }
}

//==========================================================================
// Class Definition
//==========================================================================
// MomentumFemEquationSystem - manages uvw pde system
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MomentumFemEquationSystem::MomentumFemEquationSystem(
  EquationSystems& eqSystems)
  : EquationSystem(eqSystems, "MomentumEQS","momentum"),
    velocity_(NULL),
    uTmp_(NULL),
    density_(NULL),
    viscosity_(NULL),
    cflReyAlgDriver_(new AlgorithmDriver(realm_))
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("velocity");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_MOMENTUM);
  linsys_ = LinearSystem::create(realm_, realm_.spatialDimension_, this, solver);

  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  // do not support PNG for velocity
  if ( realm_.get_consistent_mass_matrix_png("velocity") )
    throw std::runtime_error("MomentumFemEquationSystem::Error: This Eqs does not support a PNG equation system");
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
MomentumFemEquationSystem::~MomentumFemEquationSystem()
{
  delete cflReyAlgDriver_;
}

//--------------------------------------------------------------------------
//-------- initial_work ----------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumFemEquationSystem::initial_work()
{
  // call base class method (BDF2 state management, etc)
  EquationSystem::initial_work();

  // proceed with a bunch of initial work; wrap in timer
  const double timeA = NaluEnv::self().nalu_time();
  realm_.compute_vrtm();
  cflReyAlgDriver_->execute();

  const double timeB = NaluEnv::self().nalu_time();
  timerMisc_ += (timeB-timeA);
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
MomentumFemEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  // register dof; set it as a restart variable
  velocity_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity", numStates));
  stk::mesh::put_field_on_mesh(*velocity_, *part, nDim, nullptr);
  realm_.augment_restart_variable_list("velocity");

  // delta solution for linear solver
  uTmp_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "uTmp"));
  stk::mesh::put_field_on_mesh(*uTmp_, *part, nDim, nullptr);

  VectorFieldType *coordinates =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates"));
  stk::mesh::put_field_on_mesh(*coordinates, *part, nDim, nullptr);

  density_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "density", numStates));
  stk::mesh::put_field_on_mesh(*density_, *part, nullptr);
 
  viscosity_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity"));
  stk::mesh::put_field_on_mesh(*viscosity_, *part, nullptr);

  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && (!realm_.restarted_simulation() || realm_.support_inconsistent_restart()) ) {
    VectorFieldType &velocityN = velocity_->field_of_state(stk::mesh::StateN);
    VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
    
    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               &velocityNp1, &velocityN,
                               0, nDim,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlg);
  }
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
MomentumFemEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{
  // types of algorithms
  const AlgorithmType algType = INTERIOR;
    
  stk::topology partTopo = part->topology();
  auto& solverAlgMap = solverAlgDriver_->solverAlgorithmMap_;
    
  AssembleElemSolverAlgorithm* solverAlg = nullptr;
  bool solverAlgWasBuilt = false;
    
  std::tie(solverAlg, solverAlgWasBuilt) = build_or_add_part_to_solver_alg(*this, *part, solverAlgMap);
    
  ElemDataRequests& dataPreReqs = solverAlg->dataNeededByKernels_;
  auto& activeKernels = solverAlg->activeKernels_;
  
  if (solverAlgWasBuilt) {

    build_fem_topo_kernel_if_requested<MomentumAdvFemKernel>
      (partTopo, *this, activeKernels, "advection",
       realm_.bulk_data(), *realm_.solutionOptions_, velocity_, dataPreReqs);
   
    build_fem_topo_kernel_if_requested<MomentumDiffFemKernel>
      (partTopo, *this, activeKernels, "diffusion",
       realm_.bulk_data(), *realm_.solutionOptions_, velocity_, viscosity_, dataPreReqs);
    
    build_fem_topo_kernel_if_requested<MomentumMassFemKernel>
      (partTopo, *this, activeKernels, "momentum_time_derivative",
       realm_.bulk_data(), *realm_.solutionOptions_, velocity_, density_, dataPreReqs);

    build_fem_topo_kernel_if_requested<MomentumBodyForceFemKernel>
      (partTopo, *this, activeKernels, "body_force",
       realm_.bulk_data(), *realm_.solutionOptions_, dataPreReqs);
    
    report_invalid_supp_alg_names();
    report_built_supp_alg_names();
  }
  
  // non-solver CFL alg
  std::map<AlgorithmType, Algorithm *>::iterator it
    = cflReyAlgDriver_->algMap_.find(algType);
  if ( it == cflReyAlgDriver_->algMap_.end() ) {
    AssembleCourantReynoldsFemAlgorithm*theAlg
      = new AssembleCourantReynoldsFemAlgorithm(realm_, part);
    cflReyAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
MomentumFemEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*parTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{
  // push mesh part
  notProjectedPart_.push_back(part);

  // set face/element required sort
  realm_.solutionOptions_->set_consolidated_bc_solver_alg();
  
  // algorithm type
  const AlgorithmType algType = INFLOW;

  // velocity np1
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  
  stk::mesh::MetaData &meta_data = realm_.meta_data();
  const unsigned nDim = meta_data.spatial_dimension();

  // register boundary data; velocity_bc
  VectorFieldType *theBcField = &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_bc"));
  stk::mesh::put_field_on_mesh(*theBcField, *part, nDim, nullptr);
  
  // extract the value for user specified velocity and save off the AuxFunction
  InflowUserData userData = inflowBCData.userData_;
  std::string velocityName = "velocity";
  UserDataType theDataType = get_bc_data_type(userData, velocityName);

  AuxFunction *theAuxFunc = NULL;
  if ( CONSTANT_UD == theDataType ) {
    Velocity ux = userData.u_;
    std::vector<double> userSpec(nDim);
    userSpec[0] = ux.ux_;
    userSpec[1] = ux.uy_;
    if ( nDim > 2)
      userSpec[2] = ux.uz_;
    
    // new it
    theAuxFunc = new ConstantAuxFunction(0, nDim, userSpec);
    
  }
  else if ( FUNCTION_UD == theDataType ) {
    // extract the name/params
    std::string fcnName = get_bc_function_name(userData, velocityName);
    std::vector<double> theParams = get_bc_function_params(userData, velocityName);

    // switch on the name found...
    if ( fcnName == "pulse" ) {
      theAuxFunc = new PulseVelocityAuxFunction(0,nDim,theParams);
    }
    else {
      throw std::runtime_error("MomentumFemEquationSystem::register_inflow_bc: limited functions supported");
    }
  }
  else {
    throw std::runtime_error("MomentumFemEquationSystem::register_inflow_bc: only constant and user function supported");
  }
  
  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
			       theBcField, theAuxFunc,
			       stk::topology::NODE_RANK);
  
  // how to populate the field?
  if ( userData.externalData_ ) {
    // xfer will handle population; only need to populate the initial value
    realm_.initCondAlg_.push_back(auxAlg);
  }
  else {
    // put it on bcData
    bcDataAlg_.push_back(auxAlg);
  }

  // copy velocity_bc to velocity np1...
  CopyFieldAlgorithm *theCopyAlg
    = new CopyFieldAlgorithm(realm_, part,
                             theBcField, &velocityNp1,
                             0, nDim,
                             stk::topology::NODE_RANK);
  bcDataMapAlg_.push_back(theCopyAlg);
  
  // Dirichlet bc
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
    DirichletBC *theAlg
      = new DirichletBC(realm_, this, part, &velocityNp1, theBcField, 0, nDim);
    solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
  }
  else {
    itd->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_open_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumFemEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &partTopo,
  const OpenBoundaryConditionData &openBCData)
{
  // algorithm type
  const AlgorithmType algType = OPEN;

  // set face/element required sort
  realm_.solutionOptions_->set_consolidated_bc_solver_alg();

  // extract types of data
  OpenUserData userData = openBCData.userData_;
  stk::mesh::MetaData &metaData = realm_.meta_data();
  const int nDim = metaData.spatial_dimension();
  
  // register boundary data; open_velocity_bc
  VectorFieldType *theBcField 
    = &(metaData.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "open_velocity_bc"));
  stk::mesh::put_field_on_mesh(*theBcField, *part, nDim, nullptr);

  // extract the value for user specified velocity and save off the AuxFunction
  Velocity ux = userData.u_;
  std::vector<double> userSpec(nDim);
  userSpec[0] = ux.ux_;
  userSpec[1] = ux.uy_;
  if ( nDim > 2)
    userSpec[2] = ux.uz_;

  // new it
  ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, nDim, userSpec);

  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               theBcField, theAuxFunc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlg);

  // solver for momentum open
  auto& solverAlgMap = solverAlgDriver_->solverAlgorithmMap_;
  
  stk::topology elemTopo = get_elem_topo(realm_, *part);
  
  AssembleFaceElemSolverAlgorithm* faceElemSolverAlg = nullptr;
  bool solverAlgWasBuilt = false;

  std::tie(faceElemSolverAlg, solverAlgWasBuilt) 
    = build_or_add_part_to_face_elem_solver_alg(algType, *this, *part, elemTopo, solverAlgMap, "open");
  
  auto& activeKernels = faceElemSolverAlg->activeKernels_;
  
  if (solverAlgWasBuilt) {
    build_fem_face_elem_topo_kernel_automatic<MomentumOpenAdvDiffFemKernel>
      (partTopo, elemTopo, *this, activeKernels, "momentum_open",
       realm_.meta_data(), *realm_.solutionOptions_, velocity_, viscosity_,
       faceElemSolverAlg->faceDataNeeded_, faceElemSolverAlg->elemDataNeeded_);
  }
}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumFemEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &partTopo,
  const WallBoundaryConditionData &wallBCData)
{
  const AlgorithmType algType = WALL;
  
  // set face/element required sort
  realm_.solutionOptions_->set_consolidated_bc_solver_alg();
  
  // find out if this is a wall function approach
  WallUserData userData = wallBCData.userData_;
  const bool wallFunctionApproach = userData.wallFunctionApproach_;
  const bool anyWallFunctionActivated = wallFunctionApproach;

  // push mesh part
  if ( !anyWallFunctionActivated )
    notProjectedPart_.push_back(part);
  else 
    throw std::runtime_error("MomentumFemEquationSystem::register_wall_bc does not support wall functions");

  // np1 velocity
  VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  
  stk::mesh::MetaData &meta_data = realm_.meta_data();
  const unsigned nDim = meta_data.spatial_dimension();
  
  const std::string bcFieldName = anyWallFunctionActivated ? "wall_velocity_bc" : "velocity_bc";
  
  // register boundary data; velocity_bc
  VectorFieldType *theBcField = &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, bcFieldName));
  stk::mesh::put_field_on_mesh(*theBcField, *part, nDim, nullptr);
  
  // extract the value for user specified velocity and save off the AuxFunction
  AuxFunction *theAuxFunc = NULL;
  std::string velocityName = "velocity";

  if ( bc_data_specified(userData, velocityName) ) {

    UserDataType theDataType = get_bc_data_type(userData, velocityName);
    if ( CONSTANT_UD == theDataType ) {
      // constant data type specification
      Velocity ux = userData.u_;
      std::vector<double> userSpec(nDim);
      userSpec[0] = ux.ux_;
      userSpec[1] = ux.uy_;
      if ( nDim > 2)
        userSpec[2] = ux.uz_;
      theAuxFunc = new ConstantAuxFunction(0, nDim, userSpec);
    }
    else if ( FUNCTION_UD == theDataType ) {
      throw std::runtime_error("MomentumFemEquationSystem::register_wall_bc Zero user function wall support");
    }
  }
  else {
    throw std::runtime_error("Invalid Wall Data Specification; must provide const or fcn for velocity");
  }

  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               theBcField, theAuxFunc,
                               stk::topology::NODE_RANK);
  
  // check to see if this is an FSI interface to determine how we handle velocity population
  if ( userData.isFsiInterface_ ) {
    // xfer will handle population; only need to populate the initial value
    realm_.initCondAlg_.push_back(auxAlg);
  }
  else {
    bcDataAlg_.push_back(auxAlg);
  }
  
  // copy velocity_bc to velocity np1
  CopyFieldAlgorithm *theCopyAlg
    = new CopyFieldAlgorithm(realm_, part,
			     theBcField, &velocityNp1,
			     0, nDim,
			     stk::topology::NODE_RANK);
  
  // wall function activity will only set dof velocity np1 wall value as an IC
  if ( anyWallFunctionActivated )
    realm_.initCondAlg_.push_back(theCopyAlg);
  else
    bcDataMapAlg_.push_back(theCopyAlg);
    
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
    DirichletBC *theAlg
      = new DirichletBC(realm_, this, part, &velocityNp1, theBcField, 0, nDim);
    solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
  }
  else {
    itd->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_overset_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
MomentumFemEquationSystem::register_overset_bc()
{
  create_constraint_algorithm(velocity_);

  int nDim = realm_.meta_data().spatial_dimension();
  UpdateOversetFringeAlgorithmDriver* theAlg = new UpdateOversetFringeAlgorithmDriver(realm_);
  // Perform fringe updates before all equation system solves
  equationSystems_.preIterAlgDriver_.push_back(theAlg);
  theAlg->fields_.push_back(std::unique_ptr<OversetFieldData>(new OversetFieldData(velocity_,1,nDim)));
  
  if ( realm_.has_mesh_motion() ) {
    UpdateOversetFringeAlgorithmDriver* theAlgPost = new UpdateOversetFringeAlgorithmDriver(realm_,false);
    // Perform fringe updates after all equation system solves (ideally on the post_time_step)
    equationSystems_.postIterAlgDriver_.push_back(theAlgPost);
    theAlgPost->fields_.push_back(std::unique_ptr<OversetFieldData>(new OversetFieldData(velocity_,1,nDim)));
  }
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumFemEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
MomentumFemEquationSystem::reinitialize_linear_system()
{
  // delete linsys
  delete linsys_;

  // delete old solver
  const EquationType theEqID = EQ_MOMENTUM;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("velocity");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_MOMENTUM);
  linsys_ = LinearSystem::create(realm_, realm_.spatialDimension_, this, solver);

  // initialize new solver
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- predict_state ---------------------------------------------------
//--------------------------------------------------------------------------
void
MomentumFemEquationSystem::predict_state()
{
  // copy state n to state np1
  VectorFieldType &uN = velocity_->field_of_state(stk::mesh::StateN);
  VectorFieldType &uNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), uN, uNp1, realm_.get_activate_aura());
}

//--------------------------------------------------------------------------
//-------- register_initial_condition_fcn ----------------------------------
//--------------------------------------------------------------------------
void
MomentumFemEquationSystem::register_initial_condition_fcn(
  stk::mesh::Part *part,
  const std::map<std::string, std::string> &theNames,
  const std::map<std::string, std::vector<double> > &theParams)
{
  // extract nDim
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();

  // iterate map and check for name
  const std::string dofName = "velocity";
  std::map<std::string, std::string>::const_iterator iterName
    = theNames.find(dofName);
  if (iterName != theNames.end()) {
    std::string fcnName = (*iterName).second;
    
    // save off the field (np1 state)
    VectorFieldType *velocityNp1 = meta_data.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
    
    // create a few Aux things
    AuxFunction *theAuxFunc = NULL;
    AuxFunctionAlgorithm *auxAlg = NULL;
    std::vector<double> fcnParams;

    // extract the params
    std::map<std::string, std::vector<double> >::const_iterator iterParams
      = theParams.find(dofName);
    if (iterParams != theParams.end()) {
      fcnParams = (*iterParams).second;	
    }
    
    // create the aux function
    if ( fcnName == "TaylorGreen"  ) {
      theAuxFunc = new TaylorGreenVelocityAuxFunction(0,nDim); 
    }
    else if ( fcnName == "wind_energy_taylor_vortex") {
      // create the function
      theAuxFunc = new WindEnergyTaylorVortexAuxFunction(0,nDim,fcnParams);
    }
    else if ( fcnName == "OneTwoTenVelocity" ) {      
      theAuxFunc = new OneTwoTenVelocityAuxFunction(0,nDim);
    }
    else {
      throw std::runtime_error("InitialCondFunction::non-supported velocity IC"); 
    }

    // create the algorithm
    auxAlg = new AuxFunctionAlgorithm(realm_, part,
                                      velocityNp1, theAuxFunc,
                                      stk::topology::NODE_RANK);
    
    // push to ic
    realm_.initCondAlg_.push_back(auxAlg);
  }
}

//==========================================================================
// Class Definition
//==========================================================================
// ContinuityFemEquationSystem - manages p pde system
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ContinuityFemEquationSystem::ContinuityFemEquationSystem(
  EquationSystems& eqSystems)
  : EquationSystem(eqSystems, "ContinuityEQS","continuity"),
    managePNG_(realm_.get_consistent_mass_matrix_png("pressure")),
    pressure_(NULL),
    dpdx_(NULL),
    density_(NULL),
    pTmp_(NULL),
    projectedNodalGradEqs_(NULL)
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("pressure");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_CONTINUITY);
  linsys_ = LinearSystem::create(realm_, 1, this, solver);

  // push back EQ to manager
  realm_.equationSystems_.equationSystemVector_.push_back(this);
  
  // create projected nodal gradient equation system
  if ( managePNG_ ) {
    manage_projected_nodal_gradient(eqSystems);
  }
  else {
    throw std::runtime_error("ContinuityFemEquationSystem::Error: This Eqs does not support a PNG equation system");
  }
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ContinuityFemEquationSystem::~ContinuityFemEquationSystem()
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
ContinuityFemEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // register dof; set it as a restart variable
  pressure_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure"));
  stk::mesh::put_field_on_mesh(*pressure_, *part, nullptr);
  realm_.augment_restart_variable_list("pressure");

  dpdx_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dpdx"));
  stk::mesh::put_field_on_mesh(*dpdx_, *part, nDim, nullptr);
  
  const int numStates = realm_.number_of_states();
  density_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "density", numStates));
  stk::mesh::put_field_on_mesh(*density_, *part, nullptr);

  // delta solution for linear solver; share delta with other split systems
  pTmp_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "pTmp"));
  stk::mesh::put_field_on_mesh(*pTmp_, *part, nullptr);
 
  VectorFieldType *coordinates =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates"));
  stk::mesh::put_field_on_mesh(*coordinates, *part, nDim, nullptr);
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
ContinuityFemEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{
  stk::topology partTopo = part->topology();
  auto& solverAlgMap = solverAlgDriver_->solverAlgorithmMap_;
  
  AssembleElemSolverAlgorithm* solverAlg = nullptr;
  bool solverAlgWasBuilt = false;
    
  std::tie(solverAlg, solverAlgWasBuilt) = build_or_add_part_to_solver_alg(*this, *part, solverAlgMap);
    
  ElemDataRequests& dataPreReqs = solverAlg->dataNeededByKernels_;
  auto& activeKernels = solverAlg->activeKernels_;
  
  if (solverAlgWasBuilt) {

    build_fem_topo_kernel_if_requested<ContinuityAdvFemKernel>
      (partTopo, *this, activeKernels, "advection",
       realm_.bulk_data(), *realm_.solutionOptions_, dataPreReqs);
    
    report_invalid_supp_alg_names();
    report_built_supp_alg_names();
  }
}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
ContinuityFemEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &partTopo,
  const InflowBoundaryConditionData &inflowBCData)
{
  // set face/element required sort
  realm_.solutionOptions_->set_consolidated_bc_solver_alg();
  
  stk::mesh::MetaData &meta_data = realm_.meta_data();
  const unsigned nDim = meta_data.spatial_dimension();

  // register boundary data; cont_velocity_bc
  VectorFieldType *theBcField = &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "cont_velocity_bc"));
  stk::mesh::put_field_on_mesh(*theBcField, *part, nDim, nullptr);
  
  // extract the value for user specified velocity and save off the AuxFunction
  InflowUserData userData = inflowBCData.userData_;
  std::string velocityName = "velocity";
  UserDataType theDataType = get_bc_data_type(userData, velocityName);
  
  AuxFunction *theAuxFunc = NULL;
  if ( CONSTANT_UD == theDataType ) {
    Velocity ux = userData.u_;
    std::vector<double> userSpec(nDim);
    userSpec[0] = ux.ux_;
    userSpec[1] = ux.uy_;
    if ( nDim > 2)
      userSpec[2] = ux.uz_;
    
    // new it
    theAuxFunc = new ConstantAuxFunction(0, nDim, userSpec);    
  }
  else if ( FUNCTION_UD == theDataType ) {
    // extract the name/params
    std::string fcnName = get_bc_function_name(userData, velocityName);
    std::vector<double> theParams = get_bc_function_params(userData, velocityName);
    
    // switch on the name found...
    if ( fcnName == "pulse" ) {
      theAuxFunc = new PulseVelocityAuxFunction(0,nDim,theParams);
    }
    else {
      throw std::runtime_error("ContinuityFemEquationSystem::register_inflow_bc: limited functions supported");
    }
  }
  else {
    throw std::runtime_error("ContinuityFemEquationSystem::register_inflow_bc: only constant and user function supported");
  }
  
  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               theBcField, theAuxFunc,
                               stk::topology::NODE_RANK);
  
  // how to populate the field?
  if ( userData.externalData_ ) {
    // xfer will handle population; only need to populate the initial value
    realm_.initCondAlg_.push_back(auxAlg);
  }
  else {
    // put it on bcData
    bcDataAlg_.push_back(auxAlg);
  }
  
  // element-based uses consolidated approach fully
  auto& solverAlgMap = solverAlgDriver_->solverAlgorithmMap_;
  
  AssembleElemSolverAlgorithm* solverAlg = nullptr;
  bool solverAlgWasBuilt = false;
  
  std::tie(solverAlg, solverAlgWasBuilt) = build_or_add_part_to_face_bc_solver_alg(*this, *part, solverAlgMap, "inflow");
  
  ElemDataRequests& dataPreReqs = solverAlg->dataNeededByKernels_;
  auto& activeKernels = solverAlg->activeKernels_;
  
  if (solverAlgWasBuilt) {
    build_fem_face_topo_kernel_automatic<ContinuityInflowFemKernel>
      (partTopo, *this, activeKernels, "continuity_inflow",
       realm_.bulk_data(), *realm_.solutionOptions_, realm_.solutionOptions_->cvfemShiftMdot_, dataPreReqs);
  }
}

//--------------------------------------------------------------------------
//-------- register_open_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
ContinuityFemEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &partTopo,
  const OpenBoundaryConditionData &openBCData)
{
  // algorithm type
  const AlgorithmType algType = OPEN;
  
  // set face/element required sort
  realm_.solutionOptions_->set_consolidated_bc_solver_alg();
      
  // solver for continuity open
  auto& solverAlgMap = solverAlgDriver_->solverAlgorithmMap_;
  
  stk::topology elemTopo = get_elem_topo(realm_, *part);
  
  AssembleFaceElemSolverAlgorithm* faceElemSolverAlg = nullptr;
  bool solverAlgWasBuilt = false;

  std::tie(faceElemSolverAlg, solverAlgWasBuilt) 
    = build_or_add_part_to_face_elem_solver_alg(algType, *this, *part, elemTopo, solverAlgMap, "open");
  
  auto& activeKernels = faceElemSolverAlg->activeKernels_;
  
  if (solverAlgWasBuilt) {
    build_fem_face_elem_topo_kernel_automatic<ContinuityOpenFemKernel>
      (partTopo, elemTopo, *this, activeKernels, "continuity_open",
       realm_.meta_data(), *realm_.solutionOptions_,
       faceElemSolverAlg->faceDataNeeded_, faceElemSolverAlg->elemDataNeeded_);
  }
}

//--------------------------------------------------------------------------
//-------- register_overset_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
ContinuityFemEquationSystem::register_overset_bc()
{
  create_constraint_algorithm(pressure_);

  UpdateOversetFringeAlgorithmDriver* theAlg = new UpdateOversetFringeAlgorithmDriver(realm_);
  // Perform fringe updates before all equation system solves
  equationSystems_.preIterAlgDriver_.push_back(theAlg);

  // manage pressure; variable density requires a pre-timestep evaluation of independent variables
  theAlg->fields_.push_back(
    std::unique_ptr<OversetFieldData>(new OversetFieldData(pressure_,1,1)));
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
ContinuityFemEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
ContinuityFemEquationSystem::reinitialize_linear_system()
{
  // delete linsys
  delete linsys_;

  // delete old solver
  const EquationType theEqID = EQ_CONTINUITY;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("pressure");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_CONTINUITY);
  linsys_ = LinearSystem::create(realm_, 1, this, solver);

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- register_initial_condition_fcn ----------------------------------
//--------------------------------------------------------------------------
void
ContinuityFemEquationSystem::register_initial_condition_fcn(
  stk::mesh::Part *part,
  const std::map<std::string, std::string> &theNames,
  const std::map<std::string, std::vector<double> >& theParams)
{
  // iterate map and check for name
  const std::string dofName = "pressure";
  std::map<std::string, std::string>::const_iterator iterName
    = theNames.find(dofName);
  if (iterName != theNames.end()) {
    std::string fcnName = (*iterName).second;
    AuxFunction *theAuxFunc = NULL;
    if ( fcnName == "TaylorGreen") {
      // create the function
      theAuxFunc = new TaylorGreenPressureAuxFunction();      
    }
    else {
      throw std::runtime_error("ContinuityFemEquationSystem::register_initial_condition_fcn: limited functions supported");
    }
    
    // create the algorithm
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
				 pressure_, theAuxFunc,
				 stk::topology::NODE_RANK);
    
    // push to ic
    realm_.initCondAlg_.push_back(auxAlg);
  }
}

//--------------------------------------------------------------------------
//-------- manage_projected_nodal_gradient ---------------------------------
//--------------------------------------------------------------------------
void
ContinuityFemEquationSystem::manage_projected_nodal_gradient(
  EquationSystems& eqSystems)
{
  if ( NULL == projectedNodalGradEqs_ ) {
    projectedNodalGradEqs_ 
      = new ProjectedNodalGradientEquationSystem(eqSystems, EQ_PNG_P, "dpdx", "qTmp", "pressure", "PNGradPEQS", false, true);
  }
  // fill the map for expected boundary condition names...
  projectedNodalGradEqs_->set_data_map(INFLOW_BC, "pressure");
  projectedNodalGradEqs_->set_data_map(WALL_BC, "pressure");
  projectedNodalGradEqs_->set_data_map(OPEN_BC, "pressure_bc");
  projectedNodalGradEqs_->set_data_map(SYMMETRY_BC, "pressure");
}

//--------------------------------------------------------------------------
//-------- compute_projected_nodal_gradient---------------------------------
//--------------------------------------------------------------------------
void
ContinuityFemEquationSystem::compute_projected_nodal_gradient()
{
  if ( NULL != projectedNodalGradEqs_ )
    projectedNodalGradEqs_->solve_and_update_external();
}

} // namespace nalu
} // namespace Sierra
