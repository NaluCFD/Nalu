/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "MixtureFractionFemEquationSystem.h"
#include "AlgorithmDriver.h"
#include "AuxFunctionAlgorithm.h"
#include "AssembleNodalGradAlgorithmDriver.h"
#include "ConstantAuxFunction.h"
#include "CopyFieldAlgorithm.h"
#include "DirichletBC.h"
#include "EffectiveDiffFluxCoeffAlgorithm.h"
#include "EquationSystem.h"
#include "EquationSystems.h"
#include "Enums.h"
#include "FieldFunctions.h"
#include "LinearSolvers.h"
#include "LinearSolver.h"
#include "LinearSystem.h"
#include "NaluEnv.h"
#include "NaluParsing.h"
#include "ProjectedNodalGradientEquationSystem.h"
#include "Realm.h"
#include "Realms.h"
#include "Simulation.h"
#include "SolutionOptions.h"
#include "TimeIntegrator.h"
#include "SolverAlgorithmDriver.h"

// algorithm drivers
#include "AssembleCourantReynoldsFemAlgorithm.h"

// template for kernels
#include "AlgTraits.h"
#include "kernel/KernelBuilder.h"
#include "kernel/KernelBuilderLog.h"

// kernels
#include "AssembleElemSolverAlgorithm.h"
#include "kernel/ScalarMassFemKernel.h"
#include "kernel/ScalarAdvFemKernel.h"
#include "kernel/ScalarDiffFemKernel.h"

// kernels; src - n/a

// kernels; bc
#include "kernel/ScalarFluxPenaltyFemKernel.h"
#include "kernel/ScalarOpenAdvFemKernel.h"

// kernels; nso
#include "nso/ScalarNSOFemKernel.h"

// user function
#include "user_functions/WorkshopMMSMixFracAuxFunction.h"

// overset - n/a

// stk_util
#include <stk_util/parallel/Parallel.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>

// stk_io
#include <stk_io/IossBridge.hpp>

// stk_topo
#include <stk_topology/topology.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

// nalu utility
#include <utils/StkHelpers.h>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// MixtureFractionFemEquationSystem - manages z pde system for FEM
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MixtureFractionFemEquationSystem::MixtureFractionFemEquationSystem(
  EquationSystems& eqSystems,
  const bool outputClippingDiag,
  const double deltaZClip,
  const bool computePNG)
  : EquationSystem(eqSystems, "MixtureFractionFemEQS", "mixture_fraction"),
    outputClippingDiag_(outputClippingDiag),
    deltaZClip_(deltaZClip),
    computePNG_(computePNG),
    managePNG_(realm_.get_consistent_mass_matrix_png("mixture_fraction")),
    mixFrac_(NULL),
    mixFracUF_(NULL),
    zTmp_(NULL),
    Gjz_(NULL),
    velocity_(NULL),
    density_(NULL),
    visc_(NULL),
    tvisc_(NULL),
    evisc_(NULL),
    diffFluxCoeffAlgDriver_(new AlgorithmDriver(realm_)),
    cflReyAlgDriver_(new AlgorithmDriver(realm_)),
    projectedNodalGradEqs_(NULL),
    isInit_(true)
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("mixture_fraction");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_MIXTURE_FRACTION);
  linsys_ = LinearSystem::create(realm_, 1, this, solver);

  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  // advertise as fluids, non-uniform
  realm_.hasFluids_ = true;
  realm_.uniform_ = false;

  // create projected nodal gradient equation system
  if ( computePNG_ ) {
    if ( managePNG_ ) {
      manage_projected_nodal_gradient(eqSystems);
    }
    else {
      throw std::runtime_error("MixtureFractionFemEquationSystem must activate consistent mass matrix if compute_png is active");
    }
  }
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
MixtureFractionFemEquationSystem::~MixtureFractionFemEquationSystem()
{
  delete diffFluxCoeffAlgDriver_;
  delete cflReyAlgDriver_;
}

//--------------------------------------------------------------------------
//-------- populate_derived_quantities -------------------------------------
//--------------------------------------------------------------------------
void
MixtureFractionFemEquationSystem::populate_derived_quantities()
{
  // placeholder
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
MixtureFractionFemEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  // register dof; set it as a restart variable
  mixFrac_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixture_fraction", numStates));
  stk::mesh::put_field_on_mesh(*mixFrac_, *part, nullptr);
  realm_.augment_restart_variable_list("mixture_fraction");

  // for a sanity check, keep around the un-filterd/clipped field
  mixFracUF_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "uf_mixture_fraction", numStates));
  stk::mesh::put_field_on_mesh(*mixFracUF_, *part, nullptr);
 
  // delta solution for linear solver; share delta since this is a split system
  zTmp_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "pTmp"));
  stk::mesh::put_field_on_mesh(*zTmp_, *part, nullptr);
  
  // aux variables
  velocity_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity", numStates));
  stk::mesh::put_field_on_mesh(*velocity_, *part, nDim, nullptr);
  realm_.augment_restart_variable_list("velocity");

  // properties
  density_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "density", numStates));
  stk::mesh::put_field_on_mesh(*density_, *part, nullptr);
  realm_.augment_restart_variable_list("density");

  visc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity"));
  stk::mesh::put_field_on_mesh(*visc_, *part, nullptr);

  // push to property list
  realm_.augment_property_map(DENSITY_ID, density_);
  realm_.augment_property_map(VISCOSITY_ID, visc_);

  if ( realm_.is_turbulent() ) {
    tvisc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity"));
    stk::mesh::put_field_on_mesh(*tvisc_, *part, nullptr);
  }

  evisc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "effective_viscosity_z"));
  stk::mesh::put_field_on_mesh(*evisc_, *part, nullptr);

  // projected nodal gradients  
  ScalarFieldType *dualNodalVolume = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume"));
  stk::mesh::put_field_on_mesh(*dualNodalVolume, *part, nullptr);
  
  Gjz_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dzdx"));
  stk::mesh::put_field_on_mesh(*Gjz_, *part, nullptr);
 
  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && (!realm_.restarted_simulation() || realm_.support_inconsistent_restart()) ) {
    ScalarFieldType &mixFracN = mixFrac_->field_of_state(stk::mesh::StateN);
    ScalarFieldType &mixFracNp1 = mixFrac_->field_of_state(stk::mesh::StateNP1);

    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               &mixFracNp1, &mixFracN,
                               0, 1,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlg);

    // manage velocity as well (consequence of the stand-alone)
    VectorFieldType &velocityN = velocity_->field_of_state(stk::mesh::StateN);
    VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
    
    CopyFieldAlgorithm *theCopyAlgVel
      = new CopyFieldAlgorithm(realm_, part,
                               &velocityNp1, &velocityN,
                               0, nDim,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlgVel);
  }

  // manage density for now
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

}

//--------------------------------------------------------------------------
//-------- register_element_fields -------------------------------------------
//--------------------------------------------------------------------------
void
MixtureFractionFemEquationSystem::register_element_fields(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // provide mean element Peclet and Courant fields
  GenericFieldType *elemReynolds
    = &(meta_data.declare_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "element_reynolds"));
  stk::mesh::put_field_on_mesh(*elemReynolds, *part, 1, nullptr);
  GenericFieldType *elemCourant
    = &(meta_data.declare_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "element_courant"));
  stk::mesh::put_field_on_mesh(*elemCourant, *part, 1, nullptr);
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
MixtureFractionFemEquationSystem::register_interior_algorithm(
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

    build_fem_topo_kernel_if_requested<ScalarMassFemKernel>
      (partTopo, *this, activeKernels, "mixture_fraction_time_derivative",
       realm_.bulk_data(), *realm_.solutionOptions_, mixFrac_, density_, dataPreReqs);
    
    build_fem_topo_kernel_if_requested<ScalarAdvFemKernel>
      (partTopo, *this, activeKernels, "advection",
       realm_.bulk_data(), *realm_.solutionOptions_, mixFrac_, density_, velocity_, dataPreReqs);
    
    build_fem_topo_kernel_if_requested<ScalarDiffFemKernel>
      (partTopo, *this, activeKernels, "diffusion",
       realm_.bulk_data(), *realm_.solutionOptions_, mixFrac_, evisc_, dataPreReqs);
    
    build_fem_topo_kernel_if_requested<ScalarNSOFemKernel>
      (partTopo, *this, activeKernels, "nso",
       realm_.bulk_data(), *realm_.solutionOptions_, mixFrac_, dataPreReqs);
    
    report_invalid_supp_alg_names();
    report_built_supp_alg_names();
  }
  
  // effective viscosity alg
  const double lamSc = realm_.get_lam_schmidt(mixFrac_->name());
  const double turbSc = realm_.get_turb_schmidt(mixFrac_->name());
  std::map<AlgorithmType, Algorithm *>::iterator itev =
    diffFluxCoeffAlgDriver_->algMap_.find(algType);
  if ( itev == diffFluxCoeffAlgDriver_->algMap_.end() ) {
    EffectiveDiffFluxCoeffAlgorithm *theAlg
      = new EffectiveDiffFluxCoeffAlgorithm(realm_, part, visc_, tvisc_, evisc_, lamSc, turbSc);
    diffFluxCoeffAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    itev->second->partVec_.push_back(part);
  }

  // manage Courant number
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
MixtureFractionFemEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{
  // algorithm type
  const AlgorithmType algType = INFLOW;

  ScalarFieldType &mixFracNp1 = mixFrac_->field_of_state(stk::mesh::StateNP1);
 
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; mixFrac_bc
  ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixFrac_bc"));
  stk::mesh::put_field_on_mesh(*theBcField, *part, nullptr);

  // extract the value for user specified mixFrac and save off the AuxFunction
  InflowUserData userData = inflowBCData.userData_;
  std::string mixFracName = "mixture_fraction";
  UserDataType theDataType = get_bc_data_type(userData, mixFracName);

  AuxFunction *theAuxFunc = NULL;
  if ( CONSTANT_UD == theDataType ) {
    MixtureFraction mixFrac = userData.mixFrac_;
    std::vector<double> userSpec(1);
    userSpec[0] = mixFrac.mixFrac_;

    // new it
    theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);
  }
  else if ( FUNCTION_UD == theDataType ) {
    throw std::runtime_error("MixFracFemEquationSystem::register_inflow_bc: no inflow functions supported");
  }
  else {
    throw std::runtime_error("MixFracFemEquationSystem::register_inflow_bc: only constant and user function supported");
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

  // copy mixFrac_bc to mixture_fraction np1...
  CopyFieldAlgorithm *theCopyAlg
    = new CopyFieldAlgorithm(realm_, part,
                             theBcField, &mixFracNp1,
                             0, 1,
                             stk::topology::NODE_RANK);
  bcDataMapAlg_.push_back(theCopyAlg);

  // Dirichlet bc
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
    DirichletBC *theAlg
      = new DirichletBC(realm_, this, part, &mixFracNp1, theBcField, 0, 1);
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
MixtureFractionFemEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &partTopo,
  const OpenBoundaryConditionData &openBCData)
{
  // algorithm type
  const AlgorithmType algType = OPEN;
 
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; mixFrac_bc
  ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "open_mixFrac_bc"));
  stk::mesh::put_field_on_mesh(*theBcField, *part, nullptr);

  // extract the value for user specified mixFrac and save off the AuxFunction
  OpenUserData userData = openBCData.userData_;
  MixtureFraction mixFrac = userData.mixFrac_;
  std::vector<double> userSpec(1);
  userSpec[0] = mixFrac.mixFrac_;

  // new it
  ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               theBcField, theAuxFunc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlg);

  // solver for mixture fraction open
  auto& solverAlgMap = solverAlgDriver_->solverAlgorithmMap_;
    
  stk::topology elemTopo = get_elem_topo(realm_, *part);
  
  AssembleFaceElemSolverAlgorithm* faceElemSolverAlg = nullptr;
  bool solverAlgWasBuilt = false;
  
  std::tie(faceElemSolverAlg, solverAlgWasBuilt) 
    = build_or_add_part_to_face_elem_solver_alg(algType, *this, *part, elemTopo, solverAlgMap, "open");
    
  auto& activeKernels = faceElemSolverAlg->activeKernels_;
  
  if (solverAlgWasBuilt) {

    build_fem_face_elem_topo_kernel_automatic<ScalarOpenAdvFemKernel>
      (partTopo, elemTopo, *this, activeKernels, "mixture_fraction_open",
       realm_.meta_data(), *realm_.solutionOptions_, mixFrac_, theBcField, 
       faceElemSolverAlg->faceDataNeeded_, faceElemSolverAlg->elemDataNeeded_); 
  }
}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
MixtureFractionFemEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &partTopo,
  const WallBoundaryConditionData &wallBCData)
{

  // algorithm type
  const AlgorithmType algType = WALL;

  // np1
  ScalarFieldType &mixFracNp1 = mixFrac_->field_of_state(stk::mesh::StateNP1);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // extract the value for user specified mixFrac and save off the AuxFunction
  WallUserData userData = wallBCData.userData_;
  std::string mixFracName = "mixture_fraction";
  if ( bc_data_specified(userData, mixFracName) ) {

    // FIXME: Generalize for constant vs function

    // register boundary data; mixFrac_bc
    ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixFrac_bc"));
    stk::mesh::put_field_on_mesh(*theBcField, *part, nullptr);

    // extract data
    std::vector<double> userSpec(1);
    MixtureFraction mixFrac = userData.mixFrac_;
    userSpec[0] = mixFrac.mixFrac_;

    // new it
    ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

    // bc data alg
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 theBcField, theAuxFunc,
                                 stk::topology::NODE_RANK);
    bcDataAlg_.push_back(auxAlg);

    // copy mixFrac_bc to mixFrac np1...
    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               theBcField, &mixFracNp1,
                               0, 1,
                               stk::topology::NODE_RANK);
    bcDataMapAlg_.push_back(theCopyAlg);

    // use weak penalty approach for wall bc if consolidated bc approach is activated
    if ( realm_.solutionOptions_->useConsolidatedBcSolverAlg_ ) {

      auto& solverAlgMap = solverAlgDriver_->solverAlgorithmMap_;
      
      stk::topology elemTopo = get_elem_topo(realm_, *part);
      
      AssembleFaceElemSolverAlgorithm* faceElemSolverAlg = nullptr;
      bool solverAlgWasBuilt = false;
      
      std::tie(faceElemSolverAlg, solverAlgWasBuilt) 
        = build_or_add_part_to_face_elem_solver_alg(algType, *this, *part, elemTopo, solverAlgMap, "wall");
      
      auto& activeKernels = faceElemSolverAlg->activeKernels_;
      
      if (solverAlgWasBuilt) {
        build_fem_face_elem_topo_kernel_automatic<ScalarFluxPenaltyFemKernel>
          (partTopo, elemTopo, *this, activeKernels, "mixture_fraction_wall",
           realm_.meta_data(), *realm_.solutionOptions_, mixFrac_, theBcField, evisc_,
           faceElemSolverAlg->faceDataNeeded_, faceElemSolverAlg->elemDataNeeded_);
      }
    }
    else {
      // Dirichlet bc
      std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
        solverAlgDriver_->solverDirichAlgMap_.find(algType);
      if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
        DirichletBC *theAlg
          = new DirichletBC(realm_, this, part, &mixFracNp1, theBcField, 0, 1);
        solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
      }
      else {
        itd->second->partVec_.push_back(part);
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- register_symmetry_bc --------------------------------------------
//--------------------------------------------------------------------------
void
MixtureFractionFemEquationSystem::register_symmetry_bc(
  stk::mesh::Part */*part*/,
  const stk::topology &/*theTopo*/,
  const SymmetryBoundaryConditionData &/*wallBCData*/)
{
  // no-op
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
MixtureFractionFemEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
MixtureFractionFemEquationSystem::reinitialize_linear_system()
{

  // delete linsys
  delete linsys_;

  // delete old solver
  const EquationType theEqID = EQ_MIXTURE_FRACTION;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("mixture_fraction");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_MIXTURE_FRACTION);
  linsys_ = LinearSystem::create(realm_, 1, this, solver);

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- register_initial_condition_fcn ----------------------------------
//--------------------------------------------------------------------------
void
MixtureFractionFemEquationSystem::register_initial_condition_fcn(
  stk::mesh::Part *part,
  const std::map<std::string, std::string> &theNames,
  const std::map<std::string, std::vector<double> > &/*theParams*/)
{
  // iterate map and check for name
  const std::string dofName = "mixture_fraction";
  std::map<std::string, std::string>::const_iterator iterName
    = theNames.find(dofName);
  if (iterName != theNames.end()) {
    std::string fcnName = (*iterName).second;
    AuxFunction *theAuxFunc = NULL;
    if ( fcnName == "workshop_mms" ) {
      theAuxFunc = new WorkshopMMSMixFracAuxFunction();      
    }
    else {
      throw std::runtime_error("MixtureFractionFemEquationSystem::register_initial_condition_fcn: VariableDensity only supported");
    }
    
    // create the algorithm
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
				 mixFrac_, theAuxFunc,
				 stk::topology::NODE_RANK);
    
    // push to ic
    realm_.initCondAlg_.push_back(auxAlg);
    
  }
}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
MixtureFractionFemEquationSystem::solve_and_update()
{
  // nothing to do
  if ( isInit_ ) {
    isInit_ = false;
    compute_projected_nodal_gradient();
  }

  // compute effective viscosity
  diffFluxCoeffAlgDriver_->execute();

  //return;

  for ( int k = 0; k < maxIterations_; ++k ) {

    NaluEnv::self().naluOutputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << userSuppliedName_ << std::endl;

    // mixture fraction assemble, load_complete and solve
    assemble_and_solve(zTmp_);

    // update
    double timeA = NaluEnv::self().nalu_time();
    update_and_clip();
    double timeB = NaluEnv::self().nalu_time();
    timerAssemble_ += (timeB-timeA);
  }

  // projected nodal gradient (assumed decoupled from iteration above)
  compute_projected_nodal_gradient();
  
  cflReyAlgDriver_->execute();
}

//--------------------------------------------------------------------------
//-------- update_and_clip -------------------------------------------------
//--------------------------------------------------------------------------
void
MixtureFractionFemEquationSystem::update_and_clip()
{
  const double deltaZ = deltaZClip_;
  const double lowBound = 0.0-deltaZ;
  const double highBound = 1.0+deltaZ;
  size_t numClip[2] = {0,0};
  double minZ = +1.0e16;
  double maxZ = -1.0e16;

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*mixFrac_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    double *mixFrac = stk::mesh::field_data(*mixFrac_, b);
    double *mixFracUF = stk::mesh::field_data(*mixFracUF_, b);
    double *zTmp    = stk::mesh::field_data(*zTmp_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      double mixFracNp1 = mixFrac[k] + zTmp[k];
      // store un-filtered value for numerical methods development purposes
      mixFracUF[k] = mixFracNp1;
      // clip now
      if ( mixFracNp1 < lowBound ) {
        minZ = std::min(mixFracNp1, minZ);
        mixFracNp1 = lowBound;
        numClip[0]++;
      }
      else if ( mixFracNp1 > highBound ) {
        maxZ = std::max(mixFracNp1, maxZ);
        mixFracNp1 = highBound;
        numClip[1]++;
      }
      mixFrac[k] = mixFracNp1;
    }
  }

  // parallel assemble clipped value
  if ( outputClippingDiag_ ) {
    size_t g_numClip[2] = {};
    stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
    stk::all_reduce_sum(comm, numClip, g_numClip, 2);

    if ( g_numClip[0] > 0 ) {
      double g_minZ = 0;
      stk::all_reduce_min(comm, &minZ, &g_minZ, 1);
      NaluEnv::self().naluOutputP0() << "mixFrac clipped (-) " << g_numClip[0] << " times; min: " << g_minZ << std::endl;
    }
    else {
      NaluEnv::self().naluOutputP0() << "mixFrac clipped (-) zero times" << std::endl;
    }

    if ( g_numClip[1] > 0 ) {
      double g_maxZ = 0;
      stk::all_reduce_max(comm, &maxZ, &g_maxZ, 1);
      NaluEnv::self().naluOutputP0() << "mixFrac clipped (+) " << g_numClip[1] << " times; max: " << g_maxZ << std::endl;
    }
    else {
      NaluEnv::self().naluOutputP0() << "mixFrac clipped (+) zero times" << std::endl;
    }
  }
}

//--------------------------------------------------------------------------
//-------- predict_state ---------------------------------------------------
//--------------------------------------------------------------------------
void
MixtureFractionFemEquationSystem::predict_state()
{
  // copy state n to state np1
  ScalarFieldType &zN = mixFrac_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &zNp1 = mixFrac_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), zN, zNp1, realm_.get_activate_aura());

  // manage velocity (consequence of this stand-alone passive scalar state)
  VectorFieldType &uN = velocity_->field_of_state(stk::mesh::StateN);
  VectorFieldType &uNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), uN, uNp1, realm_.get_activate_aura());
}

//--------------------------------------------------------------------------
//-------- manage_projected_nodal_gradient ---------------------------------
//--------------------------------------------------------------------------
void
MixtureFractionFemEquationSystem::manage_projected_nodal_gradient(
  EquationSystems& eqSystems)
{
  if ( NULL == projectedNodalGradEqs_ ) {
    projectedNodalGradEqs_ 
      = new ProjectedNodalGradientEquationSystem(eqSystems, EQ_PNG_Z, "dzdx", "qTmp", "mixture_fraction", "PNGradZEQS", 
                                                 false, true);
  }
  // fill the map for expected boundary condition names...
  projectedNodalGradEqs_->set_data_map(INFLOW_BC, "mixture_fraction");
  projectedNodalGradEqs_->set_data_map(WALL_BC, "mixture_fraction");
  projectedNodalGradEqs_->set_data_map(OPEN_BC, "mixture_fraction");
  projectedNodalGradEqs_->set_data_map(SYMMETRY_BC, "mixture_fraction");
}

//--------------------------------------------------------------------------
//-------- compute_projected_nodal_gradient---------------------------------
//--------------------------------------------------------------------------
void
MixtureFractionFemEquationSystem::compute_projected_nodal_gradient()
{
  if ( NULL != projectedNodalGradEqs_ )
    projectedNodalGradEqs_->solve_and_update_external();
}

} // namespace nalu
} // namespace Sierra
