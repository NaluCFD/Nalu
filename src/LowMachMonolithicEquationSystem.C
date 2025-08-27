/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <LowMachMonolithicEquationSystem.h>
#include <AlgorithmDriver.h>
#include <CopyFieldAlgorithm.h>
#include <Enums.h>
#include <EquationSystem.h>
#include <EquationSystems.h>
#include <FieldFunctions.h>
#include <LinearSolver.h>
#include <LinearSolvers.h>
#include <LinearSystem.h>
#include <master_element/MasterElement.h>
#include <NaluEnv.h>
#include <NaluParsing.h>
#include <Realm.h>
#include <Realms.h>
#include <Simulation.h>
#include <SolutionOptions.h>
#include <SolverAlgorithmDriver.h>

// aux algorithms
#include <AuxFunctionAlgorithm.h>

// algorithms
#include <AssembleCourantReynoldsElemAlgorithm.h>
#include <AssembleNodalGradAlgorithmDriver.h>
#include <AssembleNodalGradUAlgorithmDriver.h>
#include <AssembleNodalGradElemAlgorithm.h>
#include "AssembleNodalGradUElemAlgorithm.h"

// user function
#include <user_functions/TaylorGreenPressureAuxFunction.h>
#include <user_functions/TaylorGreenVelocityAuxFunction.h>

// template for kernels
#include <AlgTraits.h>
#include <kernel/KernelBuilder.h>
#include <kernel/KernelBuilderLog.h>

// kernels
#include <kernel/MomentumContinuityElemKernel.h>

// stk_util
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/Comm.hpp>

// stk_io
#include <stk_io/IossBridge.hpp>

// stk_topo
#include <stk_topology/topology.hpp>

// basic c++
#include <vector>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// LowMachMonolithicEquationSystem - manage the low Mach equation system (uvw_p)
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
LowMachMonolithicEquationSystem::LowMachMonolithicEquationSystem(
  EquationSystems& eqSystems)
  : EquationSystem(eqSystems, "LowMachMonoEOS","uvwp"),
    velocity_(NULL),
    pressure_(NULL),
    dudx_(NULL),
    dpdx_(NULL),
    dpdxOld_(NULL),
    uvwp_(NULL),
    uvwpTmp_(NULL),
    coordinates_(NULL),
    density_(NULL),
    viscosity_(NULL),
    assembleNodalGradUAlgDriver_(new AssembleNodalGradUAlgorithmDriver(realm_, "dudx")),
    assembleNodalGradPAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "pressure", "dpdx")),
    cflReyAlgDriver_(new AlgorithmDriver(realm_)),
    nDim_(realm_.meta_data().spatial_dimension()),
    sizeOfSystem_(realm_.meta_data().spatial_dimension()+1),
    edgeNodalGradientU_(false),
    edgeNodalGradientP_(false),
    isInit_(true)
{

  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("uvwp");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_UVWP);
  linsys_ = LinearSystem::create(realm_, sizeOfSystem_, this, solver);

  // determine nodal gradient form; velocity
  std::map<std::string, std::string >::iterator iiU =
    realm_.solutionOptions_->nodalGradMap_.find("velocity");
  if ( iiU != realm_.solutionOptions_->nodalGradMap_.end() ) {
    if ( iiU->second == "edge" ) 
      edgeNodalGradientU_ = true;
  }

  // pressure
  std::map<std::string, std::string >::iterator iiP =
    realm_.solutionOptions_->nodalGradMap_.find("pressure");
  if ( iiP != realm_.solutionOptions_->nodalGradMap_.end() ) {
    if ( iiP->second == "edge" ) 
      edgeNodalGradientP_ = true;
  }

  // output to user
  NaluEnv::self().naluOutputP0() << "LowMachMonoEOS::Edge projected nodal gradient for velocity/pressure: " 
                                 << edgeNodalGradientU_ << "/" 
                                 << edgeNodalGradientP_ <<std::endl;
  
  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  // inform realm
  realm_.hasFluids_ = true;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
LowMachMonolithicEquationSystem::~LowMachMonolithicEquationSystem()
{
  delete assembleNodalGradPAlgDriver_;
  delete assembleNodalGradUAlgDriver_;
  delete cflReyAlgDriver_;
}

//--------------------------------------------------------------------------
//-------- initial_work ----------------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::initial_work()
{
  // call base class method (BDF2 state management, etc)
  EquationSystem::initial_work();

  // proceed with a bunch of initial work; wrap in timer
  const double timeA = NaluEnv::self().nalu_time();
  if ( !realm_.restarted_simulation() ) {
    copy_uvw_p_to_uvwp();
  }

  realm_.compute_vrtm();
  assembleNodalGradPAlgDriver_->execute();
  assembleNodalGradUAlgDriver_->execute();
  cflReyAlgDriver_->execute();

  const double timeB = NaluEnv::self().nalu_time();
  timerMisc_ += (timeB-timeA);
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();

  const int numStates = realm_.number_of_states();
  
  // dof and delta solution for linear solver; set as restart
  uvwp_ =  &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "uvwp", numStates));
  stk::mesh::put_field_on_mesh(*uvwp_, *part, sizeOfSystem_, nullptr);
  uvwpTmp_ =  &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "uvwpTmp"));
  stk::mesh::put_field_on_mesh(*uvwpTmp_, *part, sizeOfSystem_, nullptr);
  realm_.augment_restart_variable_list("uvwp");

  // register convenience dof
  velocity_ =  &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "velocity", numStates));
  stk::mesh::put_field_on_mesh(*velocity_, *part, nDim_, nullptr);
  stk::io::set_field_output_type(*velocity_, stk::io::FieldOutputType::VECTOR_3D);

  dudx_ =  &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "dudx"));
  stk::mesh::put_field_on_mesh(*dudx_, *part, nDim*nDim, nullptr);

  pressure_ =  &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "pressure"));
  stk::mesh::put_field_on_mesh(*pressure_, *part, nullptr);
  dpdx_ =  &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "dpdx"));
  stk::mesh::put_field_on_mesh(*dpdx_, *part, nDim_, nullptr);
  stk::io::set_field_output_type(*dpdx_, stk::io::FieldOutputType::VECTOR_3D);

  // old
  dpdxOld_ =  &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "dpdx_old"));
  stk::mesh::put_field_on_mesh(*dpdxOld_, *part, nDim_, nullptr);
  stk::io::set_field_output_type(*dpdxOld_, stk::io::FieldOutputType::VECTOR_3D);

  coordinates_ =  &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "coordinates"));
  stk::mesh::put_field_on_mesh(*coordinates_, *part, nDim, nullptr);
  stk::io::set_field_output_type(*coordinates_, stk::io::FieldOutputType::VECTOR_3D);

  // add properties; density needs to be a restart field
  density_ =  &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "density", numStates));
  stk::mesh::put_field_on_mesh(*density_, *part, nullptr);
  realm_.augment_restart_variable_list("density");

  viscosity_ =  &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "viscosity"));
  stk::mesh::put_field_on_mesh(*viscosity_, *part, nullptr);

  // push to property list
  realm_.augment_property_map(DENSITY_ID, density_);
  realm_.augment_property_map(VISCOSITY_ID, viscosity_);

  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && (!realm_.restarted_simulation() || realm_.support_inconsistent_restart()) ) {
    VectorFieldType &velocityN = velocity_->field_of_state(stk::mesh::StateN);
    VectorFieldType &velocityNp1 = velocity_->field_of_state(stk::mesh::StateNP1);
    
    CopyFieldAlgorithm *theCopyAlgU
      = new CopyFieldAlgorithm(realm_, part,
                               &velocityNp1, &velocityN,
                               0, nDim_,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlgU);

    ScalarFieldType &densityN = density_->field_of_state(stk::mesh::StateN);
    ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);
    
    CopyFieldAlgorithm *theCopyAlgRho
      = new CopyFieldAlgorithm(realm_, part,
                               &densityNp1, &densityN,
                               0, 1,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlgRho);
  }
}

//--------------------------------------------------------------------------
//-------- register_element_fields -------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::register_element_fields(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // provide mean element Peclet and Courant fields; always...
  GenericFieldType *elemReynolds
    = &(meta_data.declare_field<double>(stk::topology::ELEMENT_RANK, "element_reynolds"));
  stk::mesh::put_field_on_mesh(*elemReynolds, *part, nullptr);
  GenericFieldType *elemCourant
    = &(meta_data.declare_field<double>(stk::topology::ELEMENT_RANK, "element_courant"));
  stk::mesh::put_field_on_mesh(*elemCourant, *part, nullptr);
}

//--------------------------------------------------------------------------
//-------- register_edge_fields -------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::register_edge_fields(
  stk::mesh::Part *part)
{
  if ( realm_.realmUsesEdges_ ) {
    throw std::runtime_error("LowMachMonolithicEquationSystem::Error edge-based scheme not supported");
  }
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{
  // types of algorithms
  const AlgorithmType algType = INTERIOR;

  // non-solver CFL alg
  std::map<AlgorithmType, Algorithm *>::iterator it
    = cflReyAlgDriver_->algMap_.find(algType);
  if ( it == cflReyAlgDriver_->algMap_.end() ) {
    AssembleCourantReynoldsElemAlgorithm*theAlg
      = new AssembleCourantReynoldsElemAlgorithm(realm_, part);
    cflReyAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

  // Gjui
  std::map<AlgorithmType, Algorithm *>::iterator itGU
    = assembleNodalGradUAlgDriver_->algMap_.find(algType);
  if ( itGU == assembleNodalGradUAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg = new AssembleNodalGradUElemAlgorithm(realm_, part, velocity_, dudx_, edgeNodalGradientU_);
    assembleNodalGradUAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    itGU->second->partVec_.push_back(part);
  }

  // Gjp
  std::map<AlgorithmType, Algorithm *>::iterator itGP
    = assembleNodalGradPAlgDriver_->algMap_.find(algType);
  if ( itGP == assembleNodalGradPAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, pressure_, dpdx_, edgeNodalGradientP_);
    assembleNodalGradPAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    itGP->second->partVec_.push_back(part);
  }
  
  // linear solver
  stk::topology partTopo = part->topology();
  auto& solverAlgMap = solverAlgDriver_->solverAlgorithmMap_;
  
  AssembleElemSolverAlgorithm* solverAlg = nullptr;
  bool solverAlgWasBuilt = false;

  std::tie(solverAlg, solverAlgWasBuilt) = build_or_add_part_to_solver_alg
    (*this, *part, solverAlgMap);
  
  ElemDataRequests& dataPreReqs = solverAlg->dataNeededByKernels_;
  auto& activeKernels = solverAlg->activeKernels_;
  
  if (solverAlgWasBuilt) {
   
    build_topo_kernel_if_requested<MomentumContinuityElemKernel>
      (partTopo, *this, activeKernels, "uvwp_time_advection_diffusion",
       realm_.bulk_data(), *realm_.solutionOptions_, velocity_, density_, viscosity_, false, dataPreReqs);
    build_topo_kernel_if_requested<MomentumContinuityElemKernel>
      (partTopo, *this, activeKernels, "uvwp_lumped_time_advection_diffusion",
       realm_.bulk_data(), *realm_.solutionOptions_, velocity_, density_, viscosity_, true, dataPreReqs);
    report_invalid_supp_alg_names();
    report_built_supp_alg_names();
  }
}
  
//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{
  throw std::runtime_error("LowMachMonolithicEquationSystem::register_inflow_bc: not supported");
}
  
//--------------------------------------------------------------------------
//-------- register_open_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const OpenBoundaryConditionData &openBCData)
{
  throw std::runtime_error("LowMachMonolithicEquationSystem::register_open_bc: not supported");
}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const WallBoundaryConditionData &wallBCData)
{
  throw std::runtime_error("LowMachMonolithicEquationSystem::register_wall_bc: not supported");
}

//--------------------------------------------------------------------------
//-------- register_symmetry_bc --------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const SymmetryBoundaryConditionData &/*symmBCData*/)
{
  throw std::runtime_error("LowMachMonolithicEquationSystem::register_symmetry_bc: not supported");
}

//--------------------------------------------------------------------------
//-------- register_initial_condition_fcn ----------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::register_initial_condition_fcn(
  stk::mesh::Part *part,
  const std::map<std::string, std::string> &theNames,
  const std::map<std::string, std::vector<double> > &theParams)
{
  // extract nDim
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();

  // iterate map and check for velocity
  const std::string velocityName = "velocity";
  std::map<std::string, std::string>::const_iterator iterU
    = theNames.find(velocityName);
  if (iterU != theNames.end()) {
    std::string fcnName = (*iterU).second;
    
    // create a few Aux things
    AuxFunction *theAuxFunc = NULL;
    AuxFunctionAlgorithm *auxAlg = NULL;

    if ( fcnName == "TaylorGreen"  ) {
      theAuxFunc = new TaylorGreenVelocityAuxFunction(0,nDim); 
    }
    else {
      throw std::runtime_error("InitialCondFunction::non-supported velocity IC"); 
    }

    // create the algorithm
    auxAlg = new AuxFunctionAlgorithm(realm_, part,
                                      velocity_, theAuxFunc,
                                      stk::topology::NODE_RANK, true, realm_.hasPeriodic_);
    
    // push to ic
    realm_.initCondAlg_.push_back(auxAlg);
  }
  
  // iterate map and check for pressure
  const std::string pressureName = "pressure";
  std::map<std::string, std::string>::const_iterator iterP
    = theNames.find(pressureName);
  if (iterP != theNames.end()) {
    std::string fcnName = (*iterP).second;
    
    // create a few Aux things
    AuxFunction *theAuxFunc = NULL;
    AuxFunctionAlgorithm *auxAlg = NULL;

    if ( fcnName == "TaylorGreen"  ) {
      theAuxFunc = new TaylorGreenPressureAuxFunction(); 
    }
    else {
      throw std::runtime_error("InitialCondFunction::non-supported pressure IC"); 
    }

    // create the algorithm
    auxAlg = new AuxFunctionAlgorithm(realm_, part,
                                      pressure_, theAuxFunc,
                                      stk::topology::NODE_RANK, true, true);
    
    // push to ic
    realm_.initCondAlg_.push_back(auxAlg);
  }
}

//--------------------------------------------------------------------------
//-------- pre_timestep_work -----------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::pre_timestep_work()
{
  // call base class
  EquationSystem::pre_timestep_work();

  // predict state was based on uvwp 
  copy_uvwp_to_uvw_p();

  /* FIXME: 
     The state management, at present, is a bit confusing for the monolithic system. 
     Specifically, data for velocity and pressure are populated from, e.g., xfers,
     time varying data, etc. However, the dof that is predicted is uvwp. The above
     unload will destroy any boundary data in uvw or p. 
  */
}

//--------------------------------------------------------------------------
//-------- pre_iter_work ---------------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::pre_iter_work()
{
  // nothing yet...
}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::solve_and_update()
{
  // wrap timing
  double timeA, timeB;
  if ( isInit_ ) {
    const double timeA = -NaluEnv::self().nalu_time();
    assembleNodalGradPAlgDriver_->execute();
    assembleNodalGradUAlgDriver_->execute();
    timerMisc_ += (NaluEnv::self().nalu_time() + timeA);
    isInit_ = false;
  }
  
  // start the iteration loop
  for ( int k = 0; k < maxIterations_; ++k ) {

    NaluEnv::self().naluOutputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << userSuppliedName_ << std::endl;

    // momentum assemble, load_complete and solve
    assemble_and_solve(uvwpTmp_);

    // update
    timeA = NaluEnv::self().nalu_time();
    field_axpby(
      realm_.meta_data(),
      realm_.bulk_data(),
      1.0, *uvwpTmp_,
      1.0, *uvwp_,
      realm_.get_activate_aura());

    // copy uvwp to velocity(uvw) and p
    copy_uvwp_to_uvw_p();

    timeB = NaluEnv::self().nalu_time();
    timerAssemble_ += (timeB-timeA);

    // compute velocity relative to mesh with new velocity
    realm_.compute_vrtm();

    // compute projected nodal gradients
    assembleNodalGradPAlgDriver_->execute();
    assembleNodalGradUAlgDriver_->execute();
  }
  
  // process CFL/Reynolds (mdot omitted)
  cflReyAlgDriver_->execute();
}
  
//--------------------------------------------------------------------------
//-------- predict_state ---------------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::predict_state()
{
  // copy state n to state np1
  GenericFieldType &uvwpN = uvwp_->field_of_state(stk::mesh::StateN);
  GenericFieldType &uvwpNp1 = uvwp_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), uvwpN, uvwpNp1, realm_.get_activate_aura());
  copy_uvwp_to_uvw_p();
}

//--------------------------------------------------------------------------
//-------- post_converged_work ---------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::post_converged_work()
{
  // copy gradP to gradPOld
  field_copy(realm_.meta_data(), realm_.bulk_data(), *dpdx_, *dpdxOld_, realm_.get_activate_aura());
} 
  
//--------------------------------------------------------------------------
//-------- copy_uvwp_to_uvw_p ----------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::copy_uvwp_to_uvw_p()
{
  // selector and node_buckets
  stk::mesh::Selector s_all_nodes
    = stk::mesh::selectField(*uvwp_);
  stk::mesh::BucketVector const& p_node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes);
  
  // process loop
  for ( stk::mesh::BucketVector::const_iterator ib = p_node_buckets.begin() ;
        ib != p_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * velocity = stk::mesh::field_data(*velocity_, b);
    double * pressure = stk::mesh::field_data(*pressure_, b);
    const double * uvwp = stk::mesh::field_data(*uvwp_, b);
    
    for ( stk::mesh::Bucket::size_type k = 0; k < length; ++k ) {
      const size_t kNdim = k*nDim_;
      const size_t kSize = k*sizeOfSystem_;
      for ( int j = 0; j < nDim_; ++j ) {
        velocity[kNdim+j] = uvwp[kSize+j];
      }
      pressure[k] = uvwp[kSize+nDim_];
    }
  }
}

//--------------------------------------------------------------------------
//-------- copy_uvw_p_to_uvwp ----------------------------------------------
//--------------------------------------------------------------------------
void
LowMachMonolithicEquationSystem::copy_uvw_p_to_uvwp()
{
  // selector and node_buckets
  stk::mesh::Selector s_all_nodes
    = stk::mesh::selectField(*uvwp_);
  stk::mesh::BucketVector const& p_node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes);
  
  // process loop
  for ( stk::mesh::BucketVector::const_iterator ib = p_node_buckets.begin() ;
        ib != p_node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    const double * velocity = stk::mesh::field_data(*velocity_, b);
    const double * pressure = stk::mesh::field_data(*pressure_, b);
    double * uvwp = stk::mesh::field_data(*uvwp_, b);
    
    for ( stk::mesh::Bucket::size_type k = 0; k < length; ++k ) {
      const size_t kNdim = k*nDim_;
      const size_t kSize = k*sizeOfSystem_;
      for ( int j = 0; j < nDim_; ++j ) {
        uvwp[kSize+j] = velocity[kNdim+j];
      }
      uvwp[kSize+nDim_] = pressure[k];
    }
  }
}

} // namespace nalu
} // namespace Sierra
