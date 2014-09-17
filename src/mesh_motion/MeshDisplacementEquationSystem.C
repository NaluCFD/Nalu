/*------------------------------------------------------------------------*/
/*  Nalu 1.0 Copyright 2014 Sandia Corporation.                           */
/*  This software is released under the BSD license detailed              */
/*  in the file, LICENSE which is located in the top-level Nalu           */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <mesh_motion/MeshDisplacementEquationSystem.h>

#include <mesh_motion/AssembleMeshDisplacementElemSolverAlgorithm.h>
#include <mesh_motion/MeshDisplacementMassBackwardEulerNodeSuppAlg.h>

#include <user_functions/LinearRampMeshDisplacementAuxFunction.h>
#include <user_functions/SinMeshDisplacementAuxFunction.h>

#include <AlgorithmDriver.h>
#include <AssembleNodeSolverAlgorithm.h>
#include <AuxFunctionAlgorithm.h>
#include <ConstantAuxFunction.h>
#include <CopyFieldAlgorithm.h>
#include <DirichletBC.h>
#include <Enums.h>
#include <EquationSystem.h>
#include <EquationSystems.h>
#include <FieldFunctions.h>
#include <LinearSolver.h>
#include <LinearSolvers.h>
#include <LinearSystem.h>
#include <master_element/MasterElement.h>
#include <NaluParsing.h>
#include <Realm.h>
#include <Realms.h>
#include <Simulation.h>
#include <SolutionOptions.h>
#include <SolverAlgorithmDriver.h>
#include <TimeIntegrator.h>

// stk_util
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/Env.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/Comm.hpp>

// stk_io
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>

// stk_topo
#include <stk_topology/topology.hpp>

// basic c++
#include <iostream>
#include <math.h>
#include <utility>
#include <vector>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// MeshDisplacementEquationSystem - manages uvw pde system
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MeshDisplacementEquationSystem::MeshDisplacementEquationSystem(
  EquationSystems& eqSystems,
  const bool activateMass)
  : EquationSystem(eqSystems, "MeshDisplacementEQS"),
    activateMass_(activateMass),
    isInit_(false),
    meshDisplacement_(NULL),
    meshVelocity_(NULL),
    coordinates_(NULL),
    currentCoordinates_(NULL),
    dualNodalVolume_(NULL),
    density_(NULL),
    lameMu_(NULL),
    lameLambda_(NULL),
    dxTmp_(NULL)
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("mesh_displacement");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_MESH_DISPLACEMENT);
  linsys_ = LinearSystem::create(realm_, realm_.spatialDimension_, name_, solver);

  // push back EQ to manager
  realm_.equationSystems_.push_back(this);

  realm_.meshDeformation_ = true;

}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
MeshDisplacementEquationSystem::~MeshDisplacementEquationSystem()
{
  // does nothing
}

//--------------------------------------------------------------------------
//-------- initial_work ----------------------------------------------------
//--------------------------------------------------------------------------
void
MeshDisplacementEquationSystem::initial_work()
{
  // call base class method
  EquationSystem::initial_work();

}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
MeshDisplacementEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.fixture_->meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  // register dof; set it as a restart variable
  meshDisplacement_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "mesh_displacement", numStates));
  stk::mesh::put_field(*meshDisplacement_, *part, nDim);
  realm_.augment_restart_variable_list("mesh_displacement");

  // mesh velocity (used for fluids coupling)
  meshVelocity_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "mesh_velocity"));
  stk::mesh::put_field(*meshVelocity_, *part, nDim);

  // delta solution for linear solver
  dxTmp_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dxTmp"));
  stk::mesh::put_field(*dxTmp_, *part, nDim);

  // geometry
  coordinates_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates"));
  stk::mesh::put_field(*coordinates_, *part, nDim);

  currentCoordinates_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "current_coordinates"));
  stk::mesh::put_field(*currentCoordinates_, *part, nDim);

  dualNodalVolume_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume"));
  stk::mesh::put_field(*dualNodalVolume_, *part);

  // properties
  if ( activateMass_ ) {
    density_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "density"));
    stk::mesh::put_field(*density_, *part);
    realm_.augment_property_map(DENSITY_ID, density_);
  }

  lameMu_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "lame_mu"));
  stk::mesh::put_field(*lameMu_, *part);

  lameLambda_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "lame_lambda"));
  stk::mesh::put_field(*lameLambda_, *part);

  // push to property list
  realm_.augment_property_map(LAME_MU_ID, lameLambda_);
  realm_.augment_property_map(LAME_LAMBDA_ID, lameMu_);

  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && !realm_.restarted_simulation() ) {
    VectorFieldType &meshDisplacementN = meshDisplacement_->field_of_state(stk::mesh::StateN);
    VectorFieldType &meshDisplacementNp1 = meshDisplacement_->field_of_state(stk::mesh::StateNP1);

    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               &meshDisplacementNp1, &meshDisplacementN,
                               0, nDim,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlg);
  }

}


//--------------------------------------------------------------------------
//-------- register_element_fields -------------------------------------------
//--------------------------------------------------------------------------
void
MeshDisplacementEquationSystem::register_element_fields(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{
  //====================================================
  // Register element data
  //====================================================

  stk::io::StkMeshIoBroker *fixture = realm_.fixture_;
  stk::mesh::MetaData & meta_data = fixture->meta_data();

  const int numScvIp = theTopo.num_nodes();
  GenericFieldType *scVolume = &(meta_data.declare_field<GenericFieldType>(stk::topology::ELEMENT_RANK, "sc_volume"));
  stk::mesh::put_field(*scVolume, *part, numScvIp );

}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
MeshDisplacementEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{

  // types of algorithms
  const AlgorithmType algType = INTERIOR;
  const AlgorithmType algMass = MASS;

  // solver; interior elem contribution (advection + diffusion)
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
    solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    AssembleMeshDisplacementElemSolverAlgorithm *theAlg
      = new AssembleMeshDisplacementElemSolverAlgorithm(realm_, part, this);
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;
  }
  else {
    itsi->second->partVec_.push_back(part);
  }

  // solver; time contribution (lumped mass matrix)
  if ( activateMass_ ) {
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsm
      = solverAlgDriver_->solverAlgMap_.find(algMass);
    if ( itsm == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleNodeSolverAlgorithm *theAlg
        = new AssembleNodeSolverAlgorithm(realm_, part, this);
      solverAlgDriver_->solverAlgMap_[algMass] = theAlg;

      // now create the supplemental alg for mass term; generalize for BDF2
      MeshDisplacementMassBackwardEulerNodeSuppAlg *theMass
        = new MeshDisplacementMassBackwardEulerNodeSuppAlg(realm_);
      theAlg->supplementalAlg_.push_back(theMass);
    }
    else {
      itsm->second->partVec_.push_back(part);
    }
  }

}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
MeshDisplacementEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const WallBoundaryConditionData &wallBCData)
{

  // algorithm type
  const AlgorithmType algType = WALL;

  // np1 mesh displacement
  VectorFieldType &displacementNp1 = meshDisplacement_->field_of_state(stk::mesh::StateNP1);

  stk::mesh::MetaData &meta_data = realm_.fixture_->meta_data();
  const unsigned nDim = meta_data.spatial_dimension();

  // register boundary data; mesh_displacement_bc
  VectorFieldType *theBcField = &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "mesh_displacement_bc"));
  stk::mesh::put_field(*theBcField, *part, nDim);

  // extract the value for user specified velocity and save off the AuxFunction
  WallUserData userData = wallBCData.userData_;
  std::string displacementName = "mesh_displacement";

  if ( bc_data_specified(userData, displacementName) ) {

    AuxFunction *theAuxFunc = NULL;

    UserDataType theDataType = get_bc_data_type(userData, displacementName);
    if ( CONSTANT_UD == theDataType ) {
      // constant data type specification
      Velocity dx = userData.dx_;
      std::vector<double> userSpec(nDim);
      userSpec[0] = dx.ux_;
      userSpec[1] = dx.uy_;
      if ( nDim > 2)
        userSpec[2] = dx.uz_;
      theAuxFunc = new ConstantAuxFunction(0, nDim, userSpec);
    }
    else if ( FUNCTION_UD == theDataType ) {
      // extract the name
      std::string fcnName = get_bc_function_name(userData, displacementName);
      std::vector<double> theParams = get_bc_function_params(userData, displacementName);
      // switch on the name found...
      if ( fcnName == "linear" ) {
        theAuxFunc = new LinearRampMeshDisplacementAuxFunction(0,nDim, theParams);
      }
      else if ( fcnName == "sinusoidal") {
        theAuxFunc = new SinMeshDisplacementAuxFunction(0,nDim, theParams);
      }
      else {
        throw std::runtime_error("Only linear ramp user functions supported");
      }
    }

    // proceed with aux function and dirichlet setup

    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
          theBcField, theAuxFunc,
          stk::topology::NODE_RANK);
    bcDataAlg_.push_back(auxAlg);


    // copy mesh_displacement_bc to mesh_displacement np1...
    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
          theBcField, &displacementNp1,
          0, nDim,
          stk::topology::NODE_RANK);
    bcDataMapAlg_.push_back(theCopyAlg);

    // Dirichlet bc
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itd
      = solverAlgDriver_->solverDirichAlgMap_.find(algType);
    if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
      DirichletBC *theAlg
        = new DirichletBC(realm_, this, part, &displacementNp1, theBcField, 0, nDim);
        solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
    }
    else {
      itd->second->partVec_.push_back(part);
    }
  }
  else {
    Env::outputP0() << "No displacement specified: zero surface traction applied" << std::endl;
  }

}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------
//--------------------------------------------------------------------------
void
MeshDisplacementEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
MeshDisplacementEquationSystem::reinitialize_linear_system()
{

  // delete linsys
  delete linsys_;

  // delete old solver
  const EquationType theEqID = EQ_MESH_DISPLACEMENT;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("mesh_velocity");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_MESH_DISPLACEMENT);
  linsys_ = LinearSystem::create(realm_, realm_.spatialDimension_, name_, solver);

  // initialize new solver
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}


//--------------------------------------------------------------------------
//-------- predict_state ---------------------------------------------------
//--------------------------------------------------------------------------
void
MeshDisplacementEquationSystem::predict_state()
{
  stk::mesh::MetaData & meta_data = realm_.fixture_->meta_data();

  const int nDim = meta_data.spatial_dimension();

  VectorFieldType &displacementN = meshDisplacement_->field_of_state(stk::mesh::StateN);
  VectorFieldType &displacementNp1 = meshDisplacement_->field_of_state(stk::mesh::StateNP1);

  // define some common selectors; select all nodes (locally and shared)
  // where velocity is defined
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*meshDisplacement_);

  //===========================================================
  // copy state N into N+1
  //===========================================================

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * dxN = stk::mesh::field_data(displacementN, b);
    double * dxNp1 = stk::mesh::field_data(displacementNp1, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j ) {
        dxNp1[offSet+j] = dxN[offSet+j];
      }
    }
  }

}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
MeshDisplacementEquationSystem::solve_and_update()
{

  // nothing to do
  if ( isInit_ ) {
    isInit_ = false;
  }

  // start the iteration loop
  for ( int k = 0; k < maxIterations_; ++k ) {

    Env::outputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << name_ << std::endl;

    // tke assemble, load_complete and solve
    assemble_and_solve(dxTmp_);

    // update
    double timeA = stk::cpu_time();
    field_axpby(
      realm_.fixture_->meta_data(),
      realm_.fixture_->bulk_data(),
      1.0, *dxTmp_,
      1.0, meshDisplacement_->field_of_state(stk::mesh::StateNP1));
    double timeB = stk::cpu_time();
    timerAssemble_ += (timeB-timeA);

    compute_current_coordinates();
  }

}


//--------------------------------------------------------------------------
//-------- compute_current_coordinates -------------------------------------
//--------------------------------------------------------------------------
void
MeshDisplacementEquationSystem::compute_current_coordinates()
{
  stk::mesh::MetaData & meta_data = realm_.fixture_->meta_data();

  VectorFieldType &displacementN = meshDisplacement_->field_of_state(stk::mesh::StateN);
  VectorFieldType &displacementNp1 = meshDisplacement_->field_of_state(stk::mesh::StateNP1);

  const int nDim = meta_data.spatial_dimension();
  const double dt = realm_.get_time_step();

  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
      &stk::mesh::selectField(*meshDisplacement_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    const double * dxN = stk::mesh::field_data(displacementN, b);
    const double * dxNp1 = stk::mesh::field_data(displacementNp1, b);
    double * meshVelocity = stk::mesh::field_data(*meshVelocity_, b);
    const double * coordinates = stk::mesh::field_data(*coordinates_, b);
    double * currentCoordinates = stk::mesh::field_data(*currentCoordinates_, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      size_t offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j ) {
        currentCoordinates[offSet+j] = coordinates[offSet+j] + dxNp1[offSet+j];
        // hack a mesh velocity to be first order backward Euler
        meshVelocity[offSet+j] = (dxNp1[offSet+j] - dxN[offSet+j])/dt;
      }
    }
  }

}

} // namespace nalu
} // namespace Sierra
