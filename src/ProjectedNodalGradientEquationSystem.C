/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "ProjectedNodalGradientEquationSystem.h"

#include "AssemblePNGElemSolverAlgorithm.h"
#include "AssemblePNGBoundarySolverAlgorithm.h"
#include "AssemblePNGPressureBoundarySolverAlgorithm.h"
#include "AssemblePNGNonConformalSolverAlgorithm.h"
#include "EquationSystem.h"
#include "EquationSystems.h"
#include "Enums.h"
#include "FieldFunctions.h"
#include "LinearSolvers.h"
#include "LinearSolver.h"
#include "LinearSystem.h"
#include "NaluEnv.h"
#include "Realm.h"
#include "Realms.h"
#include "Simulation.h"
#include "SolutionOptions.h"
#include "SolverAlgorithmDriver.h"

// user functions
#include "user_functions/SteadyThermalContactAuxFunction.h"

//overset
#include "overset/UpdateOversetFringeAlgorithmDriver.h"

// template for kernels
#include "AlgTraits.h"
#include "kernel/KernelBuilder.h"
#include "kernel/KernelBuilderLog.h"

// kernels
#include "AssembleElemSolverAlgorithm.h"
#include "kernel/ScalarPngFemKernel.h"

// bc kernels
#include "kernel/ScalarPngBcFemKernel.h"

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
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/Comm.hpp>

// stk_io
#include <stk_io/IossBridge.hpp>

#include <stk_topology/topology.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// ProjectedNodalGradientEquationSystem - do some stuff
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ProjectedNodalGradientEquationSystem::ProjectedNodalGradientEquationSystem(
 EquationSystems& eqSystems,
 const EquationType eqType,
 const std::string dofName, 
 const std::string deltaName, 
 const std::string independentDofName,
 const std::string eqSysName,
 const bool managesSolve,
 const bool isFEM)
  : EquationSystem(eqSystems, eqSysName, dofName),
    eqType_(eqType),
    dofName_(dofName),
    deltaName_(deltaName),
    independentDofName_(independentDofName),
    eqSysName_(eqSysName),
    managesSolve_(managesSolve),
    isFEM_(isFEM),
    dqdx_(NULL),
    qTmp_(NULL)
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name(dofName);
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, eqType_);
  linsys_ = LinearSystem::create(realm_, realm_.spatialDimension_, this, solver);

  // push back EQ to manager
  realm_.push_equation_to_systems(this);
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ProjectedNodalGradientEquationSystem::~ProjectedNodalGradientEquationSystem()
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- set_data_map ----------------------------------------------------
//--------------------------------------------------------------------------
void
ProjectedNodalGradientEquationSystem::set_data_map( 
  BoundaryConditionType BC, std::string name)
{
  dataMap_[BC] = name;
}

//--------------------------------------------------------------------------
//-------- get_name_given_bc -----------------------------------------------
//--------------------------------------------------------------------------
std::string
ProjectedNodalGradientEquationSystem::get_name_given_bc( 
  BoundaryConditionType BC)
{
  std::map<BoundaryConditionType, std::string>::iterator it;
  it=dataMap_.find(BC);
  if ( it == dataMap_.end() )
    throw std::runtime_error("PNGEqSys::missing BC type specification (developer error)!");
  else
    return it->second;
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
ProjectedNodalGradientEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  dqdx_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, dofName_));
  stk::mesh::put_field_on_mesh(*dqdx_, *part, nDim, nullptr);

  // delta solution for linear solver
  qTmp_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, deltaName_));
  stk::mesh::put_field_on_mesh(*qTmp_, *part, nDim, nullptr);
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
ProjectedNodalGradientEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{
  // type of algorithms
  const AlgorithmType algType = INTERIOR;

  if (!isFEM_ ) {
    // solver
    std::map<AlgorithmType, SolverAlgorithm *>::iterator its
      = solverAlgDriver_->solverAlgMap_.find(algType);
    if ( its == solverAlgDriver_->solverAlgMap_.end() ) {
      AssemblePNGElemSolverAlgorithm *theAlg
        = new AssemblePNGElemSolverAlgorithm(realm_, part, this, independentDofName_, dofName_);
      solverAlgDriver_->solverAlgMap_[algType] = theAlg;
    }
    else {
      its->second->partVec_.push_back(part);
    }
  }
  else {
    // consolidated only
    stk::topology partTopo = part->topology();
    auto& solverAlgMap = solverAlgDriver_->solverAlgorithmMap_;
    
    AssembleElemSolverAlgorithm* solverAlg = nullptr;
    bool solverAlgWasBuilt = false;
    
    std::tie(solverAlg, solverAlgWasBuilt) = build_or_add_part_to_solver_alg(*this, *part, solverAlgMap);
    
    ElemDataRequests& dataPreReqs = solverAlg->dataNeededByKernels_;
    auto& activeKernels = solverAlg->activeKernels_;
    
    if (solverAlgWasBuilt) {
      
      build_fem_topo_kernel_if_requested<ScalarPngFemKernel>
        (partTopo, *this, activeKernels, "interior_png",
         realm_.bulk_data(), *realm_.solutionOptions_, independentDofName_, dofName_, dataPreReqs);
      
      report_invalid_supp_alg_names();
      report_built_supp_alg_names();
    } 
  }

}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
ProjectedNodalGradientEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &partTopo,
  const WallBoundaryConditionData &/*wallBCData*/)
{
  const AlgorithmType algType = WALL;
  
  // extract the field name for this bc type
  std::string fieldName = get_name_given_bc(WALL_BC);

  if ( !isFEM_ ) {        
    // create lhs/rhs [CVFEM] algorithm;
    std::map<AlgorithmType, SolverAlgorithm *>::iterator its =
      solverAlgDriver_->solverAlgMap_.find(algType);
    if ( its == solverAlgDriver_->solverAlgMap_.end() ) {
      AssemblePNGBoundarySolverAlgorithm *theAlg
        = new AssemblePNGBoundarySolverAlgorithm(realm_, part, this, fieldName);
      solverAlgDriver_->solverAlgMap_[algType] = theAlg;
    }
    else {
      its->second->partVec_.push_back(part);
    }
  }
  else {
    // element-based uses consolidated approach fully
    auto& solverAlgMap = solverAlgDriver_->solverAlgorithmMap_;
    
    AssembleElemSolverAlgorithm* solverAlg = nullptr;
    bool solverAlgWasBuilt = false;
    
    std::tie(solverAlg, solverAlgWasBuilt) = build_or_add_part_to_face_bc_solver_alg(*this, *part, solverAlgMap, "wall_png");
        
    ElemDataRequests& dataPreReqs = solverAlg->dataNeededByKernels_;
    auto& activeKernels = solverAlg->activeKernels_;
    
    if (solverAlgWasBuilt) {
      build_fem_face_topo_kernel_automatic<ScalarPngBcFemKernel>
        (partTopo, *this, activeKernels, "png_wall",
         realm_.bulk_data(), *realm_.solutionOptions_, fieldName, dataPreReqs);
      report_built_supp_alg_names();   
    }
  }
}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
ProjectedNodalGradientEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &partTopo,
  const InflowBoundaryConditionData &/*inflowBCData*/)
{
  // type of algorithm
  const AlgorithmType algType = INFLOW;
  
  // extract the field name for this bc type
  std::string fieldName = get_name_given_bc(INFLOW_BC);
  
  if ( !isFEM_ ) {
    // create lhs/rhs algorithm;
    std::map<AlgorithmType, SolverAlgorithm *>::iterator its =
      solverAlgDriver_->solverAlgMap_.find(algType);
    if ( its == solverAlgDriver_->solverAlgMap_.end() ) {
      AssemblePNGBoundarySolverAlgorithm *theAlg
        = new AssemblePNGBoundarySolverAlgorithm(realm_, part, this, fieldName);
      solverAlgDriver_->solverAlgMap_[algType] = theAlg;
    }
    else {
      its->second->partVec_.push_back(part);
    }
  }
  else {
    // element-based uses consolidated approach fully
    auto& solverAlgMap = solverAlgDriver_->solverAlgorithmMap_;
    
    AssembleElemSolverAlgorithm* solverAlg = nullptr;
    bool solverAlgWasBuilt = false;
    
    std::tie(solverAlg, solverAlgWasBuilt) = build_or_add_part_to_face_bc_solver_alg(*this, *part, solverAlgMap, "inflow_png");
    
    ElemDataRequests& dataPreReqs = solverAlg->dataNeededByKernels_;
    auto& activeKernels = solverAlg->activeKernels_;
    
    if (solverAlgWasBuilt) {
      build_fem_face_topo_kernel_automatic<ScalarPngBcFemKernel>
        (partTopo, *this, activeKernels, "png_inflow",
         realm_.bulk_data(), *realm_.solutionOptions_, fieldName, dataPreReqs);
      report_built_supp_alg_names();   
    }  
  }
}

//--------------------------------------------------------------------------
//-------- register_open_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
ProjectedNodalGradientEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &partTopo,
  const OpenBoundaryConditionData &/*openBCData*/)
{
  // type of algorithm
  const AlgorithmType algType = OPEN;

  // extract the field name for this bc type
  std::string fieldName = get_name_given_bc(OPEN_BC);

  if ( !isFEM_ ) {
    // create lhs/rhs algorithm;
    std::map<AlgorithmType, SolverAlgorithm *>::iterator its =
      solverAlgDriver_->solverAlgMap_.find(algType);
    if ( its == solverAlgDriver_->solverAlgMap_.end() ) {
      SolverAlgorithm *theSolverAlg = NULL;
      if ( fieldName == "pressure_bc" ) {
        // Gjp requires a special algorithm
        theSolverAlg = new AssemblePNGPressureBoundarySolverAlgorithm(realm_, part, this, fieldName);
      }
      else {
        theSolverAlg = new AssemblePNGBoundarySolverAlgorithm(realm_, part, this, fieldName);
      }
      solverAlgDriver_->solverAlgMap_[algType] = theSolverAlg;
    }
    else {
      its->second->partVec_.push_back(part);
    }
  }
  else {
    // element-based uses consolidated approach fully
    auto& solverAlgMap = solverAlgDriver_->solverAlgorithmMap_;
    
    AssembleElemSolverAlgorithm* solverAlg = nullptr;
    bool solverAlgWasBuilt = false;
    
    std::tie(solverAlg, solverAlgWasBuilt) = build_or_add_part_to_face_bc_solver_alg(*this, *part, solverAlgMap, "open_png");
    
    ElemDataRequests& dataPreReqs = solverAlg->dataNeededByKernels_;
    auto& activeKernels = solverAlg->activeKernels_;
    
    if (solverAlgWasBuilt) {
      build_fem_face_topo_kernel_automatic<ScalarPngBcFemKernel>
        (partTopo, *this, activeKernels, "png_open",
         realm_.bulk_data(), *realm_.solutionOptions_, fieldName, dataPreReqs);
      report_built_supp_alg_names();   
    }  
  }

}

//--------------------------------------------------------------------------
//-------- register_symmetry_bc --------------------------------------------
//--------------------------------------------------------------------------
void
ProjectedNodalGradientEquationSystem::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const SymmetryBoundaryConditionData &/*symmetryBCData*/)
{
  // not supported yet for FEM
  if ( isFEM_ ) 
    return;
  
  const AlgorithmType algType = SYMMETRY;

  // extract the field name for this bc type
  std::string fieldName = get_name_given_bc(SYMMETRY_BC);
  // create lhs/rhs algorithm;
  std::map<AlgorithmType, SolverAlgorithm *>::iterator its =
    solverAlgDriver_->solverAlgMap_.find(algType);
  if ( its == solverAlgDriver_->solverAlgMap_.end() ) {
    AssemblePNGBoundarySolverAlgorithm *theAlg
      = new AssemblePNGBoundarySolverAlgorithm(realm_, part, this, fieldName);
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;
  }
  else {
    its->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_non_conformal_bc ---------------------------------------
//--------------------------------------------------------------------------
void
ProjectedNodalGradientEquationSystem::register_non_conformal_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/)
{
  // not supported yet for FEM
  if ( isFEM_ ) 
    return;

  const AlgorithmType algType = NON_CONFORMAL;

  // create lhs/rhs algorithm;
  std::map<AlgorithmType, SolverAlgorithm *>::iterator its =
    solverAlgDriver_->solverAlgMap_.find(algType);
  if ( its == solverAlgDriver_->solverAlgMap_.end() ) {
    AssemblePNGNonConformalSolverAlgorithm *theAlg
      = new AssemblePNGNonConformalSolverAlgorithm(realm_, part, this, independentDofName_, dofName_, realm_.solutionOptions_->ncAlgPngPenalty_);
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;
  }
  else {
    its->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_overset_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
ProjectedNodalGradientEquationSystem::register_overset_bc()
{
  create_constraint_algorithm(dqdx_);

  const int nDim = realm_.meta_data().spatial_dimension();
  UpdateOversetFringeAlgorithmDriver* theAlg = new UpdateOversetFringeAlgorithmDriver(realm_);
  // Perform fringe updates before all equation system solves
  equationSystems_.preIterAlgDriver_.push_back(theAlg);
  theAlg->fields_.push_back(
    std::unique_ptr<OversetFieldData>(new OversetFieldData(dqdx_,1,nDim)));

  if ( realm_.has_mesh_motion() ) {
    UpdateOversetFringeAlgorithmDriver* theAlgPost = new UpdateOversetFringeAlgorithmDriver(realm_,false);
    // Perform fringe updates after all equation system solves (ideally on the post_time_step)
    equationSystems_.postIterAlgDriver_.push_back(theAlgPost);
    theAlgPost->fields_.push_back(std::unique_ptr<OversetFieldData>(new OversetFieldData(dqdx_,1,nDim)));
  }
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
ProjectedNodalGradientEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
ProjectedNodalGradientEquationSystem::reinitialize_linear_system()
{
  // delete linsys; set previously set parameters on linsys
  const bool provideOutput = linsys_->provideOutput_;
  delete linsys_;

  // delete old solver
  const EquationType theEqID = eqType_;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver; reset parameters
  std::string solverName = realm_.equationSystems_.get_solver_block_name(dofName_);
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, eqType_);
  linsys_ = LinearSystem::create(realm_, realm_.spatialDimension_, this, solver);
  linsys_->provideOutput_ = provideOutput;

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
ProjectedNodalGradientEquationSystem::solve_and_update()
{
  if ( managesSolve_ )
    solve_and_update_external();
}

//--------------------------------------------------------------------------
//-------- solve_and_update_external ---------------------------------------
//--------------------------------------------------------------------------
void
ProjectedNodalGradientEquationSystem::solve_and_update_external()
{
  for ( int k = 0; k < maxIterations_; ++k ) {

    // projected nodal gradient, load_complete and solve
    assemble_and_solve(qTmp_);
    
    // update
    double timeA = NaluEnv::self().nalu_time();
    field_axpby(
      realm_.meta_data(),
      realm_.bulk_data(),
      1.0, *qTmp_,
      1.0, *dqdx_, 
      realm_.get_activate_aura());
    double timeB = NaluEnv::self().nalu_time();
    timerAssemble_ += (timeB-timeA);   
  }
}

//--------------------------------------------------------------------------
//-------- deactivate_output -----------------------------------------------
//--------------------------------------------------------------------------
void
ProjectedNodalGradientEquationSystem::deactivate_output()
{
  linsys_->provideOutput_ = false;
}

} // namespace nalu
} // namespace Sierra
