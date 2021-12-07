//*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "VolumeOfFluidEquationSystem.h"
#include "AlgorithmDriver.h"
#include "AssembleNodalGradAlgorithmDriver.h"
#include "AssembleNodalGradElemAlgorithm.h"
#include "AssembleNodalGradBoundaryAlgorithm.h"
#include "AuxFunctionAlgorithm.h"
#include "ConstantAuxFunction.h"
#include "CopyFieldAlgorithm.h"
#include "DirichletBC.h"
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

// template for kernels
#include "AlgTraits.h"
#include "kernel/KernelBuilder.h"
#include "kernel/KernelBuilderLog.h"

// kernels
#include "AssembleElemSolverAlgorithm.h"
#include "kernel/VolumeOfFluidMassElemKernel.h"
#include "kernel/VolumeOfFluidScvAdvElemKernel.h"
#include "kernel/VolumeOfFluidScsAdvElemKernel.h"
#include "kernel/VolumeOfFluidDivElemKernel.h"
#include "kernel/VolumeOfFluidSucvNsoElemKernel.h"
#include "kernel/VolumeOfFluidSharpenElemKernel.h"
#include "kernel/VolumeOfFluidOpenAdvElemKernel.h"
#include "kernel/VolumeOfFluidGclElemKernel.h"

// user function
#include "user_functions/RayleighTaylorMixFracAuxFunction.h"
#include "user_functions/FixedHeightMixFracAuxFunction.h"

#include "overset/UpdateOversetFringeAlgorithmDriver.h"

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
// VolumeOfFluidEquationSystem - manages vof pde system
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
VolumeOfFluidEquationSystem::VolumeOfFluidEquationSystem(
  EquationSystems& eqSystems,
  const bool outputClippingDiag,
  const double deltaVofClip,
  const double Fo,
  const double cAlpha,
  const bool smooth,
  const int smoothIter)
  : EquationSystem(eqSystems, "VolumeOfFluidEQS", "volume_of_fluid"),
    managePNG_(realm_.get_consistent_mass_matrix_png("volume_of_fluid")),
    vof_(NULL),
    vofSmoothed_(NULL),
    smoothedRhs_(NULL),
    interfaceNormal_(NULL),
    interfaceCurvature_(NULL),
    surfaceTension_(NULL),
    dvofdx_(NULL),
    vofTmp_(NULL),
    assembleNodalGradAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "volume_of_fluid", "dvofdx")),
    projectedNodalGradEqs_(NULL),
    outputClippingDiag_(outputClippingDiag),
    deltaVofClip_(deltaVofClip),
    Fo_(Fo),
    cAlpha_(cAlpha),
    dxMin_(1.0e16),
    smooth_(smooth),
    smoothIter_(smoothIter),
    isInit_(true),
    scsAdvection_(false),
    buoyancyStab_(realm_.solutionOptions_->buoyancyPressureStab_ ? 1.0 : 0.0)
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("volume_of_fluid");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_VOLUME_OF_FLUID);
  linsys_ = LinearSystem::create(realm_, 1, this, solver);

  // determine nodal gradient form
  set_nodal_gradient("volume_of_fluid");
  if ( realm_.realmUsesEdges_ )
    NaluEnv::self().naluOutputP0() << "VolumeOfFluidEquationSystem will use CVFEM" << std::endl;
      
  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  // advertise as non uniform
  realm_.uniform_ = false;

  // create projected nodal gradient equation system
  if ( managePNG_ ) {
    manage_projected_nodal_gradient(eqSystems);
  }

  // look for scs_advection (controls an additional open bc)
  std::map<std::string, std::vector<std::string> >::iterator isrc
    = realm_.solutionOptions_->elemSrcTermsMap_.find("volume_of_fluid");
  if ( isrc != realm_.solutionOptions_->elemSrcTermsMap_.end() ) {
    std::vector<std::string> mapNameVec = isrc->second;
    for (size_t k = 0; k < mapNameVec.size(); ++k ) {
      std::string sourceName = mapNameVec[k];
      if ( sourceName == "scs_advection" ) {
        scsAdvection_ = true;
        break;
      }
    }
  }

  // output options
  if ( smooth_ ) {
    NaluEnv::self().naluOutputP0() << "SCS smoothing_iterations: " << smoothIter_ << std::endl;
    NaluEnv::self().naluOutputP0() << "fourier_number: " << Fo_ << std::endl;
  }

  if ( supp_alg_is_requested("sharpen") )
    NaluEnv::self().naluOutputP0() << "compression_constant: " << cAlpha_ << std::endl;

  if ( scsAdvection_ ) 
    NaluEnv::self().naluOutputP0() << "VOF scs_advection is active " << std::endl;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
VolumeOfFluidEquationSystem::~VolumeOfFluidEquationSystem()
{
  delete assembleNodalGradAlgDriver_;
}

//--------------------------------------------------------------------------
//-------- populate_derived_quantities -------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::populate_derived_quantities()
{
  // placeholder
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  // register dof; set it as a restart variable
  vof_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "volume_of_fluid", numStates));
  stk::mesh::put_field_on_mesh(*vof_, *part, nullptr);
  realm_.augment_restart_variable_list("volume_of_fluid");

  // smoothed field and smoothed rhs
  if ( smooth_ ) {
    vofSmoothed_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "volume_of_fluid_smoothed"));
    stk::mesh::put_field_on_mesh(*vofSmoothed_, *part, nullptr);    
    smoothedRhs_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "smoothed_rhs"));
    stk::mesh::put_field_on_mesh(*smoothedRhs_, *part, nullptr);
  }

  // normal
  interfaceNormal_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "interface_normal"));
  stk::mesh::put_field_on_mesh(*interfaceNormal_, *part, nDim, nullptr);

  // always register curvature and surface tension; let property specification dictate its importance
  interfaceCurvature_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "interface_curvature"));
  stk::mesh::put_field_on_mesh(*interfaceCurvature_, *part, nullptr);
  surfaceTension_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "surface_tension"));
  stk::mesh::put_field_on_mesh(*surfaceTension_, *part, nullptr);

  // push to property list
  realm_.augment_property_map(SURFACE_TENSION_ID, surfaceTension_);

  // projected nodal gradient
  dvofdx_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dvofdx"));
  stk::mesh::put_field_on_mesh(*dvofdx_, *part, nDim, nullptr);

  // delta solution for linear solver; share delta since this is a split system
  vofTmp_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "pTmp"));
  stk::mesh::put_field_on_mesh(*vofTmp_, *part, nullptr);

  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && (!realm_.restarted_simulation() || realm_.support_inconsistent_restart()) ) {
    ScalarFieldType &vofN = vof_->field_of_state(stk::mesh::StateN);
    ScalarFieldType &vofNp1 = vof_->field_of_state(stk::mesh::StateNP1);

    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               &vofNp1, &vofN,
                               0, 1,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlg);
  }
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{

  // types of algorithms
  const AlgorithmType algType = INTERIOR;

  ScalarFieldType &vofNp1 = (smooth_) ? vofSmoothed_->field_of_state(stk::mesh::StateNone) 
    : vof_->field_of_state(stk::mesh::StateNP1);  
  VectorFieldType &dvofdxNone = dvofdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to projected nodal gradient; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, &vofNp1, &dvofdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    } 
  }

  // Homogeneous kernel implementation
  stk::topology partTopo = part->topology();
  auto& solverAlgMap = solverAlgDriver_->solverAlgorithmMap_;
  
  AssembleElemSolverAlgorithm* solverAlg = nullptr;
  bool solverAlgWasBuilt = false;
  
  std::tie(solverAlg, solverAlgWasBuilt) = build_or_add_part_to_solver_alg(*this, *part, solverAlgMap);
  
  ElemDataRequests& dataPreReqs = solverAlg->dataNeededByKernels_;
  auto& activeKernels = solverAlg->activeKernels_;
  
  if (solverAlgWasBuilt) {
    build_topo_kernel_if_requested<VolumeOfFluidMassElemKernel>
      (partTopo, *this, activeKernels, "lumped_mass",
       realm_.bulk_data(), *realm_.solutionOptions_, vof_, dataPreReqs, true);

    build_topo_kernel_if_requested<VolumeOfFluidMassElemKernel>
      (partTopo, *this, activeKernels, "mass",
       realm_.bulk_data(), *realm_.solutionOptions_, vof_, dataPreReqs, false);

    build_topo_kernel_if_requested<VolumeOfFluidScvAdvElemKernel>
      (partTopo, *this, activeKernels, "scv_advection",
       realm_.bulk_data(), *realm_.solutionOptions_, vof_, dataPreReqs);

    build_topo_kernel_if_requested<VolumeOfFluidScsAdvElemKernel>
      (partTopo, *this, activeKernels, "scs_advection",
       realm_.bulk_data(), *realm_.solutionOptions_, vof_, dataPreReqs);

    build_topo_kernel_if_requested<VolumeOfFluidDivElemKernel>
      (partTopo, *this, activeKernels, "div_correction",
       realm_.bulk_data(), *realm_.solutionOptions_, vof_, dataPreReqs);

    build_topo_kernel_if_requested<VolumeOfFluidSucvNsoElemKernel>
      (partTopo, *this, activeKernels, "sucv_nso",
       realm_.bulk_data(), *realm_.solutionOptions_, vof_, 1.0, 1.0, 1.0, dataPreReqs);

    build_topo_kernel_if_requested<VolumeOfFluidSucvNsoElemKernel>
      (partTopo, *this, activeKernels, "sucv",
       realm_.bulk_data(), *realm_.solutionOptions_, vof_, 1.0, 0.0, 1.0, dataPreReqs);

    build_topo_kernel_if_requested<VolumeOfFluidSucvNsoElemKernel>
      (partTopo, *this, activeKernels, "nso",
       realm_.bulk_data(), *realm_.solutionOptions_, vof_, 0.0, 1.0, 0.0, dataPreReqs);
  
    build_topo_kernel_if_requested<VolumeOfFluidSharpenElemKernel>
      (partTopo, *this, activeKernels, "sharpen",
       realm_.bulk_data(), *realm_.solutionOptions_, vof_, cAlpha_, dataPreReqs);

    build_topo_kernel_if_requested<VolumeOfFluidGclElemKernel>
      (partTopo, *this, activeKernels, "lumped_gcl",
       realm_.bulk_data(), *realm_.solutionOptions_, vof_, dataPreReqs, true);

    build_topo_kernel_if_requested<VolumeOfFluidGclElemKernel>
      (partTopo, *this, activeKernels, "gcl",
       realm_.bulk_data(), *realm_.solutionOptions_, vof_, dataPreReqs, false);
      
    report_invalid_supp_alg_names();
    report_built_supp_alg_names();
  }
}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{

  // algorithm type
  const AlgorithmType algType = INFLOW;

  ScalarFieldType &vofNp1 = (smooth_) ? vofSmoothed_->field_of_state(stk::mesh::StateNone) 
    : vof_->field_of_state(stk::mesh::StateNP1);  
  VectorFieldType &dvofdxNone = dvofdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; vof_bc
  ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "volume_of_fluid_bc"));
  stk::mesh::put_field_on_mesh(*theBcField, *part, nullptr);

  // extract the value for user specified vof and save off the AuxFunction
  InflowUserData userData = inflowBCData.userData_;
  std::string vofName = "volume_of_fluid";
  UserDataType theDataType = get_bc_data_type(userData, vofName);

  AuxFunction *theAuxFunc = NULL;
  if ( CONSTANT_UD == theDataType ) {
    VolumeOfFluid vof = userData.vof_;
    std::vector<double> userSpec(1);
    userSpec[0] = vof.vof_;

    // new it
    theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);
  }
  else {
    throw std::runtime_error("VolumeOfFluidEquationSystem::register_inflow_bc: only constant supported");
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

  // copy vof_bc to volume_of_fluid np1...
  CopyFieldAlgorithm *theCopyAlg
    = new CopyFieldAlgorithm(realm_, part,
                             theBcField, &vofNp1,
                             0, 1,
                             stk::topology::NODE_RANK);
  bcDataMapAlg_.push_back(theCopyAlg);

  // non-solver; dvofdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &vofNp1, &dvofdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // Dirichlet bc
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
    DirichletBC *theAlg
      = new DirichletBC(realm_, this, part, &vofNp1, theBcField, 0, 1);
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
VolumeOfFluidEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &partTopo,
  const OpenBoundaryConditionData &openBCData)
{

  // algorithm type
  const AlgorithmType algType = OPEN;

  ScalarFieldType &vofNp1 = (smooth_) ? vofSmoothed_->field_of_state(stk::mesh::StateNone) 
    : vof_->field_of_state(stk::mesh::StateNP1);  
  VectorFieldType &dvofdxNone = dvofdx_->field_of_state(stk::mesh::StateNone);

  // check if scs advection is active. If so, we need entrainment and open bc algorithms
  if ( scsAdvection_ ) {
    
    stk::mesh::MetaData &meta_data = realm_.meta_data();

    // first, vdot at open bc; register field
    MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(partTopo);
    const int numScsBip = meFC->numIntPoints_;    
    GenericFieldType *volBip 
      = &(meta_data.declare_field<GenericFieldType>(static_cast<stk::topology::rank_t>(meta_data.side_rank()), 
                                                   "open_volume_flow_rate"));
    stk::mesh::put_field_on_mesh(*volBip, *part, numScsBip, nullptr);

    // register boundary data; volume_of_fluid at bc
    ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "open_vof_bc"));
    stk::mesh::put_field_on_mesh(*theBcField, *part, nullptr);
    
    // extract the value for user specified mixFrac and save off the AuxFunction
    OpenUserData userData = openBCData.userData_;
    VolumeOfFluid vof = userData.vof_;
    std::vector<double> userSpec(1);
    userSpec[0] = vof.vof_;
    
    // new it
    ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);
    
    // bc data alg
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 theBcField, theAuxFunc,
                                 stk::topology::NODE_RANK);
    bcDataAlg_.push_back(auxAlg);

    auto& solverAlgMap = solverAlgDriver_->solverAlgorithmMap_;
    
    stk::topology elemTopo = get_elem_topo(realm_, *part);
    
    AssembleFaceElemSolverAlgorithm* faceElemSolverAlg = nullptr;
    bool solverAlgWasBuilt = false;
    
    std::tie(faceElemSolverAlg, solverAlgWasBuilt) 
      = build_or_add_part_to_face_elem_solver_alg(algType, *this, *part, elemTopo, solverAlgMap, "open");
    
    auto& activeKernels = faceElemSolverAlg->activeKernels_;
    
    if (solverAlgWasBuilt) {
  
      build_face_elem_topo_kernel_automatic<VolumeOfFluidOpenAdvElemKernel>
        (partTopo, elemTopo, *this, activeKernels, "volume_of_fluid_open",
         realm_.meta_data(), *realm_.solutionOptions_,
         this, vof_, theBcField, 
         faceElemSolverAlg->faceDataNeeded_, faceElemSolverAlg->elemDataNeeded_);
      
    }

  }

  // non-solver; dvofdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &vofNp1, &dvofdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const WallBoundaryConditionData &wallBCData)
{

  // algorithm type
  const AlgorithmType algType = WALL;

  ScalarFieldType &vofNp1 = (smooth_) ? vofSmoothed_->field_of_state(stk::mesh::StateNone) 
    : vof_->field_of_state(stk::mesh::StateNP1);  
  VectorFieldType &dvofdxNone = dvofdx_->field_of_state(stk::mesh::StateNone);
 
  stk::mesh::MetaData & meta_data = realm_.meta_data();
 
  // extract the value for [optional] user specified vof
  WallUserData userData = wallBCData.userData_;
  std::string vofName = "volume_of_fluid";
  if ( bc_data_specified(userData, vofName) ) {

    // FIXME: Generalize for constant vs function
    
    ScalarFieldType &realVofNp1 = vof_->field_of_state(stk::mesh::StateNP1);  
    
    // register boundary data; vof_bc
    ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "vof_bc"));
    stk::mesh::put_field_on_mesh(*theBcField, *part, nullptr);

    // extract data
    std::vector<double> userSpec(1);
    VolumeOfFluid vof = userData.vof_;
    userSpec[0] = vof.vof_;

    // new it
    ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

    // bc data alg
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 theBcField, theAuxFunc,
                                 stk::topology::NODE_RANK);
    bcDataAlg_.push_back(auxAlg);

    // copy vof_bc to vof np1...
    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               theBcField, &realVofNp1,
                               0, 1,
                               stk::topology::NODE_RANK);
    bcDataMapAlg_.push_back(theCopyAlg);

    // Dirichlet bc
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
      solverAlgDriver_->solverDirichAlgMap_.find(algType);
    if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
      DirichletBC *theAlg
        = new DirichletBC(realm_, this, part, &realVofNp1, theBcField, 0, 1);
      solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
    }
    else {
      itd->second->partVec_.push_back(part);
    }
  }

  // non-solver; dvofdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &vofNp1, &dvofdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }
}
  
//--------------------------------------------------------------------------
//-------- register_symmetry_bc --------------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const SymmetryBoundaryConditionData &/*wallBCData*/)
{

  // algorithm type
  const AlgorithmType algType = SYMMETRY;

  ScalarFieldType &vofNp1 = (smooth_) ? vofSmoothed_->field_of_state(stk::mesh::StateNone) 
    : vof_->field_of_state(stk::mesh::StateNP1);  
  VectorFieldType &dvofdxNone = dvofdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; dvofdx; allow for element-based shifted
  if ( !managePNG_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &vofNp1, &dvofdxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }
}

//--------------------------------------------------------------------------
//-------- register_overset_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::register_overset_bc()
{
  create_constraint_algorithm(vof_);
  stk::mesh::MetaData &metaData = realm_.meta_data();

  const int nDim = metaData.spatial_dimension();

  UpdateOversetFringeAlgorithmDriver* theAlg = new UpdateOversetFringeAlgorithmDriver(realm_);
  // Perform fringe updates before all equation system solves
  equationSystems_.preIterAlgDriver_.push_back(theAlg);
  theAlg->fields_.push_back(
    std::unique_ptr<OversetFieldData>(new OversetFieldData(vof_,1,1)));
  theAlg->fields_.push_back(
    std::unique_ptr<OversetFieldData>(new OversetFieldData(interfaceNormal_,1,nDim)));
  theAlg->fields_.push_back(
    std::unique_ptr<OversetFieldData>(new OversetFieldData(interfaceCurvature_,1,1)));
    
  if ( realm_.has_mesh_motion() ) {
    UpdateOversetFringeAlgorithmDriver* theAlgPost = new UpdateOversetFringeAlgorithmDriver(realm_,false);
    // Perform fringe updates after all equation system solves (ideally on the post_time_step)
    equationSystems_.postIterAlgDriver_.push_back(theAlgPost);
    theAlgPost->fields_.push_back(std::unique_ptr<OversetFieldData>(new OversetFieldData(vof_,1,1)));

    if (realm_.number_of_states()>2)
    {
      ScalarFieldType &vofN = vof_->field_of_state(stk::mesh::StateN);
      theAlgPost->fields_.push_back(std::unique_ptr<OversetFieldData>(new OversetFieldData(&vofN,1,1)));
    }
  }
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::reinitialize_linear_system()
{

  // delete linsys
  delete linsys_;

  // delete old solver
  const EquationType theEqID = EQ_VOLUME_OF_FLUID;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("volume_of_fluid");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_VOLUME_OF_FLUID);
  linsys_ = LinearSystem::create(realm_, 1, this, solver);

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- register_initial_condition_fcn ----------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::register_initial_condition_fcn(
  stk::mesh::Part *part,
  const std::map<std::string, std::string> &theNames,
  const std::map<std::string, std::vector<double> > &theParams)
{
  // iterate map and check for name
  const std::string dofName = "volume_of_fluid";
  std::map<std::string, std::string>::const_iterator iterName
    = theNames.find(dofName);
  if (iterName != theNames.end()) {
    std::string fcnName = (*iterName).second;

    AuxFunction *theAuxFunc = NULL;
    std::vector<double> fcnParams;
    
    // extract the params
    std::map<std::string, std::vector<double> >::const_iterator iterParams
      = theParams.find(dofName);
    if (iterParams != theParams.end()) {
      fcnParams = (*iterParams).second;	
    }

    if ( fcnName == "RayleighTaylor" ) {
      // create the function
      theAuxFunc = new RayleighTaylorMixFracAuxFunction(fcnParams);      
    }
    else if ( fcnName == "FixedHeight" ) {
      // create the function
      theAuxFunc = new FixedHeightMixFracAuxFunction(fcnParams);      
    }
    else {
      throw std::runtime_error("VolumeOfFluidEquationSystem::register_initial_condition_fcn: limited functions supported");
    }
    
    // create the algorithm
    AuxFunctionAlgorithm *auxAlg
      = new AuxFunctionAlgorithm(realm_, part,
				 vof_, theAuxFunc,
				 stk::topology::NODE_RANK);
    
    // push to ic
    realm_.initCondAlg_.push_back(auxAlg);
    
  }
}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::solve_and_update()
{

  // compute dvof/dx
  if ( isInit_ ) {
    sharpen_interface_explicit();
    smooth_vof();
    compute_interface_normal();
    compute_interface_curvature();
    compute_projected_nodal_gradient();
    isInit_ = false;
  }

  for ( int k = 0; k < maxIterations_; ++k ) {

    NaluEnv::self().naluOutputP0() << " " << k+1 << "/" << maxIterations_
                    << std::setw(15) << std::right << userSuppliedName_ << std::endl;
    
    // vof assemble, load_complete and solve
    assemble_and_solve(vofTmp_);

    // update
    double timeA = NaluEnv::self().nalu_time();
    update_and_clip();
    double timeB = NaluEnv::self().nalu_time();
    timerAssemble_ += (timeB-timeA);
    
    // sharpen/smoothing
    sharpen_interface_explicit();
    smooth_vof();
    compute_interface_normal();
    compute_interface_curvature();
    
    // projected nodal gradient
    compute_projected_nodal_gradient();
  }

  // allow for property evaluation
  if ( realm_.solutionOptions_->balancedForce_ ) {
    NaluEnv::self().naluOutputP0() 
      << "VolumeOfFluidEquationSystem::solve_and_update(): evaluate_properties" << std::endl;
    realm_.evaluate_properties();
  }
}

//--------------------------------------------------------------------------
//-------- update_and_clip -------------------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::update_and_clip()
{
  const double deltaVof = deltaVofClip_;
  const double lowBound = 0.0-deltaVof;
  const double highBound = 1.0+deltaVof;
  size_t numClip[2] = {0,0};
  double minVof = +1.0e16;
  double maxVof = -1.0e16;

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // define some common selectors
  stk::mesh::Selector s_all_nodes
    = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
    &stk::mesh::selectField(*vof_);

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin();
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    double *vof = stk::mesh::field_data(*vof_, b);
    double *vofTmp    = stk::mesh::field_data(*vofTmp_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      double vofNp1 = vof[k] + vofTmp[k];
      // clip now
      if ( vofNp1 < lowBound ) {
        minVof = std::min(vofNp1, minVof);
        vofNp1 = lowBound;
        numClip[0]++;
      }
      else if ( vofNp1 > highBound ) {
        maxVof = std::max(vofNp1, maxVof);
        vofNp1 = highBound;
        numClip[1]++;
      }
      vof[k] = vofNp1;
    }
  }

  // parallel assemble clipped value
  if ( outputClippingDiag_ ) {
    size_t g_numClip[2] = {};
    stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
    stk::all_reduce_sum(comm, numClip, g_numClip, 2);

    if ( g_numClip[0] > 0 ) {
      double g_minVof = 0;
      stk::all_reduce_min(comm, &minVof, &g_minVof, 1);
      NaluEnv::self().naluOutputP0() << "vof clipped (-) " << g_numClip[0] << " times; min: " << g_minVof << std::endl;
    }
    else {
      NaluEnv::self().naluOutputP0() << "vof clipped (-) zero times" << std::endl;
    }

    if ( g_numClip[1] > 0 ) {
      double g_maxVof = 0;
      stk::all_reduce_max(comm, &maxVof, &g_maxVof, 1);
      NaluEnv::self().naluOutputP0() << "vof clipped (+) " << g_numClip[1] << " times; max: " << g_maxVof << std::endl;
    }
    else {
      NaluEnv::self().naluOutputP0() << "vof clipped (+) zero times" << std::endl;
    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_interface_normal ----------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::compute_interface_normal()
{
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  
  const int nDim = metaData.spatial_dimension();
  const double small = 1.0e-16;

  // extract nodal fields
  VectorFieldType *coordinates 
    = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  ScalarFieldType *dualNodalVolume 
    = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  // interface normal is generally computed based on vofSmoothed_; user defines this
  ScalarFieldType &vofNp1 = (smooth_) ? vofSmoothed_->field_of_state(stk::mesh::StateNone) 
    : vof_->field_of_state(stk::mesh::StateNP1);  
  
  // zero assembled normal
  field_fill( metaData, bulkData, 0.0, *interfaceNormal_, realm_.get_activate_aura());

  // integration point data that is fixed
  std::vector<double> dvofdxIp(nDim);
  double *p_dvofdxIp = &dvofdxIp[0];

  // fields to gather
  std::vector<double> ws_coordinates;
  std::vector<double> ws_vofNp1;
  std::vector<double> ws_dualVolume;
  std::vector<double> ws_scVolume;
  std::vector<double> ws_dndx;
  std::vector<double> ws_deriv;
  std::vector<double> ws_det_j;

  // select locally owned where vof is defined; exclude inactive block
  stk::mesh::Selector s_locally_owned_union = metaData.locally_owned_part()
    & stk::mesh::selectField(vofNp1)
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(b.topology());

    // extract master element specifics
    const int nodesPerElement = meSCV->nodesPerElement_;
    const int numScvIp = meSCV->numIntPoints_;
    const int *ipNodeMap = meSCV->ipNodeMap();

    // algorithm related
    ws_coordinates.resize(nodesPerElement*nDim);
    ws_vofNp1.resize(nodesPerElement);
    ws_dualVolume.resize(nodesPerElement);
    ws_scVolume.resize(numScvIp);
    ws_dndx.resize(nDim*numScvIp*nodesPerElement);
    ws_deriv.resize(nDim*numScvIp*nodesPerElement);
    ws_det_j.resize(numScvIp);

    // pointers
    double *p_coordinates = &ws_coordinates[0];
    double *p_vofNp1 = &ws_vofNp1[0];
    double *p_dualVolume = &ws_dualVolume[0];
    double *p_scVolume = &ws_scVolume[0];
    double *p_dndx = &ws_dndx[0];
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const * node_rels = b.begin_nodes(k);
      int num_nodes = b.num_nodes(k);

      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerElement );

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // pointers to real data
        const double * coords = stk::mesh::field_data(*coordinates, node);

        // gather scalars
        p_vofNp1[ni] = *stk::mesh::field_data(vofNp1, node);
        p_dualVolume[ni] = *stk::mesh::field_data(*dualNodalVolume, node);

        // gather vectors
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          p_coordinates[offSet+j] = coords[j];
        }
      }
      
      // compute geometry
      double scvError = 0.0;
      meSCV->determinant(1, &p_coordinates[0], &p_scVolume[0], &scvError);

      // compute dndx
      meSCV->grad_op(1, &p_coordinates[0], &p_dndx[0], &ws_deriv[0], &ws_det_j[0], &scvError);
      
      for ( int ip = 0; ip < numScvIp; ++ip ) {
        
        // zero local ip gradient
        for ( int j = 0; j < nDim; ++j ) {
          p_dvofdxIp[j] = 0.0;
        }
        
        // compute gradient
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const int offSetDnDx = nDim*nodesPerElement*ip + ic*nDim;
          const double nodalVof = p_vofNp1[ic];
          for ( int j = 0; j < nDim; ++j ) {
            p_dvofdxIp[j] += p_dndx[offSetDnDx+j]*nodalVof;
          }
        }
        
        // magnitude
        double dvofMag = small;
        for ( int j = 0; j < nDim; ++j ) {
          dvofMag += p_dvofdxIp[j]*p_dvofdxIp[j];
        }
        dvofMag = std::sqrt(dvofMag);

        // nearest node for this ip
        const int nn = ipNodeMap[ip];
        stk::mesh::Entity node = node_rels[nn];

        // pointers to real data
        double * interfaceNormal = stk::mesh::field_data(*interfaceNormal_, node);

        const double volumeFac = p_scVolume[ip]/p_dualVolume[nn];
        for ( int j = 0; j < nDim; ++j )
          interfaceNormal[j] += p_dvofdxIp[j]/dvofMag*volumeFac;
      }
    }
  }
  
  // parallel sum
  stk::mesh::parallel_sum(bulkData, {interfaceNormal_});
  
  // periodic assemble
  if ( realm_.hasPeriodic_) {
    realm_.periodic_field_update(interfaceNormal_, nDim);
  }

  // overset update
  if ( realm_.hasOverset_ ) {
    realm_.overset_constraint_node_field_update(interfaceNormal_, 1, nDim);
  }
}

//--------------------------------------------------------------------------
//-------- compute_interface_curvature -------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::compute_interface_curvature()
{
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  
  const int nDim = metaData.spatial_dimension();

  // extract nodal fields
  VectorFieldType *coordinates 
    = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  ScalarFieldType *dualNodalVolume 
    = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");
  
  // zero assembled interface curvature
  field_fill( metaData, bulkData, 0.0, *interfaceCurvature_, realm_.get_activate_aura());

  // fields to gather
  std::vector<double> ws_coordinates;
  std::vector<double> ws_interfaceNormal;
  std::vector<double> ws_dualVolume;
  std::vector<double> ws_scVolume;
  std::vector<double> ws_dndx;
  std::vector<double> ws_deriv;
  std::vector<double> ws_det_j;

  // select locally owned where vof is defined; exclude inactive block
  stk::mesh::Selector s_locally_owned_union = metaData.locally_owned_part()
    & stk::mesh::selectField(*interfaceCurvature_)
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = elem_buckets.begin();
        ib != elem_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    MasterElement *meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(b.topology());

    // extract master element specifics
    const int nodesPerElement = meSCV->nodesPerElement_;
    const int numScvIp = meSCV->numIntPoints_;
    const int *ipNodeMap = meSCV->ipNodeMap();

    // algorithm related
    ws_coordinates.resize(nodesPerElement*nDim);
    ws_interfaceNormal.resize(nodesPerElement*nDim);
    ws_dualVolume.resize(nodesPerElement);
    ws_scVolume.resize(numScvIp);
    ws_dndx.resize(nDim*numScvIp*nodesPerElement);
    ws_deriv.resize(nDim*numScvIp*nodesPerElement);
    ws_det_j.resize(numScvIp);

    // pointers
    double *p_coordinates = &ws_coordinates[0];
    double *p_interfaceNormal = &ws_interfaceNormal[0];
    double *p_dualVolume = &ws_dualVolume[0];
    double *p_scVolume = &ws_scVolume[0];
    double *p_dndx = &ws_dndx[0];
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const * node_rels = b.begin_nodes(k);
      int num_nodes = b.num_nodes(k);

      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerElement );

      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];

        // pointers to real data
        const double * coords = stk::mesh::field_data(*coordinates, node);
        const double * interfaceNormal = stk::mesh::field_data(*interfaceNormal_, node);

        // gather scalars
        p_dualVolume[ni] = *stk::mesh::field_data(*dualNodalVolume, node);

        // gather vectors
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          p_coordinates[offSet+j] = coords[j];
          p_interfaceNormal[offSet+j] = interfaceNormal[j];
        }
      }
      
      // compute geometry
      double scvError = 0.0;
      meSCV->determinant(1, &p_coordinates[0], &p_scVolume[0], &scvError);

      // compute dndx
      meSCV->grad_op(1, &p_coordinates[0], &p_dndx[0], &ws_deriv[0], &ws_det_j[0], &scvError);
      
      for ( int ip = 0; ip < numScvIp; ++ip ) {

        double divN = 0.0;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          const int offSetDnDx = nDim*nodesPerElement*ip + ic*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            divN += p_dndx[offSetDnDx+j]*p_interfaceNormal[ic*nDim+j];
          }
        }
        
        // nearest node for this ip
        const int nn = ipNodeMap[ip];
        stk::mesh::Entity node = node_rels[nn];

        // pointers to real data
        double *interfaceCurvature = stk::mesh::field_data(*interfaceCurvature_, node);

        // augment nodal curvature
        *interfaceCurvature += -divN*p_scVolume[ip]/p_dualVolume[nn];
      }
    }
  }
  
  // parallel sum
  stk::mesh::parallel_sum(bulkData, {interfaceCurvature_});
  
  // periodic assemble
  if ( realm_.hasPeriodic_) {
    realm_.periodic_field_update(interfaceCurvature_, 1);
  }

  // overset update
  if ( realm_.hasOverset_ ) {
    realm_.overset_constraint_node_field_update(interfaceCurvature_, 1, 1);
  }
}


//--------------------------------------------------------------------------
//-------- smooth_vof ------------------------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::smooth_vof()
{

  // may not want to smooth
  if ( !smooth_)
    return;

  // otherwise, proceed
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();

  // Copy vof to vofSmoothed
  field_copy(metaData, bulkData, *vof_, *vofSmoothed_, realm_.get_activate_aura());
  
  // fixed set of iterations...
  for( int k = 0; k < smoothIter_; ++k) {
    smooth_vof_execute();
    // vofSmoothed = Fo_*dxMin_*dxMin_*smoothedRhs_ + 1.0*vofSmoothed_
    field_axpby(metaData, bulkData, Fo_*dxMin_*dxMin_, *smoothedRhs_, 1.0, *vofSmoothed_, realm_.get_activate_aura());
  }

}

//--------------------------------------------------------------------------
//-------- smooth_vof_execute ----------------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::smooth_vof_execute()
{
  // compute  smoothedRhs =  (Fo*dx*dx) d^2(vof)/dxj^2 using a low-order lumped projection
  // leave off smoothing factors until axpby and fold in the dxMin_ here..
  double l_dxMin = 1.0e16;
  
  stk::mesh::MetaData & metaData = realm_.meta_data();
  stk::mesh::BulkData & bulkData = realm_.bulk_data();
  
  const int nDim = metaData.spatial_dimension();
  
  // extract nodal fields
  VectorFieldType *coordinates 
    = metaData.get_field<VectorFieldType>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  ScalarFieldType *dualNodalVolume 
    = metaData.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

  // zero it
  field_fill(metaData, bulkData, 0.0, *smoothedRhs_, realm_.get_activate_aura());
  
  // geometry related to populate
  std::vector<double> ws_vof;
  std::vector<double> ws_dualVolume;
  std::vector<double> ws_coordinates;

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = metaData.locally_owned_part()
    & stk::mesh::selectField(*vofSmoothed_) 
    & !(realm_.get_inactive_selector());
  
  stk::mesh::BucketVector const& elem_buckets =
    realm_.get_buckets( stk::topology::ELEMENT_RANK, s_locally_owned_union );
  
  std::vector<double> ws_scsAreav;
  std::vector<double> ws_dndx;
  std::vector<double> ws_deriv;
  std::vector<double> ws_detJ;
  
  for ( const stk::mesh::Bucket* bucket_ptr : elem_buckets ) {
    const stk::mesh::Bucket & b = *bucket_ptr ;
    const stk::mesh::Bucket::size_type length   = b.size();
    
    // extract master element
    MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());
    
    // extract master element specifics
    const int nodesPerElement = meSCS->nodesPerElement_;
    const int numScsIp = meSCS->numIntPoints_;
    const int *lrscv = meSCS->adjacentNodes();
    
    // algorithm related
    ws_coordinates.resize(nodesPerElement*nDim);
    ws_vof.resize(nodesPerElement);
    ws_dualVolume.resize(nodesPerElement);
    ws_scsAreav.resize(numScsIp*nDim);
    ws_dndx.resize(nDim*numScsIp*nodesPerElement);
    ws_deriv.resize(nDim*numScsIp*nodesPerElement);
    ws_detJ.resize(numScsIp);
    
    // pointers
    double *p_coordinates = &ws_coordinates[0];
    double *p_vof = &ws_vof[0];
    double *p_dualVolume = &ws_dualVolume[0];
    double *p_scsAreav = &ws_scsAreav[0];
    double *p_dndx = &ws_dndx[0];
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      stk::mesh::Entity const *  node_rels = b.begin_nodes(k);
      int num_nodes = b.num_nodes(k);
      
      // sanity check on num nodes
      ThrowAssert( num_nodes == nodesPerElement );
      
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = node_rels[ni];
        
        // pointers to real data
        const double * coords = stk::mesh::field_data(*coordinates, node );
        
        // gather scalars
        p_vof[ni] = *stk::mesh::field_data(*vofSmoothed_, node );
        p_dualVolume[ni] = *stk::mesh::field_data(*dualNodalVolume, node );
        
        // gather vectors
        const int niNdim = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          p_coordinates[niNdim+j] = coords[j];
        }
      }
      
      // compute geometry
      double scsError = 0.0;
      meSCS->determinant(1, &p_coordinates[0], &p_scsAreav[0], &scsError);
      meSCS->grad_op(1, &p_coordinates[0], &p_dndx[0], &ws_deriv[0], &ws_detJ[0], &scsError);
      
      for ( int ip = 0; ip < numScsIp; ++ip ) {
        
        // left and right nodes for this ip
        const int il = lrscv[2*ip];
        const int ir = lrscv[2*ip+1];
        
        // determine dx; edge distance magnitude
        double dx = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double dxj = p_coordinates[ir*nDim+j] - p_coordinates[il*nDim+j];
          dx += dxj*dxj;
        }
        l_dxMin = std::min(l_dxMin, std::sqrt(dx));
        
        stk::mesh::Entity nodeL = node_rels[il];
        stk::mesh::Entity nodeR = node_rels[ir];
        
        // pointer to fields to assemble
        double *smoothedRhsL = stk::mesh::field_data(*smoothedRhs_, nodeL );
        double *smoothedRhsR = stk::mesh::field_data(*smoothedRhs_, nodeR );
        
        double qDiff = 0.0;
        for ( int ic = 0; ic < nodesPerElement; ++ic ) {
          
          double lhsfacDiff = 0.0;
          const int offSetDnDx = nDim*nodesPerElement*ip + ic*nDim;
          for ( int j = 0; j < nDim; ++j ) {
            lhsfacDiff += -p_dndx[offSetDnDx+j]*p_scsAreav[ip*nDim+j];
          }
          qDiff += lhsfacDiff*p_vof[ic];
        }
        
        *smoothedRhsL -= qDiff/ws_dualVolume[il];
        *smoothedRhsR += qDiff/ws_dualVolume[ir];
      }
    }
  }
  
  // parallel sum
  stk::mesh::parallel_sum(bulkData, {smoothedRhs_});

  // periodic update
  if ( realm_.hasPeriodic_) {
    realm_.periodic_field_update(smoothedRhs_, 1);
  }

  // overset update
  if ( realm_.hasOverset_ ) {
    realm_.overset_constraint_node_field_update(smoothedRhs_, 1, 1);
  }
  
  // compute min
  stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
  stk::all_reduce_min(comm, &l_dxMin, &dxMin_, 1);
}

//--------------------------------------------------------------------------
//-------- sharpen_interface_explicit --------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::sharpen_interface_explicit()
{
  // not yet supported
}

//--------------------------------------------------------------------------
//-------- predict_state ---------------------------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::predict_state()
{
  // copy state n to state np1
  ScalarFieldType &vofN = vof_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &vofNp1 = vof_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), vofN, vofNp1, realm_.get_activate_aura());
}

//--------------------------------------------------------------------------
//-------- manage_projected_nodal_gradient ---------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::manage_projected_nodal_gradient(
  EquationSystems& eqSystems)
{
  throw std::runtime_error("VofEquationSystem::manage_projected_nodal_gradient: Not yet supported");
}

//--------------------------------------------------------------------------
//-------- compute_projected_nodal_gradient---------------------------------
//--------------------------------------------------------------------------
void
VolumeOfFluidEquationSystem::compute_projected_nodal_gradient()
{
  // projected nodal gradient is generally using vofSmoothed_
  if ( !managePNG_ ) {
    const double timeA = -NaluEnv::self().nalu_time();
    assembleNodalGradAlgDriver_->execute();
    timerMisc_ += (NaluEnv::self().nalu_time() + timeA);
  }
  else {
    projectedNodalGradEqs_->solve_and_update_external();
  }
}

} // namespace nalu
} // namespace Sierra
