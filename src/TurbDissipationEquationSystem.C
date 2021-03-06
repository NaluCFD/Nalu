/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "TurbDissipationEquationSystem.h"
#include "AlgorithmDriver.h"
#include "AssembleScalarEdgeOpenSolverAlgorithm.h"
#include "AssembleScalarEdgeSolverAlgorithm.h"
#include "AssembleScalarElemSolverAlgorithm.h"
#include "AssembleScalarElemOpenSolverAlgorithm.h"
#include "AssembleScalarNonConformalSolverAlgorithm.h"
#include "AssembleNodeSolverAlgorithm.h"
#include "AssembleNodalGradAlgorithmDriver.h"
#include "AssembleNodalGradEdgeAlgorithm.h"
#include "AssembleNodalGradElemAlgorithm.h"
#include "AssembleNodalGradBoundaryAlgorithm.h"
#include "AssembleNodalGradNonConformalAlgorithm.h"
#include "AuxFunctionAlgorithm.h"
#include "ComputeWallModelTurbDissipationWallAlgorithm.h"
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
#include "Realm.h"
#include "Realms.h"
#include "ScalarGclNodeSuppAlg.h"
#include "ScalarMassBackwardEulerNodeSuppAlg.h"
#include "ScalarMassBDF2NodeSuppAlg.h"
#include "ScalarMassElemSuppAlgDep.h"
#include "Simulation.h"
#include "SolutionOptions.h"
#include "TimeIntegrator.h"
#include "TurbDissipationKEpsilonNodeSourceSuppAlg.h"
#include "SolverAlgorithmDriver.h"

// template for supp algs
#include "AlgTraits.h"
#include "kernel/KernelBuilder.h"
#include "kernel/KernelBuilderLog.h"

// consolidated
#include "AssembleElemSolverAlgorithm.h"
#include "kernel/ScalarMassElemKernel.h"
#include "kernel/ScalarAdvDiffElemKernel.h"
#include "kernel/ScalarUpwAdvDiffElemKernel.h"
#include "kernel/TurbDissipationKEpsilonSrcElemKernel.h"

// nso
#include "nso/ScalarNSOElemKernel.h"
#include "nso/ScalarNSOKeElemSuppAlg.h"
#include "nso/ScalarNSOElemSuppAlgDep.h"

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

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// TurbDissipationEquationSystem - manages epsilon pde system
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
TurbDissipationEquationSystem::TurbDissipationEquationSystem(
  EquationSystems& eqSystems)
  : EquationSystem(eqSystems, "TurbDissEQS","turbulent_dissipation"),
    managePNG_(realm_.get_consistent_mass_matrix_png("turbulent_dissipation")),
    eps_(NULL),
    dedx_(NULL),
    eTmp_(NULL),
    visc_(NULL),
    tvisc_(NULL),
    evisc_(NULL),
    epsWallBc_(NULL),
    assembledWallEps_(NULL),
    assembledWallArea_(NULL),
    assembleNodalGradAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "turbulent_dissipation", "dedx")),
    diffFluxCoeffAlgDriver_(new AlgorithmDriver(realm_))
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("turbulent_dissipation");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_TURBULENT_DISS);
  linsys_ = LinearSystem::create(realm_, 1, this, solver);

  // determine nodal gradient form
  set_nodal_gradient("turbulent_dissipation");
  NaluEnv::self().naluOutputP0() << "Edge projected nodal gradient for turbulent dissipation: " << edgeNodalGradient_ <<std::endl;

  // push back EQ to manager
  realm_.push_equation_to_systems(this);

  // create projected nodal gradient equation system
  if ( managePNG_ )
    throw std::runtime_error("TurbDissipationEquationSystem::Error managePNG is not supported");
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
TurbDissipationEquationSystem::~TurbDissipationEquationSystem()
{
  delete assembleNodalGradAlgDriver_;
  delete diffFluxCoeffAlgDriver_;
  std::vector<Algorithm *>::iterator ii;
  for( ii=wallModelAlg_.begin(); ii!=wallModelAlg_.end(); ++ii )
    delete *ii;
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
TurbDissipationEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  // register dof; set it as a restart variable
  eps_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_dissipation", numStates));
  stk::mesh::put_field_on_mesh(*eps_, *part, nullptr);
  realm_.augment_restart_variable_list("turbulent_dissipation");

  dedx_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dedx"));
  stk::mesh::put_field_on_mesh(*dedx_, *part, nDim, nullptr);

  // delta solution for linear solver; share delta since this is a split system
  eTmp_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "eTmp"));
  stk::mesh::put_field_on_mesh(*eTmp_, *part, nullptr);

  visc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity"));
  stk::mesh::put_field_on_mesh(*visc_, *part, nullptr);

  tvisc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity"));
  stk::mesh::put_field_on_mesh(*tvisc_, *part, nullptr);

  evisc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "effective_viscosity_eps"));
  stk::mesh::put_field_on_mesh(*evisc_, *part, nullptr);

  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && (!realm_.restarted_simulation() || realm_.support_inconsistent_restart()) ) {
    ScalarFieldType &epsN = eps_->field_of_state(stk::mesh::StateN);
    ScalarFieldType &epsNp1 = eps_->field_of_state(stk::mesh::StateNP1);

    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               &epsNp1, &epsN,
                               0, 1,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlg);
  }
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
TurbDissipationEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{

  // types of algorithms
  const AlgorithmType algType = INTERIOR;

  ScalarFieldType &epsNp1 = eps_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dedxNone = dedx_->field_of_state(stk::mesh::StateNone);

  // non-solver, dedx; allow for element-based shifted
  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleNodalGradAlgDriver_->algMap_.find(algType);
  if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg = NULL;
    if ( edgeNodalGradient_ && realm_.realmUsesEdges_ ) {
      theAlg = new AssembleNodalGradEdgeAlgorithm(realm_, part, &epsNp1, &dedxNone);
    }
    else {
      theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, &epsNp1, &dedxNone, edgeNodalGradient_);
    }
    assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

  // solver; interior contribution (advection + diffusion)
  if (!realm_.solutionOptions_->useConsolidatedSolverAlg_) {

    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
      = solverAlgDriver_->solverAlgMap_.find(algType);
    if (itsi == solverAlgDriver_->solverAlgMap_.end()) {
      SolverAlgorithm* theAlg = NULL;
      if (realm_.realmUsesEdges_) {
        theAlg = new AssembleScalarEdgeSolverAlgorithm(realm_, part, this, eps_, dedx_, evisc_);
      }
      else {
        theAlg = new AssembleScalarElemSolverAlgorithm(realm_, part, this, eps_, dedx_, evisc_);
      }
      solverAlgDriver_->solverAlgMap_[algType] = theAlg;

      // look for fully integrated source terms
      std::map<std::string, std::vector<std::string> >::iterator isrc
        = realm_.solutionOptions_->elemSrcTermsMap_.find("turbulent_dissipation");
      if (isrc != realm_.solutionOptions_->elemSrcTermsMap_.end()) {

        if (realm_.realmUsesEdges_)
          throw std::runtime_error("TurbDissipationElemSrcTerms::Error can not use element source terms for an edge-based scheme");

        std::vector<std::string> mapNameVec = isrc->second;
        for (size_t k = 0; k < mapNameVec.size(); ++k) {
          std::string sourceName = mapNameVec[k];
          SupplementalAlgorithm* suppAlg = NULL;
          if (sourceName == "NSO_2ND_ALT") {
            suppAlg = new ScalarNSOElemSuppAlgDep(realm_, eps_, dedx_, evisc_, 0.0, 1.0);
          }
          else if (sourceName == "NSO_4TH_ALT") {
            suppAlg = new ScalarNSOElemSuppAlgDep(realm_, eps_, dedx_, evisc_, 1.0, 1.0);
          }
          else if (sourceName == "NSO_2ND_KE") {
            const double turbSc = realm_.get_turb_schmidt(eps_->name());
            suppAlg = new ScalarNSOKeElemSuppAlg(realm_, eps_, dedx_, turbSc, 0.0);
          }
          else if (sourceName == "NSO_4TH_KE") {
            const double turbSc = realm_.get_turb_schmidt(eps_->name());
            suppAlg = new ScalarNSOKeElemSuppAlg(realm_, eps_, dedx_, turbSc, 1.0);
          }
          else if (sourceName == "turbulent_dissipation_time_derivative" ) {
            suppAlg = new ScalarMassElemSuppAlgDep(realm_, eps_, false);
          }
          else if (sourceName == "lumped_turbulent_dissipation_time_derivative" ) {
            suppAlg = new ScalarMassElemSuppAlgDep(realm_, eps_, true);
          }
          else {
            throw std::runtime_error("TurbDissipationElemSrcTerms::Error Source term is not supported: " + sourceName);
          }
          NaluEnv::self().naluOutputP0() << "TurbDissipationElemSrcTerms::added() " << sourceName << std::endl;
          theAlg->supplementalAlg_.push_back(suppAlg);
        }
      }
    }
    else {
      itsi->second->partVec_.push_back(part);
    }

    // time term; src; both nodally lumped
    const AlgorithmType algMass = MASS;
    // Check if the user has requested CMM or LMM algorithms; if so, do not
    // include Nodal Mass algorithms
    std::vector<std::string> checkAlgNames = {
      "turbulent_dissipation_time_derivative",
      "lumped_turbulent_dissipation_time_derivative"};
    bool elementMassAlg = supp_alg_is_requested(checkAlgNames);
    std::map<AlgorithmType, SolverAlgorithm*>::iterator itsm =
      solverAlgDriver_->solverAlgMap_.find(algMass);
    if (itsm == solverAlgDriver_->solverAlgMap_.end()) {
      // create the solver alg
      AssembleNodeSolverAlgorithm *theAlg
        = new AssembleNodeSolverAlgorithm(realm_, part, this);
      solverAlgDriver_->solverAlgMap_[algMass] = theAlg;

      // now create the supplemental alg for mass term
      if (!elementMassAlg) {
        if (realm_.number_of_states() == 2) {
          ScalarMassBackwardEulerNodeSuppAlg *theMass
            = new ScalarMassBackwardEulerNodeSuppAlg(realm_, eps_);
          theAlg->supplementalAlg_.push_back(theMass);
        }
        else {
          ScalarMassBDF2NodeSuppAlg *theMass
            = new ScalarMassBDF2NodeSuppAlg(realm_, eps_);
          theAlg->supplementalAlg_.push_back(theMass);
        }
      }

      // now create the src alg for eps source
      SupplementalAlgorithm *theSrc = NULL;
      switch(realm_.solutionOptions_->turbulenceModel_) {
      case KEPS:
        {
          theSrc = new TurbDissipationKEpsilonNodeSourceSuppAlg(realm_);
        }
        break;
      default:
        throw std::runtime_error("Unsupported turbulence model in TurbDissipation rate: only k_espison supported");
      }
      theAlg->supplementalAlg_.push_back(theSrc);

      // Add nodal src term supp alg...; limited number supported
      std::map<std::string, std::vector<std::string> >::iterator isrc
        = realm_.solutionOptions_->srcTermsMap_.find("turbulent_dissipation");
      if (isrc != realm_.solutionOptions_->srcTermsMap_.end()) {
        std::vector<std::string> mapNameVec = isrc->second;
        for (size_t k = 0; k < mapNameVec.size(); ++k) {
          std::string sourceName = mapNameVec[k];
          SupplementalAlgorithm* suppAlg = NULL;
          if (sourceName == "gcl") {
            suppAlg = new ScalarGclNodeSuppAlg(eps_, realm_);
          }
          else {
            throw std::runtime_error("TurbDissipationNodalSrcTerms::Error Source term is not supported: " + sourceName);
          }
          NaluEnv::self().naluOutputP0() << "TurbDissipationNodalSrcTerms::added() " << sourceName << std::endl;
          theAlg->supplementalAlg_.push_back(suppAlg);
        }
      }
    }
    else {
      itsm->second->partVec_.push_back(part);
    }
  }
  else {
    // Homogeneous kernel implementation
    if (realm_.realmUsesEdges_)
      throw std::runtime_error("TurbDissipationEquationSystem::Error can not use element source terms for an edge-based scheme");

    stk::topology partTopo = part->topology();
    auto& solverAlgMap = solverAlgDriver_->solverAlgorithmMap_;

    AssembleElemSolverAlgorithm* solverAlg = nullptr;
    bool solverAlgWasBuilt = false;

    std::tie(solverAlg, solverAlgWasBuilt) =
      build_or_add_part_to_solver_alg(*this, *part, solverAlgMap);

    ElemDataRequests& dataPreReqs = solverAlg->dataNeededByKernels_;
    auto& activeKernels = solverAlg->activeKernels_;

    if (solverAlgWasBuilt) {
      build_topo_kernel_if_requested<ScalarMassElemKernel>
        (partTopo, *this, activeKernels, "turbulent_dissipation_time_derivative",
         realm_.bulk_data(), *realm_.solutionOptions_, eps_, dataPreReqs, false);

      build_topo_kernel_if_requested<ScalarMassElemKernel>
        (partTopo, *this, activeKernels,  "lumped_turbulent_dissipation_time_derivative",
         realm_.bulk_data(), *realm_.solutionOptions_, eps_, dataPreReqs, true);

      build_topo_kernel_if_requested<ScalarAdvDiffElemKernel>
        (partTopo, *this, activeKernels, "advection_diffusion",
         realm_.bulk_data(), *realm_.solutionOptions_, eps_, evisc_, dataPreReqs);

      build_topo_kernel_if_requested<ScalarUpwAdvDiffElemKernel>
        (partTopo, *this, activeKernels, "upw_advection_diffusion",
        realm_.bulk_data(), *realm_.solutionOptions_, this, eps_, dedx_, evisc_, dataPreReqs);

      build_topo_kernel_if_requested<TurbDissipationKEpsilonSrcElemKernel>
        (partTopo, *this, activeKernels, "k_epsilon",
         realm_.bulk_data(), *realm_.solutionOptions_, dataPreReqs, true);

      build_topo_kernel_if_requested<ScalarNSOElemKernel>
        (partTopo, *this, activeKernels, "NSO_2ND",
         realm_.bulk_data(), *realm_.solutionOptions_, eps_, dedx_, evisc_, 0.0, 0.0, dataPreReqs);

      build_topo_kernel_if_requested<ScalarNSOElemKernel>
        (partTopo, *this, activeKernels, "NSO_2ND_ALT",
         realm_.bulk_data(), *realm_.solutionOptions_, eps_, dedx_, evisc_, 0.0, 1.0, dataPreReqs);

      build_topo_kernel_if_requested<ScalarNSOElemKernel>
        (partTopo, *this, activeKernels, "NSO_4TH",
         realm_.bulk_data(), *realm_.solutionOptions_, eps_, dedx_, evisc_, 1.0, 0.0, dataPreReqs);

      build_topo_kernel_if_requested<ScalarNSOElemKernel>
        (partTopo, *this, activeKernels, "NSO_4TH_ALT",
         realm_.bulk_data(), *realm_.solutionOptions_, eps_, dedx_, evisc_, 1.0, 1.0, dataPreReqs);

      report_invalid_supp_alg_names();
      report_built_supp_alg_names();
    }
  }

  // effective diffusive flux coefficient alg for tubrulent dissipation
  std::map<AlgorithmType, Algorithm *>::iterator itev =
    diffFluxCoeffAlgDriver_->algMap_.find(algType);
  if ( itev == diffFluxCoeffAlgDriver_->algMap_.end() ) {
    const double lamSc = realm_.get_lam_schmidt(eps_->name());
    const double turbSc = realm_.get_turb_schmidt(eps_->name());
    EffectiveDiffFluxCoeffAlgorithm *effDiffAlg 
      = new EffectiveDiffFluxCoeffAlgorithm(realm_, part, visc_, tvisc_, evisc_, lamSc, turbSc);
    diffFluxCoeffAlgDriver_->algMap_[algType] = effDiffAlg;
  }
  else {
    itev->second->partVec_.push_back(part);
  }

}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
TurbDissipationEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{

  // algorithm type
  const AlgorithmType algType = INFLOW;

  ScalarFieldType &epsNp1 = eps_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dedxNone = dedx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; eps_bc
  ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_dissipation_bc"));
  stk::mesh::put_field_on_mesh(*theBcField, *part, nullptr);

  // extract the value for user specified tke and save off the AuxFunction
  InflowUserData userData = inflowBCData.userData_;
  TurbDiss eps = userData.eps_;
  std::vector<double> userSpec(1);
  userSpec[0] = eps.turbDiss_;

  // new it
  ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

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

  // copy eps_bc to dissipation np1...
  CopyFieldAlgorithm *theCopyAlg
    = new CopyFieldAlgorithm(realm_, part,
                             theBcField, &epsNp1,
                             0, 1,
                             stk::topology::NODE_RANK);
  bcDataMapAlg_.push_back(theCopyAlg);

  // non-solver; dedx; allow for element-based shifted
  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleNodalGradAlgDriver_->algMap_.find(algType);
  if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg 
      = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &epsNp1, &dedxNone, edgeNodalGradient_);
    assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

  // Dirichlet bc
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
    DirichletBC *theAlg
      = new DirichletBC(realm_, this, part, &epsNp1, theBcField, 0, 1);
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
TurbDissipationEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const OpenBoundaryConditionData &openBCData)
{

  // algorithm type
  const AlgorithmType algType = OPEN;

  ScalarFieldType &epsNp1 = eps_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dedxNone = dedx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; eps_bc
  ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "open_eps_bc"));
  stk::mesh::put_field_on_mesh(*theBcField, *part, nullptr);

  // extract the value for user specified tke and save off the AuxFunction
  OpenUserData userData = openBCData.userData_;
  TurbDiss eps = userData.eps_;
  std::vector<double> userSpec(1);
  userSpec[0] = eps.turbDiss_;

  // new it
  ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               theBcField, theAuxFunc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlg);

  // non-solver; dedx; allow for element-based shifted
  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleNodalGradAlgDriver_->algMap_.find(algType);
  if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg 
      = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &epsNp1, &dedxNone, edgeNodalGradient_);
    assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

  // solver open; lhs
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    SolverAlgorithm *theAlg = NULL;
    if ( realm_.realmUsesEdges_ ) {
      theAlg = new AssembleScalarEdgeOpenSolverAlgorithm(realm_, part, this, eps_, theBcField, &dedxNone, evisc_);
    }
    else {
      theAlg = new AssembleScalarElemOpenSolverAlgorithm(realm_, part, this, eps_, theBcField, &dedxNone, evisc_);
    }
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;
  }
  else {
    itsi->second->partVec_.push_back(part);
  }

}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
TurbDissipationEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const WallBoundaryConditionData &wallBCData)
{

  // algorithm type
  const AlgorithmType algType = WALL;

  // np1
  ScalarFieldType &epsNp1 = eps_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dedxNone = dedx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; eps_bc
  epsWallBc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "eps_bc"));
  stk::mesh::put_field_on_mesh(*epsWallBc_, *part, nullptr);

  // need to register the assembles wall value for eps; can not share with eps_bc
  assembledWallEps_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "wall_model_eps_bc"));
  stk::mesh::put_field_on_mesh(*assembledWallEps_, *part, nullptr);

  assembledWallArea_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_wall_area_eps"));
  stk::mesh::put_field_on_mesh(*assembledWallArea_, *part, nullptr);

  // wall function support only 
  WallUserData userData = wallBCData.userData_;
  bool anyWallFunctionActivated = userData.wallFunctionApproach_ || userData.wallFunctionProjectedApproach_;
  if ( !anyWallFunctionActivated )
    throw std::runtime_error("KEpsilon only supports a wall function approach");
  
  // create proper algorithms to fill nodal eps and assembled wall area; utau managed by momentum
  ComputeWallModelTurbDissipationWallAlgorithm *wallAlg 
    = new ComputeWallModelTurbDissipationWallAlgorithm(realm_, part);
  wallModelAlg_.push_back(wallAlg);
  
  // Dirichlet bc
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
    solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
    DirichletBC *theAlg =
      new DirichletBC(realm_, this, part, &epsNp1, epsWallBc_, 0, 1);
    solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
  }
  else {
    itd->second->partVec_.push_back(part);
  }

  // non-solver; dedx; allow for element-based shifted
  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleNodalGradAlgDriver_->algMap_.find(algType);
  if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg 
      = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &epsNp1, &dedxNone, edgeNodalGradient_);
    assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_symmetry_bc --------------------------------------------
//--------------------------------------------------------------------------
void
TurbDissipationEquationSystem::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const SymmetryBoundaryConditionData &symmetryBCData)
{

  // algorithm type
  const AlgorithmType algType = SYMMETRY;

  // np1
  ScalarFieldType &epsNp1 = eps_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dedxNone = dedx_->field_of_state(stk::mesh::StateNone);

  // non-solver; dedx; allow for element-based shifted
  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleNodalGradAlgDriver_->algMap_.find(algType);
  if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg 
      = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &epsNp1, &dedxNone, edgeNodalGradient_);
    assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

}

//--------------------------------------------------------------------------
//-------- register_non_conformal_bc ---------------------------------------
//--------------------------------------------------------------------------
void
TurbDissipationEquationSystem::register_non_conformal_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/)
{

  const AlgorithmType algType = NON_CONFORMAL;

  // np1
  ScalarFieldType &epsNp1 = eps_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dedxNone = dedx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to dedx; DG algorithm decides on locations for integration points
  if ( edgeNodalGradient_ ) {    
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &epsNp1, &dedxNone, edgeNodalGradient_);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }
  else {
    // proceed with DG
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      AssembleNodalGradNonConformalAlgorithm *theAlg 
        = new AssembleNodalGradNonConformalAlgorithm(realm_, part, &epsNp1, &dedxNone);
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }
  }

  // solver; lhs; same for edge and element-based scheme
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
    solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    AssembleScalarNonConformalSolverAlgorithm *theAlg
      = new AssembleScalarNonConformalSolverAlgorithm(realm_, part, this, eps_, evisc_);
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;
  }
  else {
    itsi->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_overset_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
TurbDissipationEquationSystem::register_overset_bc()
{
  create_constraint_algorithm(eps_);

  UpdateOversetFringeAlgorithmDriver* theAlg = new UpdateOversetFringeAlgorithmDriver(realm_);
  // Perform fringe updates before all equation system solves
  equationSystems_.preIterAlgDriver_.push_back(theAlg);

  theAlg->fields_.push_back(
    std::unique_ptr<OversetFieldData>(new OversetFieldData(eps_,1,1)));

  if ( realm_.has_mesh_motion() ) {
    UpdateOversetFringeAlgorithmDriver* theAlgPost = new UpdateOversetFringeAlgorithmDriver(realm_,false);
    // Perform fringe updates after all equation system solves (ideally on the post_time_step)
    equationSystems_.postIterAlgDriver_.push_back(theAlgPost);
    theAlgPost->fields_.push_back(std::unique_ptr<OversetFieldData>(new OversetFieldData(eps_,1,1)));
  }
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
TurbDissipationEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
TurbDissipationEquationSystem::reinitialize_linear_system()
{

  // delete linsys
  delete linsys_;

  // delete old solver
  const EquationType theEqID = EQ_TURBULENT_DISS;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("turbulent_dissipation");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_TURBULENT_DISS);
  linsys_ = LinearSystem::create(realm_, 1, this, solver);

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- compute_projected_nodal_gradient() ------------------------------
//--------------------------------------------------------------------------
void
TurbDissipationEquationSystem::compute_projected_nodal_gradient()
{
  const double timeA = -NaluEnv::self().nalu_time();
  assembleNodalGradAlgDriver_->execute();
  timerMisc_ += (NaluEnv::self().nalu_time() + timeA);
}

//--------------------------------------------------------------------------
//-------- compute_effective_flux_coeff() ----------------------------------
//--------------------------------------------------------------------------
void
TurbDissipationEquationSystem::compute_effective_diff_flux_coeff()
{
  const double timeA = -NaluEnv::self().nalu_time();
  diffFluxCoeffAlgDriver_->execute();
  timerMisc_ += (NaluEnv::self().nalu_time() + timeA);
}

//--------------------------------------------------------------------------
//-------- compute_wall_model_parameters() ---------------------------------
//--------------------------------------------------------------------------
void
TurbDissipationEquationSystem::compute_wall_model_parameters()
{

  // check if we need to process anything
  if ( wallModelAlg_.size() == 0 )
    return;

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // selector; all nodes that have a epsilon-specific nodal field registered
  stk::mesh::Selector s_all_nodes
     = (meta_data.locally_owned_part() | meta_data.globally_shared_part())
     &stk::mesh::selectField(*assembledWallArea_);
  stk::mesh::BucketVector const& node_buckets
    = realm_.get_buckets( stk::topology::NODE_RANK, s_all_nodes );

  // zero the fields
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
      ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length  = b.size();
    double * assembledWallEps = stk::mesh::field_data(*assembledWallEps_, b);
    double * assembledWallArea = stk::mesh::field_data(*assembledWallArea_, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      assembledWallEps[k] = 0.0;
      assembledWallArea[k] = 0.0;
    }
  }

  // process the algorithm(s)
  for ( size_t k = 0; k < wallModelAlg_.size(); ++k ) {
    wallModelAlg_[k]->execute();
  }

  // parallel assemble
  stk::mesh::parallel_sum(bulk_data, {assembledWallEps_, assembledWallArea_});

  // periodic assemble
  if ( realm_.hasPeriodic_) {
    const unsigned fieldSize = 1;
    const bool bypassFieldCheck = false;
    realm_.periodic_field_update(assembledWallEps_, fieldSize, bypassFieldCheck);
    realm_.periodic_field_update(assembledWallArea_, fieldSize, bypassFieldCheck);
  }

  // normalize and set assembled eps to eps bc
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
      ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length  = b.size();
    double * eps = stk::mesh::field_data(*eps_, b);
    double * epsWallBc = stk::mesh::field_data(*epsWallBc_, b);
    double * assembledWallEps = stk::mesh::field_data(*assembledWallEps_, b);
    double * assembledWallArea = stk::mesh::field_data(*assembledWallArea_, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const double epsBnd = assembledWallEps[k]/assembledWallArea[k];
      epsWallBc[k] = epsBnd;
      assembledWallEps[k] = epsBnd;
      // make sure that the next matrix assembly uses the proper eps value
      eps[k] = epsBnd;
    }
  }
}

//--------------------------------------------------------------------------
//-------- predict_state() -------------------------------------------------
//--------------------------------------------------------------------------
void
TurbDissipationEquationSystem::predict_state()
{
  // copy state n to state np1
  ScalarFieldType &epsN = eps_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &epsNp1 = eps_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), epsN, epsNp1, realm_.get_activate_aura());
}

} // namespace nalu
} // namespace Sierra
