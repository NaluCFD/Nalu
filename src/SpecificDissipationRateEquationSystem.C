/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <SpecificDissipationRateEquationSystem.h>
#include <AlgorithmDriver.h>
#include <AssembleScalarEdgeContactSolverAlgorithm.h>
#include <AssembleScalarEdgeOpenSolverAlgorithm.h>
#include <AssembleScalarEdgeSolverAlgorithm.h>
#include <AssembleScalarElemSolverAlgorithm.h>
#include <AssembleScalarElemOpenSolverAlgorithm.h>
#include <AssembleScalarNonConformalSolverAlgorithm.h>
#include <AssembleNodeSolverAlgorithm.h>
#include <AssembleNodalGradAlgorithmDriver.h>
#include <AssembleNodalGradEdgeAlgorithm.h>
#include <AssembleNodalGradElemAlgorithm.h>
#include <AssembleNodalGradBoundaryAlgorithm.h>
#include <AssembleNodalGradEdgeContactAlgorithm.h>
#include <AssembleNodalGradElemContactAlgorithm.h>
#include <AssembleNodalGradNonConformalAlgorithm.h>
#include <AuxFunctionAlgorithm.h>
#include <ComputeLowReynoldsSDRWallAlgorithm.h>
#include <ComputeWallModelSDRWallAlgorithm.h>
#include <ConstantAuxFunction.h>
#include <CopyFieldAlgorithm.h>
#include <DirichletBC.h>
#include <EffectiveSSTDiffFluxCoeffAlgorithm.h>
#include <EquationSystem.h>
#include <EquationSystems.h>
#include <Enums.h>
#include <FieldFunctions.h>
#include <LinearSolvers.h>
#include <LinearSolver.h>
#include <LinearSystem.h>
#include <NaluEnv.h>
#include <NaluParsing.h>
#include <Realm.h>
#include <Realms.h>
#include <ScalarGclNodeSuppAlg.h>
#include <ScalarMassBackwardEulerNodeSuppAlg.h>
#include <ScalarMassBDF2NodeSuppAlg.h>
#include <ScalarMassElemSuppAlg.h>
#include <ScalarKeNSOElemSuppAlg.h>
#include <ScalarNSOElemSuppAlg.h>
#include <Simulation.h>
#include <SolutionOptions.h>
#include <TimeIntegrator.h>
#include <SpecificDissipationRateSSTNodeSourceSuppAlg.h>
#include <SolverAlgorithmDriver.h>

// stk_util
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/CPUTime.hpp>

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
// SpecificDissipationRateEquationSystem - manages sdr pde system
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
SpecificDissipationRateEquationSystem::SpecificDissipationRateEquationSystem(
  EquationSystems& eqSystems)
  : EquationSystem(eqSystems, "SpecDissRateEQS"),
    sdr_(NULL),
    dwdx_(NULL),
    wTmp_(NULL),
    visc_(NULL),
    tvisc_(NULL),
    evisc_(NULL),
    sdrWallBc_(NULL),
    assembledWallSdr_(NULL),
    assembledWallArea_(NULL),
    assembleNodalGradAlgDriver_(new AssembleNodalGradAlgorithmDriver(realm_, "specific_dissipation_rate", "dwdx")),
    diffFluxCoeffAlgDriver_(new AlgorithmDriver(realm_))
{
  // extract solver name and solver object
  std::string solverName = realm_.equationSystems_.get_solver_block_name("specific_dissipation_rate");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_SPEC_DISS_RATE);
  linsys_ = LinearSystem::create(realm_, 1, name_, solver);

  // determine nodal gradient form
  set_nodal_gradient("specific_dissipation_rate");
  NaluEnv::self().naluOutputP0() << "Edge projected nodal gradient for specific_dissipation_rate: " << edgeNodalGradient_ <<std::endl;

  // push back EQ to manager
  realm_.push_equation_to_systems(this);
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
SpecificDissipationRateEquationSystem::~SpecificDissipationRateEquationSystem()
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
SpecificDissipationRateEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  // register dof; set it as a restart variable
  sdr_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "specific_dissipation_rate", numStates));
  stk::mesh::put_field(*sdr_, *part);
  realm_.augment_restart_variable_list("specific_dissipation_rate");

  dwdx_ =  &(meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "dwdx"));
  stk::mesh::put_field(*dwdx_, *part, nDim);

  // delta solution for linear solver; share delta since this is a split system
  wTmp_ =  &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "wTmp"));
  stk::mesh::put_field(*wTmp_, *part);

  visc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity"));
  stk::mesh::put_field(*visc_, *part);

  tvisc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "turbulent_viscosity"));
  stk::mesh::put_field(*tvisc_, *part);

  evisc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "effective_viscosity_sdr"));
  stk::mesh::put_field(*evisc_, *part);

  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 && (!realm_.restarted_simulation() || realm_.support_inconsistent_restart()) ) {
    ScalarFieldType &sdrN = sdr_->field_of_state(stk::mesh::StateN);
    ScalarFieldType &sdrNp1 = sdr_->field_of_state(stk::mesh::StateNP1);

    CopyFieldAlgorithm *theCopyAlg
      = new CopyFieldAlgorithm(realm_, part,
                               &sdrNp1, &sdrN,
                               0, 1,
                               stk::topology::NODE_RANK);
    copyStateAlg_.push_back(theCopyAlg);
  }
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
SpecificDissipationRateEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{

  // types of algorithms
  const AlgorithmType algType = INTERIOR;

  ScalarFieldType &sdrNp1 = sdr_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dwdxNone = dwdx_->field_of_state(stk::mesh::StateNone);

  // non-solver, dwdx; allow for element-based shifted
  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleNodalGradAlgDriver_->algMap_.find(algType);
  if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg = NULL;
    if ( edgeNodalGradient_ && realm_.realmUsesEdges_ ) {
      theAlg = new AssembleNodalGradEdgeAlgorithm(realm_, part, &sdrNp1, &dwdxNone);
    }
    else {
      theAlg = new AssembleNodalGradElemAlgorithm(realm_, part, &sdrNp1, &dwdxNone, edgeNodalGradient_);
    }
    assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

  // solver; interior contribution (advection + diffusion)
  bool useCMM = false;
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi
    = solverAlgDriver_->solverAlgMap_.find(algType);
  if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
    SolverAlgorithm *theAlg = NULL;
    if ( realm_.realmUsesEdges_ ) {
      theAlg = new AssembleScalarEdgeSolverAlgorithm(realm_, part, this, sdr_, dwdx_, evisc_);
    }
    else {
      theAlg = new AssembleScalarElemSolverAlgorithm(realm_, part, this, sdr_, dwdx_, evisc_);
    }
    solverAlgDriver_->solverAlgMap_[algType] = theAlg;

    // look for src
    std::map<std::string, std::vector<std::string> >::iterator isrc 
      = realm_.solutionOptions_->elemSrcTermsMap_.find("specific_dissipation_rate");
    if ( isrc != realm_.solutionOptions_->elemSrcTermsMap_.end() ) {
      std::vector<std::string> mapNameVec = isrc->second;
      for (size_t k = 0; k < mapNameVec.size(); ++k ) {
        std::string sourceName = mapNameVec[k];
        SupplementalAlgorithm *suppAlg = NULL;
        if (sourceName == "NSO_2ND_ALT" ) {
          suppAlg = new ScalarNSOElemSuppAlg(realm_, sdr_, dwdx_, evisc_, 0.0, 1.0);
        }
        else if (sourceName == "NSO_4TH_ALT" ) {
          suppAlg = new ScalarNSOElemSuppAlg(realm_, sdr_, dwdx_, evisc_, 1.0, 1.0);
        }
        else if (sourceName == "NSO_KE_2ND" ) {
          const double turbSc = realm_.get_turb_schmidt(sdr_->name());
          suppAlg = new ScalarKeNSOElemSuppAlg(realm_, sdr_, dwdx_, turbSc, 0.0);
        }
        else if (sourceName == "NSO_KE_4TH" ) {
          const double turbSc = realm_.get_turb_schmidt(sdr_->name());
          suppAlg = new ScalarKeNSOElemSuppAlg(realm_, sdr_, dwdx_, turbSc, 1.0);
        }
        else if (sourceName == "specific_dissipation_rate_time_derivative" ) {
          useCMM = true;
          suppAlg = new ScalarMassElemSuppAlg(realm_, sdr_);
        }
        else {
          throw std::runtime_error("SpecificDissipationElemSrcTerms::Error Source term is not supported: " + sourceName);
        }     
      }
    }
  }
  else {
    itsi->second->partVec_.push_back(part);
  }

  // time term; src; both nodally lumped
  const AlgorithmType algMass = MASS;
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itsm =
    solverAlgDriver_->solverAlgMap_.find(algMass);
  if ( itsm == solverAlgDriver_->solverAlgMap_.end() ) {
    // create the solver alg
    AssembleNodeSolverAlgorithm *theAlg
      = new AssembleNodeSolverAlgorithm(realm_, part, this);
    solverAlgDriver_->solverAlgMap_[algMass] = theAlg;

    // now create the supplemental alg for mass term
    if ( !useCMM ) {
      if ( realm_.number_of_states() == 2 ) {
        ScalarMassBackwardEulerNodeSuppAlg *theMass
          = new ScalarMassBackwardEulerNodeSuppAlg(realm_, sdr_);
        theAlg->supplementalAlg_.push_back(theMass);
      }
      else {
        ScalarMassBDF2NodeSuppAlg *theMass
          = new ScalarMassBDF2NodeSuppAlg(realm_, sdr_);
        theAlg->supplementalAlg_.push_back(theMass);
      }
    }

    // now create the src alg for sdr source
    SpecificDissipationRateSSTNodeSourceSuppAlg *theSrc
      = new SpecificDissipationRateSSTNodeSourceSuppAlg(realm_);
    theAlg->supplementalAlg_.push_back(theSrc);

    // Add src term supp alg...; limited number supported
    std::map<std::string, std::vector<std::string> >::iterator isrc 
      = realm_.solutionOptions_->srcTermsMap_.find("specific_dissipation_rate");
    if ( isrc != realm_.solutionOptions_->srcTermsMap_.end() ) {
      std::vector<std::string> mapNameVec = isrc->second;   
      for (size_t k = 0; k < mapNameVec.size(); ++k ) {
        std::string sourceName = mapNameVec[k];
        SupplementalAlgorithm *suppAlg = NULL;
        if ( sourceName == "gcl" ) {
          suppAlg = new ScalarGclNodeSuppAlg(sdr_,realm_);
        }
        else {
          throw std::runtime_error("SpecificDissipationRateNodalSrcTerms::Error Source term is not supported: " + sourceName);
        }
        // add supplemental algorithm
        theAlg->supplementalAlg_.push_back(suppAlg);
      }
    }

  }
  else {
    itsm->second->partVec_.push_back(part);
  }

  // effective diffusive flux coefficient alg for SST
  std::map<AlgorithmType, Algorithm *>::iterator itev =
    diffFluxCoeffAlgDriver_->algMap_.find(algType);
  if ( itev == diffFluxCoeffAlgDriver_->algMap_.end() ) {
    const double sigmaWOne = realm_.get_turb_model_constant(TM_sigmaWOne);
    const double sigmaWTwo = realm_.get_turb_model_constant(TM_sigmaWTwo);
    EffectiveSSTDiffFluxCoeffAlgorithm *effDiffAlg
      = new EffectiveSSTDiffFluxCoeffAlgorithm(realm_, part, visc_, tvisc_, evisc_, sigmaWOne, sigmaWTwo);
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
SpecificDissipationRateEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{

  // algorithm type
  const AlgorithmType algType = INFLOW;

  ScalarFieldType &sdrNp1 = sdr_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dwdxNone = dwdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; sdr_bc
  ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "sdr_bc"));
  stk::mesh::put_field(*theBcField, *part);

  // extract the value for user specified tke and save off the AuxFunction
  InflowUserData userData = inflowBCData.userData_;
  SpecDissRate sdr = userData.sdr_;
  std::vector<double> userSpec(1);
  userSpec[0] = sdr.specDissRate_;

  // new it
  ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               theBcField, theAuxFunc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlg);

  // copy sdr_bc to specific_dissipation_rate np1...
  CopyFieldAlgorithm *theCopyAlg
    = new CopyFieldAlgorithm(realm_, part,
                             theBcField, &sdrNp1,
                             0, 1,
                             stk::topology::NODE_RANK);
  bcDataMapAlg_.push_back(theCopyAlg);

  // non-solver; dwdx; allow for element-based shifted
  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleNodalGradAlgDriver_->algMap_.find(algType);
  if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg 
      = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &sdrNp1, &dwdxNone, edgeNodalGradient_);
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
      = new DirichletBC(realm_, this, part, &sdrNp1, theBcField, 0, 1);
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
SpecificDissipationRateEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const OpenBoundaryConditionData &openBCData)
{

  // algorithm type
  const AlgorithmType algType = OPEN;

  ScalarFieldType &sdrNp1 = sdr_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dwdxNone = dwdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; sdr_bc
  ScalarFieldType *theBcField = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "open_sdr_bc"));
  stk::mesh::put_field(*theBcField, *part);

  // extract the value for user specified tke and save off the AuxFunction
  OpenUserData userData = openBCData.userData_;
  SpecDissRate sdr = userData.sdr_;
  std::vector<double> userSpec(1);
  userSpec[0] = sdr.specDissRate_;

  // new it
  ConstantAuxFunction *theAuxFunc = new ConstantAuxFunction(0, 1, userSpec);

  // bc data alg
  AuxFunctionAlgorithm *auxAlg
    = new AuxFunctionAlgorithm(realm_, part,
                               theBcField, theAuxFunc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlg);

  // non-solver; dwdx; allow for element-based shifted
  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleNodalGradAlgDriver_->algMap_.find(algType);
  if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg 
      = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &sdrNp1, &dwdxNone, edgeNodalGradient_);
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
      theAlg = new AssembleScalarEdgeOpenSolverAlgorithm(realm_, part, this, sdr_, theBcField, &dwdxNone, evisc_);
    }
    else {
      theAlg = new AssembleScalarElemOpenSolverAlgorithm(realm_, part, this, sdr_, theBcField, &dwdxNone, evisc_);
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
SpecificDissipationRateEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const WallBoundaryConditionData &wallBCData)
{

  // algorithm type
  const AlgorithmType algType = WALL;

  // np1
  ScalarFieldType &sdrNp1 = sdr_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dwdxNone = dwdx_->field_of_state(stk::mesh::StateNone);

  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register boundary data; sdr_bc
  sdrWallBc_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "sdr_bc"));
  stk::mesh::put_field(*sdrWallBc_, *part);

  // need to register the assembles wall value for sdr; can not share with sdr_bc
  assembledWallSdr_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "wall_model_sdr_bc"));
  stk::mesh::put_field(*assembledWallSdr_, *part);

  assembledWallArea_ = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "assembled_wall_area_sdr"));
  stk::mesh::put_field(*assembledWallArea_, *part);

  // are we using wall functions or is this a low Re model?
  WallUserData userData = wallBCData.userData_;
  bool wallFunctionApproach = userData.wallFunctionApproach_;

  // create proper algorithms to fill nodal omega and assembled wall area; utau managed by momentum
  Algorithm *wallAlg = NULL;
  if ( wallFunctionApproach ) {
    wallAlg = new ComputeWallModelSDRWallAlgorithm(realm_, part, realm_.realmUsesEdges_);
  }
  else {
    wallAlg = new ComputeLowReynoldsSDRWallAlgorithm(realm_, part, realm_.realmUsesEdges_);
  }
  wallModelAlg_.push_back(wallAlg);

  // Dirichlet bc
  std::map<AlgorithmType, SolverAlgorithm *>::iterator itd =
      solverAlgDriver_->solverDirichAlgMap_.find(algType);
  if ( itd == solverAlgDriver_->solverDirichAlgMap_.end() ) {
    DirichletBC *theAlg =
        new DirichletBC(realm_, this, part, &sdrNp1, sdrWallBc_, 0, 1);
    solverAlgDriver_->solverDirichAlgMap_[algType] = theAlg;
  }
  else {
    itd->second->partVec_.push_back(part);
  }

  // non-solver; dwdx; allow for element-based shifted
  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleNodalGradAlgDriver_->algMap_.find(algType);
  if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg 
      = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &sdrNp1, &dwdxNone, edgeNodalGradient_);
    assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

}

//--------------------------------------------------------------------------
//-------- register_contact_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
SpecificDissipationRateEquationSystem::register_contact_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const ContactBoundaryConditionData &contactBCData) {

  const AlgorithmType algType = CONTACT;

  ScalarFieldType &sdrNp1 = sdr_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dwdxNone = dwdx_->field_of_state(stk::mesh::StateNone);
  if ( realm_.realmUsesEdges_ ) {

    // register halo_sdr if using the element-based projected nodal gradient
    ScalarFieldType *haloSdr = NULL;
    if ( !edgeNodalGradient_ ) {
      stk::mesh::MetaData &meta_data = realm_.meta_data();
      haloSdr = &(meta_data.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "halo_sdr"));
      stk::mesh::put_field(*haloSdr, *part);
    }

    // non-solver; contribution to dwdx
    std::map<AlgorithmType, Algorithm *>::iterator it =
      assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg = NULL;
      if ( edgeNodalGradient_ ) {
        theAlg = new AssembleNodalGradEdgeContactAlgorithm(realm_, part, &sdrNp1, &dwdxNone);
      }
      else {
        theAlg = new AssembleNodalGradElemContactAlgorithm(realm_, part, &sdrNp1, &dwdxNone, haloSdr);
      }
      assembleNodalGradAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it->second->partVec_.push_back(part);
    }

    // solver; lhs
    std::map<AlgorithmType, SolverAlgorithm *>::iterator itsi =
      solverAlgDriver_->solverAlgMap_.find(algType);
    if ( itsi == solverAlgDriver_->solverAlgMap_.end() ) {
      AssembleScalarEdgeContactSolverAlgorithm *theAlg
        = new AssembleScalarEdgeContactSolverAlgorithm(realm_, part, this, sdr_, dwdx_, evisc_);
      solverAlgDriver_->solverAlgMap_[algType] = theAlg;
    }
    else {
      itsi->second->partVec_.push_back(part);
    }
  }
  else {
    throw std::runtime_error("Sorry, element-based contact not supported");
  }
}

//--------------------------------------------------------------------------
//-------- register_symmetry_bc --------------------------------------------
//--------------------------------------------------------------------------
void
SpecificDissipationRateEquationSystem::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const SymmetryBoundaryConditionData &symmetryBCData)
{

  // algorithm type
  const AlgorithmType algType = SYMMETRY;

  // np1
  ScalarFieldType &sdrNp1 = sdr_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dwdxNone = dwdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; dwdx; allow for element-based shifted
  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleNodalGradAlgDriver_->algMap_.find(algType);
  if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg 
      = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &sdrNp1, &dwdxNone, edgeNodalGradient_);
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
SpecificDissipationRateEquationSystem::register_non_conformal_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/)
{

  const AlgorithmType algType = NON_CONFORMAL;

  // np1
  ScalarFieldType &sdrNp1 = sdr_->field_of_state(stk::mesh::StateNP1);
  VectorFieldType &dwdxNone = dwdx_->field_of_state(stk::mesh::StateNone);

  // non-solver; contribution to dwdx; DG algorithm decides on locations for integration points
  if ( edgeNodalGradient_ ) {    
    std::map<AlgorithmType, Algorithm *>::iterator it
      = assembleNodalGradAlgDriver_->algMap_.find(algType);
    if ( it == assembleNodalGradAlgDriver_->algMap_.end() ) {
      Algorithm *theAlg 
        = new AssembleNodalGradBoundaryAlgorithm(realm_, part, &sdrNp1, &dwdxNone, edgeNodalGradient_);
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
        = new AssembleNodalGradNonConformalAlgorithm(realm_, part, &sdrNp1, &dwdxNone);
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
      = new AssembleScalarNonConformalSolverAlgorithm(realm_, part, this, sdr_, evisc_);
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
SpecificDissipationRateEquationSystem::register_overset_bc()
{
  create_constraint_algorithm(sdr_);
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
SpecificDissipationRateEquationSystem::initialize()
{
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
SpecificDissipationRateEquationSystem::reinitialize_linear_system()
{

  // delete linsys
  delete linsys_;

  // delete old solver
  const EquationType theEqID = EQ_SPEC_DISS_RATE;
  LinearSolver *theSolver = NULL;
  std::map<EquationType, LinearSolver *>::const_iterator iter
    = realm_.root()->linearSolvers_->solvers_.find(theEqID);
  if (iter != realm_.root()->linearSolvers_->solvers_.end()) {
    theSolver = (*iter).second;
    delete theSolver;
  }

  // create new solver
  std::string solverName = realm_.equationSystems_.get_solver_block_name("specific_dissipation_rate");
  LinearSolver *solver = realm_.root()->linearSolvers_->create_solver(solverName, EQ_SPEC_DISS_RATE);
  linsys_ = LinearSystem::create(realm_, 1, name_, solver);

  // initialize
  solverAlgDriver_->initialize_connectivity();
  linsys_->finalizeLinearSystem();
}

//--------------------------------------------------------------------------
//-------- assemble_nodal_gradient() ---------------------------------------
//--------------------------------------------------------------------------
void
SpecificDissipationRateEquationSystem::assemble_nodal_gradient()
{
  const double timeA = -stk::cpu_time();
  assembleNodalGradAlgDriver_->execute();
  timerMisc_ += (stk::cpu_time() + timeA);
}

//--------------------------------------------------------------------------
//-------- compute_effective_flux_coeff() ----------------------------------
//--------------------------------------------------------------------------
void
SpecificDissipationRateEquationSystem::compute_effective_diff_flux_coeff()
{
  const double timeA = -stk::cpu_time();
  diffFluxCoeffAlgDriver_->execute();
  timerMisc_ += (stk::cpu_time() + timeA);
}

//--------------------------------------------------------------------------
//-------- compute_wall_model_parameters() ---------------------------------
//--------------------------------------------------------------------------
void
SpecificDissipationRateEquationSystem::compute_wall_model_parameters()
{

  // check if we need to process anything
  if ( wallModelAlg_.size() == 0 )
    return;

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  // selector; all nodes that have a SST-specific nodal field registered
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
    double * assembledWallSdr = stk::mesh::field_data(*assembledWallSdr_, b);
    double * assembledWallArea = stk::mesh::field_data(*assembledWallArea_, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      assembledWallSdr[k] = 0.0;
      assembledWallArea[k] = 0.0;
    }
  }

  // process the algorithm(s)
  for ( size_t k = 0; k < wallModelAlg_.size(); ++k ) {
    wallModelAlg_[k]->execute();
  }

  // parallel assemble
  std::vector<stk::mesh::FieldBase*> fields;
  fields.push_back(assembledWallSdr_);
  fields.push_back(assembledWallArea_);
  stk::mesh::parallel_sum(bulk_data, fields);

  // periodic assemble
  if ( realm_.hasPeriodic_) {
    const unsigned fieldSize = 1;
    const bool bypassFieldCheck = false; // fields are not defined at all slave/master node pairs
    realm_.periodic_field_update(assembledWallSdr_, fieldSize, bypassFieldCheck);
    realm_.periodic_field_update(assembledWallArea_, fieldSize, bypassFieldCheck);
  }

  // normalize and set assembled sdr to sdr bc
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
      ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length  = b.size();
    double * sdr = stk::mesh::field_data(*sdr_, b);
    double * sdrWallBc = stk::mesh::field_data(*sdrWallBc_, b);
    double * assembledWallSdr = stk::mesh::field_data(*assembledWallSdr_, b);
    double * assembledWallArea = stk::mesh::field_data(*assembledWallArea_, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const double sdrBnd = assembledWallSdr[k]/assembledWallArea[k];
      sdrWallBc[k] = sdrBnd;
      assembledWallSdr[k] = sdrBnd;
      // make sure that the next matrix assembly uses the proper sdr value
      sdr[k] = sdrBnd;
    }
  }
}

//--------------------------------------------------------------------------
//-------- predict_state() -------------------------------------------------
//--------------------------------------------------------------------------
void
SpecificDissipationRateEquationSystem::predict_state()
{
  // copy state n to state np1
  ScalarFieldType &sdrN = sdr_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &sdrNp1 = sdr_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), sdrN, sdrNp1, realm_.get_activate_aura());
}

} // namespace nalu
} // namespace Sierra
