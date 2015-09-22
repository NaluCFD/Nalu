/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <Realm.h>
#include <Simulation.h>
#include <NaluEnv.h>

// percept
#if defined (NALU_USES_PERCEPT)
#include <adapt/AdaptedMeshVerifier.hpp>
#include <adapt/AdaptHelperFunctions.hpp>
#include <percept/PerceptMesh.hpp>
#include <Adapter.h>
#endif

#include <AuxFunction.h>
#include <AuxFunctionAlgorithm.h>
#include <ComputeGeometryAlgorithmDriver.h>
#include <ComputeGeometryBoundaryAlgorithm.h>
#include <ComputeGeometryExtrusionBoundaryAlgorithm.h>
#include <ComputeGeometryInteriorAlgorithm.h>
#include <ConstantAuxFunction.h>
#include <ContactInfo.h>
#include <ContactManager.h>
#include <Enums.h>
#include <EquationSystem.h>
#include <EquationSystems.h>
#include <ErrorIndicatorAlgorithmDriver.h>
#include <ExtrusionMeshDistanceBoundaryAlgorithm.h>
#include <FieldTypeDef.h>
#include <LinearSystem.h>
#include <master_element/MasterElement.h>
#include <MaterialPropertys.h>
#include <NaluParsing.h>
#include <NonConformalManager.h>
#include <NonConformalInfo.h>
#include <OutputInfo.h>
#include <AveragingInfo.h>
#include <PostProcessingInfo.h>
#include <PostProcessingData.h>
#include <SolutionNormPostProcessing.h>
#include <PeriodicManager.h>
#include <Realms.h>
#include <TurbulenceAveragingAlgorithm.h>
#include <SolutionOptions.h>
#include <TimeIntegrator.h>

// overset
#include <overset/OversetManager.h>

// props; algs, evaluators and data
#include <property_evaluator/GenericPropAlgorithm.h>
#include <property_evaluator/HDF5TablePropAlgorithm.h>
#include <property_evaluator/InverseDualVolumePropAlgorithm.h>
#include <property_evaluator/InversePropAlgorithm.h>
#include <property_evaluator/TemperaturePropAlgorithm.h>
#include <property_evaluator/LinearPropAlgorithm.h>
#include <property_evaluator/ConstantPropertyEvaluator.h>
#include <property_evaluator/EnthalpyPropertyEvaluator.h>
#include <property_evaluator/IdealGasPropertyEvaluator.h>
#include <property_evaluator/PropertyEvaluator.h>
#include <property_evaluator/ReferencePropertyData.h>
#include <property_evaluator/SpecificHeatPropertyEvaluator.h>
#include <property_evaluator/SutherlandsPropertyEvaluator.h>
#include <property_evaluator/WaterPropertyEvaluator.h>
#include <property_evaluator/MaterialPropertyData.h>

// tables
#include <tabular_props/HDF5FilePtr.h>

// transfer
#include <xfer/Transfer.h>

// stk_util
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/CPUTime.hpp>
#include <stk_util/environment/perf_util.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/SkinMesh.hpp>

// stk_io
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>
#include <Ioss_SubSystem.h>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

// yaml for parsing..
#include <yaml-cpp/yaml.h>
#include <NaluParsing.h>
#include <NaluParsingHelper.h>

// boost
#include <boost/lexical_cast.hpp>

// basic c++
#include <map>
#include <cmath>
#include <utility>
#include <stdint.h>

#define USE_NALU_PERFORMANCE_TESTING_CALLGRIND 0
#if USE_NALU_PERFORMANCE_TESTING_CALLGRIND
#include "/usr/netpub/valgrind-3.8.1/include/valgrind/callgrind.h"
#endif

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// Realm - do some stuff
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
Realm::Realm(Realms& realms)
  : realms_(realms),
    inputDBName_("input_unknown"),
    spatialDimension_(3u),  // for convenience; can always get it from meta data
    realmUsesEdges_(false),
    solveFrequency_(1),
    isTurbulent_(false),
    needsEnthalpy_(false),
    l2Scaling_(1.0),
    metaData_(NULL),
    bulkData_(NULL),
    ioBroker_(NULL),
    resultsFileIndex_(99),
    restartFileIndex_(99),
    computeGeometryAlgDriver_(0),
    extrusionMeshDistanceAlgDriver_(0),
    errorIndicatorAlgDriver_(0),
#if defined (NALU_USES_PERCEPT)
    adapter_(0),
#endif
    numInitialElements_(0),
    postConvergedAlgDriver_(0),
    timeIntegrator_(0),
    boundaryConditions_(*this),
    initialConditions_(*this),
    materialPropertys_(*this),
    equationSystems_(*this),
    maxCourant_(0.0),
    maxReynolds_(0.0),
    targetCourant_(1.0),
    timeStepChangeFactor_(1.25),
    currentNonlinearIteration_(1),
    solutionOptions_(new SolutionOptions()),
    outputInfo_(new OutputInfo()),
    averagingInfo_(new AveragingInfo()),
    postProcessingInfo_(new PostProcessingInfo()),
    solutionNormPostProcessing_(new SolutionNormPostProcessing(*this)),
    nodeCount_(0),
    estimateMemoryOnly_(false),
    availableMemoryPerCoreGB_(0),
    timerReadMesh_(0.0),
    timerOutputFields_(0.0),
    timerCreateEdges_(0.0),
    timerContact_(0.0),
    timerInitializeEqs_(0.0),
    timerPropertyEval_(0.0),
    timerAdapt_(0.0),
    timerTransferSearch_(0.0),
    contactManager_(NULL),
    nonConformalManager_(NULL),
    oversetManager_(NULL),
    hasContact_(false),
    hasNonConformal_(false),
    hasOverset_(false),
    hasTransfer_(false),
    periodicManager_(NULL),
    hasPeriodic_(false),
    globalParameters_(),
    exposedBoundaryPart_(0),
    edgesPart_(0),
    checkForMissingBcs_(false),
    isothermalFlow_(true),
    uniformFlow_(true),
    provideEntityCount_(false),
    HDF5ptr_(NULL),
    autoDecompType_("None"),
    activateAura_(false),
    activateMemoryDiagnostic_(false),
    supportInconsistentRestart_(false)
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Realm::~Realm()
{

  delete bulkData_;
  delete metaData_;
  delete ioBroker_;

  delete computeGeometryAlgDriver_;
  if ( NULL != extrusionMeshDistanceAlgDriver_ )
    delete extrusionMeshDistanceAlgDriver_;

  if ( NULL != errorIndicatorAlgDriver_)
    delete errorIndicatorAlgDriver_;

#if defined (NALU_USES_PERCEPT)
  if ( NULL != adapter_)
    delete adapter_;
#endif

  if ( NULL != postConvergedAlgDriver_ )
    delete postConvergedAlgDriver_;

  // prop algs
  std::vector<Algorithm *>::iterator ii;
  for( ii=initCondAlg_.begin(); ii!=initCondAlg_.end(); ++ii )
    delete *ii;
  for( ii=propertyAlg_.begin(); ii!=propertyAlg_.end(); ++ii ) {
    delete *ii;
  }

  // any bc data
  std::vector<AuxFunctionAlgorithm *>::iterator iaux;
  for( iaux=bcDataAlg_.begin(); iaux!=bcDataAlg_.end(); ++iaux )
    delete *iaux;

  // post converged algs
  std::vector<Algorithm *>::iterator ipc;
  for( ipc=postConvergedAlg_.begin(); ipc!=postConvergedAlg_.end(); ++ipc )
    delete *ipc;

  // delete master elements that were saved off; surface
  std::map<stk::topology, MasterElement *>::iterator it;
  for ( it = surfaceMeMap_.begin(); it!= surfaceMeMap_.end(); ++it ) {
    MasterElement *theElem = it->second;
    delete theElem;
  }
  // volume
  for ( it = volumeMeMap_.begin(); it!= volumeMeMap_.end(); ++it ) {
    MasterElement *theElem = it->second;
    delete theElem;
  }

  delete solutionOptions_;
  delete outputInfo_;
  delete averagingInfo_;
  delete postProcessingInfo_;
  delete solutionNormPostProcessing_;

  // delete contact related things
  if ( NULL != contactManager_ )
    delete contactManager_;

  // delete non-conformal related things
  if ( NULL != nonConformalManager_ )
    delete nonConformalManager_;

  // delete periodic related things
  if ( NULL != periodicManager_ )
    delete periodicManager_;

  // delete HDF5 file ptr
  if ( NULL != HDF5ptr_ )
    delete HDF5ptr_;
}

void
Realm::breadboard()
{
  computeGeometryAlgDriver_ = new ComputeGeometryAlgorithmDriver(*this);
  equationSystems_.breadboard();
}

bool
Realm::debug() const
{
  return root()->debug_;
}

//--------------------------------------------------------------------------
//-------- get_activate_memory_diagnostic ----------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_activate_memory_diagnostic()
{
  return activateMemoryDiagnostic_;
}

//--------------------------------------------------------------------------
//-------- provide_memory_summary ------------------------------------------
//--------------------------------------------------------------------------
void
Realm::provide_memory_summary() 
{
  size_t now, hwm;
  stk::get_memory_usage(now, hwm);
  // min, max, sum
  size_t global_now[3] = {now,now,now};
  size_t global_hwm[3] = {hwm,hwm,hwm};
  
  stk::all_reduce(NaluEnv::self().parallel_comm(), stk::ReduceSum<1>( &global_now[2] ) );
  stk::all_reduce(NaluEnv::self().parallel_comm(), stk::ReduceMin<1>( &global_now[0] ) );
  stk::all_reduce(NaluEnv::self().parallel_comm(), stk::ReduceMax<1>( &global_now[1] ) );
  
  stk::all_reduce(NaluEnv::self().parallel_comm(), stk::ReduceSum<1>( &global_hwm[2] ) );
  stk::all_reduce(NaluEnv::self().parallel_comm(), stk::ReduceMin<1>( &global_hwm[0] ) );
  stk::all_reduce(NaluEnv::self().parallel_comm(), stk::ReduceMax<1>( &global_hwm[1] ) );
  
  NaluEnv::self().naluOutputP0() << "Memory Overview: " << std::endl;
  NaluEnv::self().naluOutputP0() << "nalu memory: total (over all cores) current/high-water mark= "
                                 << std::setw(15) << convert_bytes(global_now[2])
                                 << std::setw(15) << convert_bytes(global_hwm[2])
                                 << std::endl;
  
  NaluEnv::self().naluOutputP0() << "nalu memory:   min (over all cores) current/high-water mark= "
                                 << std::setw(15) << convert_bytes(global_now[0])
                                 << std::setw(15) << convert_bytes(global_hwm[0])
                                 << std::endl;
  
  NaluEnv::self().naluOutputP0() << "nalu memory:   max (over all cores) current/high-water mark= "
                                  << std::setw(15) << convert_bytes(global_now[1])
                                  << std::setw(15) << convert_bytes(global_hwm[1])
                                  << std::endl;
}

//--------------------------------------------------------------------------
//-------- convert_bytes ---------------------------------------------------
//--------------------------------------------------------------------------
std::string 
Realm::convert_bytes(double bytes)
{
  const double K = 1024;
  const double M = K*1024;
  const double G = M*1024;

  std::ostringstream out;
  if (bytes < K) {
    out << bytes << " B";
  } else if (bytes < M) {
    bytes /= K;
    out << bytes << " K";
  } else if (bytes < G) {
    bytes /= M;
    out << bytes << " M";
  } else {
    bytes /= G;
    out << bytes << " G";
  }
  return out.str();
}

//--------------------------------------------------------------------------
//-------- initialize -----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::initialize()
{

  // example of how to check for something special
  sample_look_ahead();

  // initialize adaptivity - note: must be done before field registration
  setup_adaptivity();

  // field registration
  setup_nodal_fields();
  setup_edge_fields();
  setup_element_fields();

  // property maps and evaluation algorithms
  setup_property();

  // interior algorithm creation
  setup_interior_algorithms();

  // create boundary conditions
  setup_bc();

  // post processing algorithm creation
  setup_post_processing_algorithms();

  // create initial conditions
  setup_initial_conditions();

  // set global variables that have not yet been set
  initialize_global_variables();

  // Populate_mesh fills in the entities (nodes/elements/etc) and
  // connectivities, but no field-data. Field-data is not allocated yet.
  ioBroker_->populate_mesh();

  // If we want to create all internal edges, we want to do it before
  // field-data is allocated because that allows better performance in
  // the create-edges code.
  if (realmUsesEdges_ )
    create_edges();

  // output entity counts including max/min
  if ( provideEntityCount_ )
    provide_entity_count();

  // Now the mesh is fully populated, so we're ready to populate
  // field-data including coordinates, and attributes and/or distribution factors
  // if those exist on the input mesh file.
  ioBroker_->populate_field_data();

  // manage NaluGlobalId for linear system
  set_global_id();

  // check that all bcs are covering exposed surfaces
  if ( checkForMissingBcs_ )
    enforce_bc_on_exposed_faces();

  // output and restart files
  create_output_mesh();
  create_restart_mesh();

  // variables that may come from the initial mesh
  input_variables_from_mesh();

  populate_boundary_data();

  if ( solutionOptions_->meshMotion_ )
    process_mesh_motion();

  if ( hasPeriodic_ )
    periodicManager_->build_constraints();

  if ( has_mesh_deformation() )
    init_current_coordinates();

  compute_geometry();

  if ( hasContact_ )
    initialize_contact();

  if ( hasNonConformal_ )
    initialize_non_conformal();

  if ( hasOverset_ )
    initialize_overset();

  compute_l2_scaling();

  equationSystems_.initialize();

  // check job run size after mesh creation, linear system initialization
  check_job(false);
}

//--------------------------------------------------------------------------
//-------- sample_look_ahead -----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::sample_look_ahead()
{
  // example of how we might query all the nodes
  bool demo_parser_look_ahead=false;
  if (demo_parser_look_ahead) {
    const YAML::Node& root_node = root()->m_root_node;
    std::vector<const YAML::Node *> found_nodes;
    NaluParsingHelper::find_nodes_given_key("wall_boundary_condition", root_node, found_nodes);
    for (unsigned ii=0; ii < found_nodes.size(); ++ii)
    {
      const YAML::Node *val = found_nodes[ii]->FindValue("special_condition_XXX");
      if (val)
      {
        // do something...
        std::string out;
        *val >> out;
        NaluEnv::self().naluOutputP0() << "found special_condition_XXX = " << out << std::endl;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::load(const YAML::Node & node)
{

  //======================================
  // realm commands first
  //======================================

  node["name"] >> name_;
  node["mesh"] >> inputDBName_;

  // provide a high level banner
  NaluEnv::self().naluOutputP0() << std::endl;
  NaluEnv::self().naluOutputP0() << "Realm Options Review: " << name_ << std::endl;
  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;

  get_if_present(node, "estimate_memory_only", estimateMemoryOnly_, false);
  get_if_present(node, "available_memory_per_core_GB", availableMemoryPerCoreGB_, 0.0);

  // exposed bc check
  get_if_present(node, "check_for_missing_bcs", checkForMissingBcs_, checkForMissingBcs_);

  // entity count
  get_if_present(node, "provide_entity_count", provideEntityCount_, provideEntityCount_);

  // determine if edges are required and whether or not stk handles this
  node["use_edges"] >> realmUsesEdges_;

  // let everyone know about core algorithm
  if ( realmUsesEdges_ ) {
    NaluEnv::self().naluOutputP0() << "Edge-based scheme will be activated" << std::endl;
  }
  else {
    NaluEnv::self().naluOutputP0() <<"Element-based scheme will be activated" << std::endl;
  }

  // how often is the realm solved..
  get_if_present(node, "solve_frequency", solveFrequency_, solveFrequency_);

  // automatic decomposition
  get_if_present(node, "automatic_decomposition_type", autoDecompType_, autoDecompType_);

  // activate aura
  get_if_present(node, "activate_aura", activateAura_, activateAura_);
  if ( activateAura_ )
    NaluEnv::self().naluOutputP0() << "Nalu will activate aura ghosting" << std::endl;
  else
    NaluEnv::self().naluOutputP0() << "Nalu will deactivate aura ghosting" << std::endl;

  // memory diagnostic
  get_if_present(node, "activate_memory_diagnostic", activateMemoryDiagnostic_, activateMemoryDiagnostic_);
  if ( activateMemoryDiagnostic_ )
    NaluEnv::self().naluOutputP0() << "Nalu will activate detailed memory pulse" << std::endl;
  
  // allow for inconsistent restart (fields are missing)
  get_if_present(node, "support_inconsistent_multi_state_restart", supportInconsistentRestart_, supportInconsistentRestart_);

  // time step control
  const bool dtOptional = true;
  const YAML::Node *y_time_step = expect_map(node,"time_step_control", dtOptional);
  if ( y_time_step ) {
    get_if_present(*y_time_step, "target_courant", targetCourant_, targetCourant_);
    get_if_present(*y_time_step, "time_step_change_factor", timeStepChangeFactor_, timeStepChangeFactor_);
  }

  //======================================
  // now other commands/actions
  //======================================

  // load output first so we can check for serializing i/o
  outputInfo_->load(node);
  if (root()->serializedIOGroupSize_ == 0)
  {
    // only set from input file if command-line didn't set it
    root()->setSerializedIOGroupSize(outputInfo_->serializedIOGroupSize_);
  }

  // solution options - loaded before create_mesh since we need to know if
  // adaptivity is on to create the proper MetaData
  solutionOptions_->load(node);

  // once we know the mesh name, we can open the meta data, and set spatial dimension
  create_mesh();
  spatialDimension_ = metaData_->spatial_dimension();

  // averaging
  averagingInfo_->load(node);

  // post processing
  postProcessingInfo_->load(node);

  // norms
  solutionNormPostProcessing_->load(node);

  // boundary, init, material and equation systems "load"
  NaluEnv::self().naluOutputP0() << std::endl;
  NaluEnv::self().naluOutputP0() << "Boundary Condition Review: " << std::endl;
  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;
  boundaryConditions_.load(node);
  NaluEnv::self().naluOutputP0() << std::endl;
  NaluEnv::self().naluOutputP0() << "Initial Condition Review:  " << std::endl;
  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;
  initialConditions_.load(node);
  NaluEnv::self().naluOutputP0() << std::endl;
  NaluEnv::self().naluOutputP0() << "Material Prop Review:      " << std::endl;
  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;
  materialPropertys_.load(node);
  NaluEnv::self().naluOutputP0() << std::endl;
  NaluEnv::self().naluOutputP0() << "EqSys/options Review:      " << std::endl;
  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;
  equationSystems_.load(node);

  // set number of nodes, check job run size
  check_job(true);

}

Simulation *Realm::root() { return parent()->root(); }
Simulation *Realm::root() const { return parent()->root(); }
Realms *Realm::parent() { return &realms_; }
Realms *Realm::parent() const { return &realms_; }


//--------------------------------------------------------------------------
//-------- setup_adaptivity ------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_adaptivity()
{
#if defined (NALU_USES_PERCEPT)

  if ((solutionOptions_->activateUniformRefinement_ ||
       (solutionOptions_->useAdapter_ && solutionOptions_->maxRefinementLevel_ > 0)) && NULL == adapter_) {
      adapter_ = new Adapter(*this);
      // get the part that holds the "parent" elements, that is, elements that have been refined
      //   and thus have child elements - we want to avoid them in each bucket loop
      stk::mesh::EntityRank part_ranks[] = {stk::topology::ELEMENT_RANK, metaData_->side_rank()};
      unsigned nranks = 2;
      for (unsigned irank=0; irank < nranks; ++irank) {
        std::ostringstream inactive_part_name;
        inactive_part_name << "refine_inactive_elements_part_" << static_cast<unsigned int>(part_ranks[irank]);
        //stk::mesh::Part* child_elements_part = m_eMesh.get_non_const_part(active_part_name);
        stk::mesh::Part* parent_elements_part = metaData_->get_part(inactive_part_name.str());
        if (!parent_elements_part)
          throw std::runtime_error("error - no parent_elements_part can be found");
        adapterSelector_[part_ranks[irank]] = !stk::mesh::Selector(*parent_elements_part);
      }
  }

  // fields
  if (solutionOptions_->useMarker_)
    {
      percept::RefineFieldType *refineField= &(metaData_->declare_field<percept::RefineFieldType>(stk::topology::ELEMENT_RANK, "refine_field"));
      stk::mesh::put_field(*refineField, metaData_->universal_part());
      percept::RefineFieldType *refineFieldOrig= &(metaData_->declare_field<percept::RefineFieldType>(stk::topology::ELEMENT_RANK, "refine_field_orig"));
      stk::mesh::put_field(*refineFieldOrig, metaData_->universal_part());
      percept::RefineLevelType& refine_level = metaData_->declare_field<percept::RefineLevelType>(stk::topology::ELEMENT_RANK, "refine_level");
      stk::mesh::put_field( refine_level , metaData_->universal_part());
    }
#endif
}

//--------------------------------------------------------------------------
//-------- setup_nodal_fields ----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_nodal_fields()
{
  // register global id and rank fields on all parts
  const stk::mesh::PartVector parts = metaData_->get_parts();
  for ( size_t ipart = 0; ipart < parts.size(); ++ipart ) {
    naluGlobalId_ = &(metaData_->declare_field<GlobalIdFieldType>(stk::topology::NODE_RANK, "nalu_global_id"));
    stk::mesh::put_field(*naluGlobalId_, *parts[ipart]);
  }

  // loop over all material props targets and register nodal fields
  std::vector<std::string> targetNames = materialPropertys_.targetNames_;
  equationSystems_.register_nodal_fields(targetNames);

  // check for norm nodal fields
  if ( NULL != solutionNormPostProcessing_ )
    solutionNormPostProcessing_->setup(targetNames);
}

//--------------------------------------------------------------------------
//-------- setup_edge_fields -----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_edge_fields()
{
  // loop over all material props targets and register edge fields
  std::vector<std::string> targetNames = materialPropertys_.targetNames_;
  equationSystems_.register_edge_fields(targetNames);
}
//--------------------------------------------------------------------------
//-------- setup_element_fields --------------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_element_fields()
{
  // loop over all material props targets and register element fields
  std::vector<std::string> targetNames = materialPropertys_.targetNames_;
  equationSystems_.register_element_fields(targetNames);
}

//--------------------------------------------------------------------------
//-------- setup_interior_algorithms ---------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_interior_algorithms()
{

  // create adaptivity error algorithm
  if ( solutionOptions_->activateAdaptivity_) {
    if ( NULL == errorIndicatorAlgDriver_ )
      errorIndicatorAlgDriver_ = new ErrorIndicatorAlgorithmDriver(*this);
  }
  // loop over all material props targets and register interior algs
  std::vector<std::string> targetNames = materialPropertys_.targetNames_;
  equationSystems_.register_interior_algorithm(targetNames);
}

//--------------------------------------------------------------------------
//-------- setup_post_processing_algorithms --------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_post_processing_algorithms()
{
  // get a pointer to the post processing data vector
  std::vector<PostProcessingData* > &ppDataVec = postProcessingInfo_->ppDataVec_;

  // iterate and set-up
  std::vector<PostProcessingData *>::const_iterator ii;
  for( ii=ppDataVec.begin(); ii!=ppDataVec.end(); ++ii ) {

    PostProcessingData &theData = *(*ii);
    // type
    std::string theType = theData.type_;
    NaluEnv::self().naluOutputP0() << "the post processing type is " << theType << std::endl;

    // output name
    std::string theFile = theData.outputFileName_;
    NaluEnv::self().naluOutputP0() << "the post processing file name: " << theFile << std::endl;

    // physics
    std::string thePhysics = theData.physics_;
    NaluEnv::self().naluOutputP0() << "the post processing physics name: " << thePhysics << std::endl;

    // target
    std::vector<std::string> targetNames = theData.targetNames_;
    for ( size_t in = 0; in < targetNames.size(); ++in)
      NaluEnv::self().naluOutputP0() << "Target name(s): " << targetNames[in] << std::endl;

    // params
    std::vector<double> parameters = theData.parameters_;
    for ( size_t in = 0; in < parameters.size(); ++in)
      NaluEnv::self().naluOutputP0() << "Parameters used are: " << parameters[in] << std::endl;

    // call through to the Eqsys
    if ( theType == "surface" ) {
      equationSystems_.register_surface_pp_algorithm(theData);
    }
    else {
      throw std::runtime_error("Post Processing Error: only  surface-based is supported");
    }

  }
}

//--------------------------------------------------------------------------
//-------- setup_bc --------------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_bc()
{
  // loop over all bcs and register
  for (size_t ibc = 0; ibc < boundaryConditions_.size(); ++ibc) {
    BoundaryCondition& bc = *boundaryConditions_[ibc];
    switch(bc.theBcType_) {
      case WALL_BC:
        equationSystems_.register_wall_bc(bc.targetName_, *reinterpret_cast<const WallBoundaryConditionData *>(&bc));
        break;
      case INFLOW_BC:
        equationSystems_.register_inflow_bc(bc.targetName_, *reinterpret_cast<const InflowBoundaryConditionData *>(&bc));
        break;
      case OPEN_BC:
        equationSystems_.register_open_bc(bc.targetName_, *reinterpret_cast<const OpenBoundaryConditionData *>(&bc));
        break;
      case CONTACT_BC:
        equationSystems_.register_contact_bc(bc.targetName_, *reinterpret_cast<const ContactBoundaryConditionData *>(&bc));
        break;
      case SYMMETRY_BC:
        equationSystems_.register_symmetry_bc(bc.targetName_, *reinterpret_cast<const SymmetryBoundaryConditionData *>(&bc));
        break;
      case PERIODIC_BC:
        equationSystems_.register_periodic_bc(
          (*reinterpret_cast<const PeriodicBoundaryConditionData *>(&bc)).masterSlave_.master_,
          (*reinterpret_cast<const PeriodicBoundaryConditionData *>(&bc)).masterSlave_.slave_,
          *reinterpret_cast<const PeriodicBoundaryConditionData *>(&bc));
        break;
      case NON_CONFORMAL_BC:
        equationSystems_.register_non_conformal_bc(*reinterpret_cast<const NonConformalBoundaryConditionData *>(&bc));
        break;
      case OVERSET_BC:
        equationSystems_.register_overset_bc(*reinterpret_cast<const OversetBoundaryConditionData *>(&bc));
        break;
      default:
        throw std::runtime_error("unknown bc");
    }
  }
}

//--------------------------------------------------------------------------
//-------- enforce_bc_on_exposed_faces  ------------------------------------
//--------------------------------------------------------------------------
void
Realm::enforce_bc_on_exposed_faces()
{

  // first, skin mesh and, therefore, populate
  stk::mesh::Selector activePart = metaData_->locally_owned_part() | metaData_->globally_shared_part();
  stk::mesh::PartVector partVec;
  partVec.push_back(exposedBoundaryPart_);
  stk::mesh::skin_mesh(*bulkData_, activePart, partVec);

  stk::mesh::Selector selectRule = stk::mesh::Selector(*exposedBoundaryPart_)
    & !stk::mesh::selectUnion(bcPartVec_);

  stk::mesh::BucketVector const& face_buckets = bulkData_->get_buckets(metaData_->side_rank(), selectRule);

  if (!face_buckets.empty()) {
    NaluEnv::self().naluOutputP0() << "Exposed surfaces found without a boundary condition applied" << std::endl;

    // proceed to show the problem faces
    for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
            ib != face_buckets.end() ; ++ib )
    {
      stk::mesh::Bucket & b = **ib ;
      const stk::mesh::Bucket::size_type length   = b.size();
      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        // report offending set of faces; okay if to P0
        NaluEnv::self().naluOutputP0() << "Face Id: " << bulkData_->identifier(b[k]) << " is not properly covered" << std::endl;
      }
    }
    throw std::runtime_error("Realm::Error: Please aply bc to problematic exposed surfaces ");
  }

}

//--------------------------------------------------------------------------
//-------- setup_initial_conditions ----------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_initial_conditions()
{
  // loop over all ics and register
  for (size_t j_ic = 0; j_ic < initialConditions_.size(); ++j_ic) {
    InitialCondition& initCond = *initialConditions_[j_ic];

    const std::vector<std::string> targetNames = initCond.targetNames_;

    for (size_t itarget=0; itarget < targetNames.size(); ++itarget) {

      // target need not be subsetted since nothing below will depend on topo
      stk::mesh::Part *targetPart = metaData_->get_part(targetNames[itarget]);

      switch(initCond.theIcType_) {

        case CONSTANT_UD:
        {
          const ConstantInitialConditionData& genIC = *reinterpret_cast<const ConstantInitialConditionData *>(&initCond);
          ThrowAssert(genIC.data_.size() == genIC.fieldNames_.size());
          for (size_t ifield = 0; ifield < genIC.fieldNames_.size(); ++ifield) {

            std::vector<double>  genSpec = genIC.data_[ifield];
            stk::mesh::FieldBase *field = stk::mesh::get_field_by_name(genIC.fieldNames_[ifield], *metaData_);
            ThrowAssert(field);
      
            stk::mesh::FieldBase *fieldWithState = ( field->number_of_states() > 1 )
              ? field->field_state(stk::mesh::StateNP1)
              : field->field_state(stk::mesh::StateNone);

            std::vector<double> userGen = genSpec;
            ConstantAuxFunction *theGenFunc = new ConstantAuxFunction(0, genSpec.size(), userGen);
            AuxFunctionAlgorithm *auxGen
              = new AuxFunctionAlgorithm( *this, targetPart,
                                          fieldWithState, theGenFunc, stk::topology::NODE_RANK);
            initCondAlg_.push_back(auxGen);

          }
        }
        break;

        case FUNCTION_UD:
        {
          const UserFunctionInitialConditionData& fcnIC = *reinterpret_cast<const UserFunctionInitialConditionData *>(&initCond);
          equationSystems_.register_initial_condition_fcn(targetPart, fcnIC);
        }
        break;

        case USER_SUB_UD:
          throw std::runtime_error("Realm::setup_initial_conditions: USER_SUB not supported: ");

        case UserDataType_END:
          break;

        default:
          NaluEnv::self().naluOutputP0() << "Realm::setup_initial_conditions: unknown type: " << initCond.theIcType_ << std::endl;
          throw std::runtime_error("Realm::setup_initial_conditions: unknown type:");
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- setup_property --------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_property()
{
  // loop over all target names
  const std::vector<std::string> targetNames = materialPropertys_.targetNames_;
  for (size_t itarget=0; itarget < targetNames.size(); ++itarget) {

    // target need not be subsetted since nothing below will depend on topo
    stk::mesh::Part *targetPart = metaData_->get_part(targetNames[itarget]);

    // loop over propertyMap
    std::map<PropertyIdentifier, ScalarFieldType *>::iterator ii;
    for ( ii=propertyMap_.begin(); ii!=propertyMap_.end(); ++ii ) {

      // extract property id and field pointer
      PropertyIdentifier thePropId = (*ii).first;
      ScalarFieldType *thePropField = (*ii).second;

      // find the material property data object
      MaterialPropertyData *matData = NULL;
      std::map<PropertyIdentifier, MaterialPropertyData*>::iterator itf =
        materialPropertys_.propertyDataMap_.find(thePropId);
      if ( itf == materialPropertys_.propertyDataMap_.end() ) {
        // will need to throw
        NaluEnv::self().naluOutputP0() << "issue with property: " << PropertyIdentifierNames[thePropId] << std::endl;
        throw std::runtime_error("Please add property specification ");
      }
      else {
        matData = (*itf).second;
      }

      switch( matData->type_) {

        case CONSTANT_MAT:
        {

          // for constant specific heat, proceed in specialty code; create the appropriate enthalpy evaluator
          if (thePropId == SPEC_HEAT_ID && needsEnthalpy_) {
            // extract reference temperature
            double tRef = 300.0;
            extract_universal_constant("reference_temperature", tRef, true);

            // set up evaluators required for all cases
            PropertyEvaluator *theCpPropEval = NULL;
            PropertyEvaluator *theEnthPropEval = NULL;

            // check for species-based cp
            if ( matData->cpConstMap_.size() > 0.0 ) {
              if ( uniformFlow_ ) {
                throw std::runtime_error("uniform flow cp should simply use the single-valued constant");
              }
              else {
                // props computed based on local mass fractions, however, constant per species k
                theCpPropEval = new SpecificHeatConstCpkPropertyEvaluator(
                    matData->cpConstMap_, *metaData_);
                theEnthPropEval = new EnthalpyConstCpkPropertyEvaluator(
                    matData->cpConstMap_, matData->hfConstMap_, *metaData_, tRef);
              }

              // create the algorithm to compute Cp; EnthalpyEqs manages h population, i.e., no alg required
              TemperaturePropAlgorithm *auxAlg
                = new TemperaturePropAlgorithm( *this, targetPart, thePropField, theCpPropEval);
              propertyAlg_.push_back(auxAlg);
              
            }
            else {
              // single constant value
              double specificHeatValue = matData->constValue_;
              theCpPropEval = new ConstantPropertyEvaluator(specificHeatValue);
              theEnthPropEval = new EnthalpyConstSpecHeatPropertyEvaluator(specificHeatValue, tRef);

              // set the default begin/end
              int theBegin = 0;
              int theEnd = 1;

              // create everything
              std::vector<double> userConstData(1);
              userConstData[0] = matData->constValue_;
              ConstantAuxFunction *theAuxFunc
                = new ConstantAuxFunction(theBegin, theEnd, userConstData);
              AuxFunctionAlgorithm *auxAlg
                = new AuxFunctionAlgorithm( *this, targetPart,
					    thePropField, theAuxFunc, stk::topology::NODE_RANK);
              propertyAlg_.push_back(auxAlg);

            }
            // push to prop eval
            materialPropertys_.propertyEvalMap_[SPEC_HEAT_ID] = theCpPropEval;
            materialPropertys_.propertyEvalMap_[ENTHALPY_ID]  = theEnthPropEval;
          }
          else {

            // set the default begin/end
            int theBegin = 0;
            int theEnd = 1;

            // create everything
            std::vector<double> userConstData(1);
            userConstData[0] = matData->constValue_;
            ConstantAuxFunction *theAuxFunc
              = new ConstantAuxFunction(theBegin, theEnd, userConstData);
            AuxFunctionAlgorithm *auxAlg
              = new AuxFunctionAlgorithm( *this, targetPart,
					  thePropField, theAuxFunc, stk::topology::NODE_RANK);
            propertyAlg_.push_back(auxAlg);

          }
        }
        break;

        case MIXFRAC_MAT:
        {
          // extract the mixture fraction field
          ScalarFieldType *mixFrac = metaData_->get_field<ScalarFieldType>(stk::topology::NODE_RANK, "mixture_fraction");

          // primary and secondary
          const double propPrim = matData->primary_;
          const double propSec = matData->secondary_;

          // density requires inverse weighting
          if ( DENSITY_ID == thePropId ) {
            // create the inverse mix frac property algorithm
            InversePropAlgorithm *auxAlg
              = new InversePropAlgorithm( *this, targetPart, thePropField, mixFrac, propPrim, propSec);
            propertyAlg_.push_back(auxAlg);
          }
          else {
            // all else need linear weighting
            LinearPropAlgorithm *auxAlg
              = new LinearPropAlgorithm( *this, targetPart, thePropField, mixFrac, propPrim, propSec);
            propertyAlg_.push_back(auxAlg);
          }
        }
        break;

        case POLYNOMIAL_MAT:
        {

          // switch on property id
          switch( thePropId ) {

            case VISCOSITY_ID:
            {
              PropertyEvaluator *viscPropEval = NULL;
              
              if ( isothermalFlow_ ) {

                // all props will use Tref
                double tRef = 0.0;
                extract_universal_constant("reference_temperature", tRef, true);

                if ( uniformFlow_ ) {
                  // props computed based on YkRef and Tref
                  throw std::runtime_error("Realm::setup_property: Sorry, polynomial visc Ykref and Tref is not supported");
                }
                else {
                  // props computed based on Yk and Tref
                  viscPropEval = new SutherlandsYkTrefPropertyEvaluator(
                    matData->polynomialCoeffsMap_, *metaData_, tRef);
                }
                // create the GenericPropAlgorithm; push it back
                GenericPropAlgorithm *auxAlg
                  = new GenericPropAlgorithm( *this, targetPart, thePropField, viscPropEval);
                propertyAlg_.push_back(auxAlg);
              }
              else {
                // all props will use [transported] T
                if ( uniformFlow_ ) {
                  // props computed based on YkRef and T
                  viscPropEval = new SutherlandsPropertyEvaluator(
                    materialPropertys_.referencePropertyDataMap_, matData->polynomialCoeffsMap_ );
                }
                else {
                  // props computed based on Yk and T
                  viscPropEval = new SutherlandsYkPropertyEvaluator(
                    matData->polynomialCoeffsMap_, *metaData_);
                }
                // create the TemperaturePropAlgorithm; push it back
                TemperaturePropAlgorithm *auxAlg
                  = new TemperaturePropAlgorithm( *this, targetPart, thePropField, viscPropEval);
                propertyAlg_.push_back(auxAlg);
              }

              // create the property alg and push to evalmap
              materialPropertys_.propertyEvalMap_[thePropId] = viscPropEval;

            }
            break;

            case ENTHALPY_ID:
            {
              NaluEnv::self().naluOutputP0() << "Enthalpy specification is not required as Cp is sufficient";
            }
            break;

            case SPEC_HEAT_ID:
            {
              // R
              double universalR = 8314.4621;
              extract_universal_constant("universal_gas_constant", universalR, true);

              // create the property alg and push to evalmap
              PropertyEvaluator *theCpPropEval = NULL;
              PropertyEvaluator *theEnthPropEval = NULL;
              if ( uniformFlow_ ) {
                // props computed based on reference values
                theCpPropEval = new SpecificHeatPropertyEvaluator(
                    materialPropertys_.referencePropertyDataMap_, matData->lowPolynomialCoeffsMap_,
                    matData->highPolynomialCoeffsMap_, universalR);
                theEnthPropEval = new EnthalpyPropertyEvaluator(
                    materialPropertys_.referencePropertyDataMap_, matData->lowPolynomialCoeffsMap_,
                    matData->highPolynomialCoeffsMap_, universalR);
              }
              else {
                // props computed based on transported Yk values
                theCpPropEval = new SpecificHeatTYkPropertyEvaluator(
                    materialPropertys_.referencePropertyDataMap_, matData->lowPolynomialCoeffsMap_,
                    matData->highPolynomialCoeffsMap_, universalR, *metaData_);
                theEnthPropEval = new EnthalpyTYkPropertyEvaluator(
                   materialPropertys_.referencePropertyDataMap_, matData->lowPolynomialCoeffsMap_,
                   matData->highPolynomialCoeffsMap_, universalR, *metaData_);
              }

              // create the algorithm to compute Cp; EnthalpyEqs manages h population, i.e., no alg required
              TemperaturePropAlgorithm *auxAlg
                = new TemperaturePropAlgorithm( *this, targetPart, thePropField, theCpPropEval);
              propertyAlg_.push_back(auxAlg);

              // set property maps...
              materialPropertys_.propertyEvalMap_[SPEC_HEAT_ID] = theCpPropEval;
              materialPropertys_.propertyEvalMap_[ENTHALPY_ID] = theEnthPropEval;
            }
            break;

            default:
              throw std::runtime_error("Realm::setup_property: polynomial supports mu, Cp and h:");
          }
        }
        break;

        case IDEAL_GAS_T_MAT: case IDEAL_GAS_T_P_MAT:
        {
          if ( DENSITY_ID == thePropId ) {
        
            // everyone will require R
            double universalR = 8314.4621;
            extract_universal_constant("universal_gas_constant", universalR, true);
            
            // create the property evaluator
            PropertyEvaluator *rhoPropEval = NULL;
            if ( uniformFlow_ ) {
              // load mw and reference species
              std::vector<std::pair<double, double> > mwMassFracVec;
              std::map<std::string, ReferencePropertyData*>::const_iterator itrp;
              for ( itrp = materialPropertys_.referencePropertyDataMap_.begin();
                    itrp!= materialPropertys_.referencePropertyDataMap_.end(); ++itrp) {
                ReferencePropertyData *propData = (*itrp).second;
                std::pair<double,double> thePair;
                thePair = std::make_pair(propData->mw_,propData->massFraction_);
                mwMassFracVec.push_back(thePair);
              }
              if ( IDEAL_GAS_T_MAT == matData->type_ ) {
                double pRef = 101325.0;
                extract_universal_constant("reference_pressure", pRef, true);
                rhoPropEval = new IdealGasTPropertyEvaluator(pRef, universalR, mwMassFracVec);
              }
              else {
                rhoPropEval = new IdealGasTPPropertyEvaluator(universalR, mwMassFracVec, *metaData_);
              }
            }
            else {
              // load mw
              std::vector<double> mwVec;
              std::map<std::string, ReferencePropertyData*>::const_iterator itrp;
              for ( itrp = materialPropertys_.referencePropertyDataMap_.begin();
                    itrp!= materialPropertys_.referencePropertyDataMap_.end(); ++itrp) {
                ReferencePropertyData *propData = (*itrp).second;
                mwVec.push_back(propData->mw_);
              }
              if ( IDEAL_GAS_T_MAT == matData->type_ ) {
                double pRef = 101325.0;
                extract_universal_constant("reference_pressure", pRef, true);
                rhoPropEval = new IdealGasTYkPropertyEvaluator(pRef, universalR, mwVec, *metaData_);
              }
              else {
                throw std::runtime_error("Realm::setup_property: ideal_gas_tp only supported for uniform flow:");
              }
            }

            // push back property evaluator to map
            materialPropertys_.propertyEvalMap_[thePropId] = rhoPropEval;

            // create the property algorithm
            TemperaturePropAlgorithm *auxAlg
              = new TemperaturePropAlgorithm( *this, targetPart, thePropField, rhoPropEval);
            propertyAlg_.push_back(auxAlg);
          }
          else {
            throw std::runtime_error("Realm::setup_property: ideal_gas_t only supported for density:");
          }

        }
        break;

        case IDEAL_GAS_YK_MAT:
        {
          if ( DENSITY_ID == thePropId ) {
        
            // pRef, tRef and R
            double pRef = 101325.0;
            double tRef = 300.0;
            double universalR = 8314.4621;
            extract_universal_constant("reference_pressure", pRef, true);
            extract_universal_constant("reference_temperature", tRef, true);
            extract_universal_constant("universal_gas_constant", universalR, true);

            // load mw
            std::vector<double> mwVec;
            std::map<std::string, ReferencePropertyData*>::const_iterator itrp;
            for ( itrp = materialPropertys_.referencePropertyDataMap_.begin();
                  itrp!= materialPropertys_.referencePropertyDataMap_.end(); ++itrp) {
              ReferencePropertyData *propData = (*itrp).second;
              mwVec.push_back(propData->mw_);
            }

            // create the property evaluator
            PropertyEvaluator *rhoPropEval = new IdealGasYkPropertyEvaluator(pRef, tRef, universalR, mwVec, *metaData_);

            // push back property evaluator to map
            materialPropertys_.propertyEvalMap_[thePropId] = rhoPropEval;

            // create the property algorithm
            GenericPropAlgorithm *auxAlg
              = new GenericPropAlgorithm( *this, targetPart, thePropField, rhoPropEval);
            propertyAlg_.push_back(auxAlg);

          }
          else {
            throw std::runtime_error("Realm::setup_property: ideal_gas_yk only supported for density:");
          }
        }
        break;

        case GEOMETRIC_MAT:
        {
          // propery is a function of inverse dual volume
          InverseDualVolumePropAlgorithm *auxAlg
            = new InverseDualVolumePropAlgorithm( *this, targetPart, thePropField);
          propertyAlg_.push_back(auxAlg);
        }
        break;

        case HDF5_TABLE_MAT:
        {
	  if ( HDF5ptr_ == NULL ) {
	    HDF5ptr_ = new HDF5FilePtr( materialPropertys_.propertyTableName_ );
	  }

 	  // create the new TablePropAlgorithm that knows how to read from HDF5 file
 	  HDF5TablePropAlgorithm * auxAlg = new HDF5TablePropAlgorithm(*this, 
								       targetPart, 
								       HDF5ptr_->get_H5IO(),
								       thePropField, 
								       matData->tablePropName_, 
								       matData->indVarName_, 
								       matData->indVarTableName_,
								       *metaData_ );
          propertyAlg_.push_back(auxAlg);

	  NaluEnv::self().naluOutputP0() << "With " << matData->tablePropName_ << " also read table for auxVarName " <<matData->auxVarName_  << std::endl;
	  
	  //TODO : need to make auxVarName_ and tableAuxVarName_ into vectors and loop over them to create a set of new auxVar's and algorithms

          // auxVariable	  
          std::string auxVarName = matData->auxVarName_;
          if ( "na" != auxVarName ) {
            // register and put the field; assume a scalar for now; species extraction will complicate the matter
            ScalarFieldType *auxVar =  &(metaData_->declare_field<ScalarFieldType>(stk::topology::NODE_RANK, auxVarName));
            stk::mesh::put_field(*auxVar, *targetPart);
            // create the algorithm to populate it from an HDF5 file
	    HDF5TablePropAlgorithm * auxVarAlg = new HDF5TablePropAlgorithm(*this, 
									 targetPart, 
									 HDF5ptr_->get_H5IO(),
									 auxVar, 
									 matData->tableAuxVarName_, 
									 matData->indVarName_, 
									 matData->indVarTableName_,
									 *metaData_ );
            propertyAlg_.push_back(auxVarAlg);
          }

	}
	break;

      case GENERIC: 
        { 
          // default property evaluator
          PropertyEvaluator *propEval = NULL;
          
          // extract the property evaluator name
          std::string propEvalName = matData->genericPropertyEvaluatorName_;

          if ( propEvalName == "water_viscosity_T" ) {
            propEval = new WaterViscosityTPropertyEvaluator(*metaData_);
          }
          else if ( propEvalName == "water_density_T" ) {
            propEval = new WaterDensityTPropertyEvaluator(*metaData_);
          }
          else if ( propEvalName == "water_specific_heat_T" ) {
            propEval = new WaterSpecHeatTPropertyEvaluator(*metaData_);
            // create the enthalpy prop evaluator and store
            WaterEnthalpyTPropertyEvaluator *theEnthPropEval 
              = new WaterEnthalpyTPropertyEvaluator(*metaData_);
            materialPropertys_.propertyEvalMap_[ENTHALPY_ID] = theEnthPropEval;
          }
          else if ( propEvalName == "water_thermal_conductivity_T" ) {
            propEval = new WaterThermalCondTPropertyEvaluator(*metaData_);
          }
          else {
            throw std::runtime_error("Realm::setup_property: unknown GENERIC type: " + propEvalName);
          }
          
          // for now, all of the above are TempPropAlgs; push it back
          TemperaturePropAlgorithm *auxAlg
            = new TemperaturePropAlgorithm( *this, targetPart, thePropField, propEval);
          propertyAlg_.push_back(auxAlg);

          // push back property evaluator to map
          materialPropertys_.propertyEvalMap_[thePropId] = propEval;
        }
        break;

        case MaterialPropertyType_END:
          break;

        default:
          throw std::runtime_error("Realm::setup_property: unknown type:");
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- extract_universal_constant --------------------------------------
//--------------------------------------------------------------------------
void
Realm::extract_universal_constant( 
  const std::string name, double &value, const bool useDefault)
{
  std::map<std::string, double >::iterator it
    = materialPropertys_.universalConstantMap_.find(name);
  if ( it == materialPropertys_.universalConstantMap_.end() ) {
    // not found
    if ( useDefault ) {
      NaluEnv::self().naluOutputP0() << "WARNING: Reference value for " << name << " not found " << " using " << value << std::endl;
    }
    else {
      throw std::runtime_error("Realm::setup_property: reference value not found: " + name);
    }
  }
  else {
    value = (*it).second;
  }
}

//--------------------------------------------------------------------------
//-------- initialize_global_variables -------------------------------------
//--------------------------------------------------------------------------
void
Realm::initialize_global_variables()
{
  // other variables created on the fly during Eqs registration
  const bool needInOutput = false;
  const bool needInRestart = true;
  globalParameters_.set_param("timeStepNm1", 1.0, needInOutput, needInRestart);
  globalParameters_.set_param("timeStepCount", 1, needInOutput, needInRestart);

  // consider pushing this parameter to some higher level design
  if ( averagingInfo_->processAveraging_ )
    globalParameters_.set_param("currentTimeFilter", 0.0, needInOutput, needInRestart);
}

//--------------------------------------------------------------------------
//-------- augment_property_map --------------------------------------------
//--------------------------------------------------------------------------
void
Realm::augment_property_map(
  PropertyIdentifier propID,
  ScalarFieldType *theField)
{
  propertyMap_[propID] = theField;
}

//--------------------------------------------------------------------------
//-------- makeSureNodesHaveValidTopology ----------------------------------
//--------------------------------------------------------------------------
void 
Realm::makeSureNodesHaveValidTopology()
{
  //To make sure nodes have valid topology, we have to make sure they are in a part that has NODE topology.
  //So first, let's obtain the node topology part:
  stk::mesh::Part& nodePart = bulkData_->mesh_meta_data().get_cell_topology_root_part(stk::mesh::get_cell_topology(stk::topology::NODE));
  stk::mesh::Selector nodesNotInNodePart = !nodePart & bulkData_->mesh_meta_data().locally_owned_part();

  //get all the nodes that are *NOT* in nodePart
  std::vector<stk::mesh::Entity> nodes_vector;
  stk::mesh::get_selected_entities(nodesNotInNodePart, bulkData_->buckets(stk::topology::NODE_RANK), nodes_vector);
  // now we require all nodes are in proper node part
  if (nodes_vector.size()) std::cout << "nodes_vector= " << nodes_vector.size() << std::endl;
  ThrowRequire(0 == nodes_vector.size());
}

//--------------------------------------------------------------------------
//-------- pre_timestep_work -----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::pre_timestep_work()
{

  if ( solutionOptions_->activateUniformRefinement_) {
    static stk::diag::Timer timerUniformRefine_("UniformRefinement", Simulation::rootTimer());
    static stk::diag::Timer timerComputeGeom_("ComputeGeom", timerUniformRefine_);
    static stk::diag::Timer timerCreateEdgesLocal_("CreateEdgesAfterAdapt", timerUniformRefine_);
    static stk::diag::Timer timerDeleteEdgesLocal_("DeleteEdgesBeforeAdapt", timerUniformRefine_);
    static stk::diag::Timer timerReInitLinSys_("ReInitLinSys", timerUniformRefine_);

    stk::diag::TimeBlock tbTimerUR_(timerUniformRefine_);
  
    for (unsigned i = 0; i < solutionOptions_->refineAt_.size(); ++i) {
      if ( solutionOptions_->refineAt_[i] == get_time_step_count() ||
           (solutionOptions_->refineAt_[i] == 0 && get_time_step_count() == 1) ) {

        NaluEnv::self().naluOutputP0() << "UniformRefinement: at step= " << get_time_step_count() << std::endl;

        if (realmUsesEdges_ ) {
          stk::diag::TimeBlock tbDeleteEdges_(timerDeleteEdgesLocal_);
          delete_edges();
        }

#if defined (NALU_USES_PERCEPT)
        adapter_->do_uniform_refine();
#endif

        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts( *bulkData_ , counts);

        NaluEnv::self().naluOutputP0() << "UniformRefine: after uniform refine, mesh has  "
                        << counts[0] << " nodes, "
                        << counts[1] << " edges, "
                        << counts[2] << " faces, "
                        << counts[3] << " elements" << std::endl;

        //call this temporary function to correct the part membership of any nodes that
        //are not in the node-topology part. Remove this function as soon as percept/adapt
        //is altered to make sure that all newly-created nodes are placed in that part.
        makeSureNodesHaveValidTopology();

        if (realmUsesEdges_ ) {
          stk::diag::TimeBlock tbCreateEdges_(timerCreateEdgesLocal_);
          create_edges();
        }

        {
          stk::diag::TimeBlock tbComputeGeom_(timerComputeGeom_);
          compute_geometry();
        }

        // now re-initialize linear system
        stk::diag::TimeBlock tbReInit_(timerReInitLinSys_);
        equationSystems_.reinitialize_linear_system();
        
      }
    }
  }

  // compute error indicator
  if ( solutionOptions_->activateAdaptivity_) {
    static stk::diag::Timer timerAdaptRealm_("AdaptRealm", Simulation::rootTimer());
    static stk::diag::Timer timerComputeGeom_("ComputeGeom", timerAdaptRealm_);
    static stk::diag::Timer timerCreateEdgesLocal_("CreateEdgesAfterAdapt", timerAdaptRealm_);
    static stk::diag::Timer timerDeleteEdgesLocal_("DeleteEdgesBeforeAdapt", timerAdaptRealm_);
    static stk::diag::Timer timerReInitLinSys_("ReInitLinSys", timerAdaptRealm_);

    stk::diag::TimeBlock tbTimerAdapt_(timerAdaptRealm_);
    double time = -stk::cpu_time();
    if ( process_adaptivity() ) {

#if defined (NALU_USES_PERCEPT)
      // mesh counts
      std::vector<size_t> counts;
      stk::mesh::comm_mesh_counts( *bulkData_ , counts);
      if (0 == numInitialElements_) {
        numInitialElements_ = counts[3];
      }

      errorIndicatorAlgDriver_->execute();

      if (solutionOptions_->useAdapter_ && solutionOptions_->maxRefinementLevel_)
        {
          NaluEnv::self().naluOutputP0() << "Adapt: running adapter..." << std::endl;
          NaluEnv::self().naluOutputP0() << "Adapt: before adapt, mesh has  "
                          << counts[0] << " nodes, "
                          << counts[1] << " edges, "
                          << counts[2] << " faces, "
                          << counts[3] << " elements" << std::endl;

#if USE_NALU_PERFORMANCE_TESTING_CALLGRIND
          CALLGRIND_START_INSTRUMENTATION;
          CALLGRIND_TOGGLE_COLLECT;
#endif

          // delete edges first
          if (realmUsesEdges_ ) {
            stk::diag::TimeBlock tbDeleteEdges_(timerDeleteEdgesLocal_);
            delete_edges();
          }

          adapter_->do_adapt(ADAPT_REFINE);

          counts.resize(0);
          stk::mesh::comm_mesh_counts( *bulkData_ , counts);

          NaluEnv::self().naluOutputP0() << "Adapt: after refine, mesh has  "
                          << counts[0] << " nodes, "
                          << counts[1] << " edges, "
                          << counts[2] << " faces, "
                          << counts[3] << " elements" << std::endl;

          // refinement destroys/recreates some elements, so we lose initial marking on some elements,
          //   thus, re-run the error indicator/marker before the unrefine stage.
          errorIndicatorAlgDriver_->execute();

          adapter_->do_adapt(ADAPT_UNREFINE);

          //call this temporary function to correct the part membership of any nodes that
          //are not in the node-topology part. Remove this function as soon as percept/adapt
          //is altered to make sure that all newly-created nodes are placed in that part.
          makeSureNodesHaveValidTopology();

#if USE_NALU_PERFORMANCE_TESTING_CALLGRIND
  CALLGRIND_TOGGLE_COLLECT;
  CALLGRIND_STOP_INSTRUMENTATION;
#endif

          counts.resize(0);
          stk::mesh::comm_mesh_counts( *bulkData_ , counts);

          NaluEnv::self().naluOutputP0() << "Adapt: after unrefine, mesh has  "
                          << counts[0] << " nodes, "
                          << counts[1] << " edges, "
                          << counts[2] << " faces, "
                          << counts[3] << " elements" << std::endl;


          if (realmUsesEdges_ ) {
            stk::diag::TimeBlock tbCreateEdges_(timerCreateEdgesLocal_);
            create_edges();
          }

          {
            stk::diag::TimeBlock tbComputeGeom_(timerComputeGeom_);
            compute_geometry();
          }

          // now re-initialize linear system
          stk::diag::TimeBlock tbReInit_(timerReInitLinSys_);
          equationSystems_.reinitialize_linear_system();
          
          // process speciality methods for adaptivity
          NaluEnv::self().naluOutputP0() << std::endl;
          NaluEnv::self().naluOutputP0() << "Post Adapt Work:" << std::endl;
          NaluEnv::self().naluOutputP0() <<"===========================" << std::endl;
          equationSystems_.post_adapt_work();
          NaluEnv::self().naluOutputP0() <<"===========================" << std::endl;
          NaluEnv::self().naluOutputP0() << std::endl;

          outputInfo_->meshAdapted_ = true;
        }
#endif
    }
    time += stk::cpu_time();
    timerAdapt_ += time;
  }

  // check for mesh motion
  if ( solutionOptions_->meshMotion_ ) {

    process_mesh_motion();
    compute_geometry();

    // check for contact
    if ( hasContact_ )
      initialize_contact();

    // and non-conformal algorithm
    if ( hasNonConformal_ )
      initialize_non_conformal();

    // and overset algorithm
    if ( hasOverset_ )
      initialize_overset();

    // now re-initialize linear system
    equationSystems_.reinitialize_linear_system();

  }

  // deal with non-topology changes, however, moving mesh
  if ( has_mesh_deformation() ) {
    // extract target parts for this physics
    if ( solutionOptions_->externalMeshDeformation_ ) {
      std::vector<std::string> targetNames = materialPropertys_.targetNames_;
      for ( size_t itarget = 0; itarget < targetNames.size(); ++itarget ) {
        stk::mesh::Part *targetPart = metaData_->get_part(targetNames[itarget]);
        set_current_coordinates(targetPart);
      }
    }
    compute_geometry();
  }

  // ask the equation system to do some work
  equationSystems_.pre_timestep_work();
}

//--------------------------------------------------------------------------
//-------- evaluate_properties ---------------------------------------------
//--------------------------------------------------------------------------
void
Realm::evaluate_properties()
{
  double start_time = stk::cpu_time();
  for ( size_t k = 0; k < propertyAlg_.size(); ++k ) {
    propertyAlg_[k]->execute();
  }
  equationSystems_.evaluate_properties();
  double end_time = stk::cpu_time();
  timerPropertyEval_ += (end_time - start_time);
}

//--------------------------------------------------------------------------
//-------- advance_time_step -----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::advance_time_step()
{
  // leave if we do not need to solve
  const int timeStepCount = get_time_step_count();
  const bool advanceMe = (timeStepCount % solveFrequency_ ) == 0 ? true : false;
  if ( !advanceMe )
    return;

  NaluEnv::self().naluOutputP0() << "NLI"
                  << std::setw(8) << std::right << "Name"
                  << std::setw(22) << std::right << "Linear Iter"
                  << std::setw(16) << std::right << "Linear Res"
                  << std::setw(16) << std::right << "NLinear Res"
                  << std::setw(14) << std::right << "Scaled NLR" << std::endl;

  NaluEnv::self().naluOutputP0() << "---"
                  << std::setw(8) << std::right << "----"
                  << std::setw(22) << std::right << "-----------"
                  << std::setw(16) << std::right << "----------"
                  << std::setw(16) << std::right << "-----------"
                  << std::setw(14) << std::right << "----------" << std::endl;

  // evaluate new geometry based on latest mesh motion geometry state (provided that external is active)
  if ( solutionOptions_->externalMeshDeformation_ )
    compute_geometry();

  // evaluate properties based on latest state including boundary and and possible xfer
  evaluate_properties();

  // compute velocity relative to mesh
  compute_vrtm();

  const int numNonLinearIterations = equationSystems_.maxIterations_;
  for ( int i = 0; i < numNonLinearIterations; ++i ) {
    currentNonlinearIteration_ = i+1;
    NaluEnv::self().naluOutputP0()
      << currentNonlinearIteration_
      << "/" << numNonLinearIterations
      << std::setw(29) << std::right << "Equation System Iteration" << std::endl;

    const bool isConverged = equationSystems_.solve_and_update();

    // evaluate properties based on latest np1 solution
    evaluate_properties();

    if ( isConverged ) {
      NaluEnv::self().naluOutputP0() << "norm convergence criteria met for all equation systems: " << std::endl;
      NaluEnv::self().naluOutputP0() << "max scaled norm is: " << equationSystems_.provide_system_norm() << std::endl;
      break;
    }
  }

}

//--------------------------------------------------------------------------
//-------- output_converged_results ----------------------------------------
//--------------------------------------------------------------------------
void
Realm::output_converged_results()
{
  provide_output();
  provide_restart_output();
}

//--------------------------------------------------------------------------
//-------- compute_adaptive_time_step --------------------------------------
//--------------------------------------------------------------------------
double
Realm::compute_adaptive_time_step()
{
  // extract current time
  const double dtN = get_time_step();

  // ratio of how off we are
  const double factorOff = targetCourant_/maxCourant_;

  // scaling for dt and candidate
  const double dtScaling = ( targetCourant_ < maxCourant_ )
    ? std::max(factorOff, 1.0/timeStepChangeFactor_)
    : std::min(factorOff, timeStepChangeFactor_);
  const double candidateDt = dtN*dtScaling;

  return candidateDt;
}

//--------------------------------------------------------------------------
//-------- commit ----------------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::commit()
{
  //====================================================
  // Commit the meta data
  //====================================================
  metaData_->commit();
}

//--------------------------------------------------------------------------
//-------- create_mesh() ---------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::create_mesh()
{
  double start_time = stk::cpu_time();

  stk::ParallelMachine pm = NaluEnv::self().parallel_comm();
  
  // news for mesh constructs
  metaData_ = new stk::mesh::MetaData();
  bulkData_ = new stk::mesh::BulkData(*metaData_, pm, activateAura_ ? stk::mesh::BulkData::AUTO_AURA : stk::mesh::BulkData::NO_AUTO_AURA);
  ioBroker_ = new stk::io::StkMeshIoBroker( pm );
  ioBroker_->set_bulk_data(*bulkData_);

  // allow for automatic decomposition
  if (autoDecompType_ != "None") 
    ioBroker_->property_add(Ioss::Property("DECOMPOSITION_METHOD", autoDecompType_));
  
  // for adaptivity we need an additional rank to store parent/child relations
  if (solutionOptions_->useAdapter_ || solutionOptions_->activateUniformRefinement_) {
    std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
    entity_rank_names.push_back("FAMILY_TREE");
    ioBroker_->set_rank_name_vector(entity_rank_names);
  }

  // Initialize meta data (from exodus file); can possibly be a restart file..
  ioBroker_->add_mesh_database( inputDBName_,
      restarted_simulation() ? stk::io::READ_RESTART : stk::io::READ_MESH );
  ioBroker_->create_input_mesh();

  // declare an exposed part for later bc coverage check
  if ( checkForMissingBcs_ ) {
    exposedBoundaryPart_ = &metaData_->declare_part("exposed_boundary_part",metaData_->side_rank());
  }

  // declare a part to hold new edges
  if (realmUsesEdges_) {
    edgesPart_ = &metaData_->declare_part("create_edges_part", stk::topology::EDGE_RANK);
  }

  const double end_time = stk::cpu_time();

  // set mesh reading
  timerReadMesh_ = (end_time - start_time);

}

//--------------------------------------------------------------------------
//-------- create_output_mesh() --------------------------------------------
//--------------------------------------------------------------------------
void
Realm::create_output_mesh()
{
  // exodus output file creation
  if (outputInfo_->hasOutputBlock_ ) {

    if (outputInfo_->outputFreq_ == 0)
      return;

    // if we are adapting, skip when no I/O happens before first adapt step
    if (solutionOptions_->useAdapter_ && outputInfo_->meshAdapted_ == false &&
        solutionOptions_->adaptivityFrequency_ <= outputInfo_->outputFreq_) {
      return;
    }

    std::string oname =  outputInfo_->outputDBName_ ;
    if (solutionOptions_->useAdapter_ && solutionOptions_->maxRefinementLevel_) {
      static int fileid = 0;
      std::ostringstream fileid_ss;
      fileid_ss << std::setfill('0') << std::setw(4) << (fileid+1);
      if (fileid++ > 0) oname += "-s" + fileid_ss.str();
    }

    resultsFileIndex_ = ioBroker_->create_output_mesh( oname, stk::io::WRITE_RESULTS );

#if defined (NALU_USES_PERCEPT)

  if (solutionOptions_->useAdapter_ && outputInfo_->meshAdapted_) {
    stk::mesh::Selector selectRule =
      metaData_->locally_owned_part() |
      metaData_->globally_shared_part();

    activePartForIO_ = Teuchos::rcp(new stk::mesh::Selector(percept::make_active_part_selector(*metaData_, selectRule)));

    ioBroker_->set_subset_selector(resultsFileIndex_, activePartForIO_);
  }

#endif

    // Tell stk_io how to output element block nodal fields:
    // if 'true' passed to function, then output them as nodeset fields;
    // if 'false', then output as nodal fields (on all nodes of the mesh, zero-filled)
    // The option is provided since some post-processing/visualization codes do not
    // correctly handle nodeset fields.
    ioBroker_->use_nodeset_for_part_nodes_fields(resultsFileIndex_, outputInfo_->outputNodeSet_);

    // FIXME: add_field can take user-defined output name, not just varName
    for ( std::set<std::string>::iterator itorSet = outputInfo_->outputFieldNameSet_.begin();
        itorSet != outputInfo_->outputFieldNameSet_.end(); ++itorSet ) {
      std::string varName = *itorSet;
      stk::mesh::FieldBase *theField = stk::mesh::get_field_by_name(varName, *metaData_);
      if ( NULL == theField ) {
        NaluEnv::self().naluOutputP0() << " Sorry, no field by the name " << varName << std::endl;
      }
      else {
        // 'varName' is the name that will be written to the database
        // For now, just using the name of the stk field
        ioBroker_->add_field(resultsFileIndex_, *theField, varName);
      }
    }

    // reset this flag
    outputInfo_->meshAdapted_ = false;
  }
}

//--------------------------------------------------------------------------
//-------- create_restart_mesh() --------------------------------------------
//--------------------------------------------------------------------------
void
Realm::create_restart_mesh()
{
  // exodus restart file creation
  if (outputInfo_->hasRestartBlock_ ) {

    if (outputInfo_->restartFreq_ == 0)
      return;

    restartFileIndex_ = ioBroker_->create_output_mesh(outputInfo_->restartDBName_, stk::io::WRITE_RESTART);

    // loop over restart variable field names supplied by Eqs
    for ( std::set<std::string>::iterator itorSet = outputInfo_->restartFieldNameSet_.begin();
        itorSet != outputInfo_->restartFieldNameSet_.end(); ++itorSet ) {
      std::string varName = *itorSet;
      stk::mesh::FieldBase *theField = stk::mesh::get_field_by_name(varName,*metaData_);
      if ( NULL == theField ) {
        NaluEnv::self().naluOutputP0() << " Sorry, no field by the name " << varName << std::endl;
      }
      else {
        // add the field for a restart output
        ioBroker_->add_field(restartFileIndex_, *theField, varName);
        // if this is a restarted simulation, we will need input
        if ( restarted_simulation() )
          ioBroker_->add_input_field(stk::io::MeshField(*theField, varName));
      }
    }

    // now global params
    stk::util::ParameterMapType::const_iterator i = globalParameters_.begin();
    stk::util::ParameterMapType::const_iterator iend = globalParameters_.end();
    for (; i != iend; ++i) {
      std::string parameterName = (*i).first;
      stk::util::Parameter parameter = (*i).second;
      if(parameter.toRestartFile) {
        ioBroker_->add_global(restartFileIndex_, parameterName, parameter.value, parameter.type);
      }
    }

    // set max size for restart data base
    ioBroker_->get_output_io_region(restartFileIndex_)->get_database()->set_cycle_count(outputInfo_->restartMaxDataBaseStepSize_);

  }

}

//--------------------------------------------------------------------------
//-------- input_variables_from_mesh() --------------------------------------------
//--------------------------------------------------------------------------
void
Realm::input_variables_from_mesh()
{
  // no variables from an input mesh if this is a restart
  if ( !restarted_simulation()) {
    std::map<std::string, std::string>::const_iterator iter;
    for ( iter = solutionOptions_->inputVarFromFileMap_.begin();
          iter != solutionOptions_->inputVarFromFileMap_.end(); ++iter) {

      std::string varName = iter->first;
      std::string userName = iter->second;

      stk::mesh::FieldBase *theField = stk::mesh::get_field_by_name(varName,*metaData_);
      if ( NULL == theField ) {
        NaluEnv::self().naluOutputP0() << " Sorry, no field by the name " << varName << std::endl;
      }
      else {
        ioBroker_->add_input_field(stk::io::MeshField(*theField, userName));
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- augment_output_variable_list() ----------------------------------
//--------------------------------------------------------------------------
void
Realm::augment_output_variable_list(
    const std::string fieldName)
{
  outputInfo_->outputFieldNameSet_.insert(fieldName);
}

//--------------------------------------------------------------------------
//-------- augment_restart_variable_list -----------------------------------
//--------------------------------------------------------------------------
void
Realm::augment_restart_variable_list(
  std::string restartFieldName)
{
  outputInfo_->restartFieldNameSet_.insert(restartFieldName);
}

//--------------------------------------------------------------------------
//-------- create_edges -----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::create_edges()
{
  static stk::diag::Timer timerCE_("CreateEdges", Simulation::rootTimer());
  stk::diag::TimeBlock tbCreateEdges_(timerCE_);

  double start_time = stk::cpu_time();
  if (solutionOptions_->useAdapter_ && solutionOptions_->maxRefinementLevel_ > 0 ) {
    stk::mesh::create_edges(*bulkData_, adapterSelector_[stk::topology::ELEMENT_RANK], edgesPart_);
  }
  else {
    stk::mesh::create_edges(*bulkData_, metaData_->universal_part(), edgesPart_);
  }
  double stop_time = stk::cpu_time();

  // timer close-out
  const double total_edge_time = stop_time - start_time;
  timerCreateEdges_ += total_edge_time;

}

//--------------------------------------------------------------------------
//-------- provide_entity_count() ------------------------------------------
//--------------------------------------------------------------------------
void
Realm::provide_entity_count() {

  std::vector<size_t> counts;
  std::vector<size_t> minCounts;
  std::vector<size_t> maxCounts;
  stk::mesh::comm_mesh_counts( *bulkData_ , counts, minCounts, maxCounts);

  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;
  NaluEnv::self().naluOutputP0() << "Realm::provide_entity_count:   " << std::endl
		  << "nodes,    " << counts[0] << " min/max: " << minCounts[0] << "/" << maxCounts[0] << std::endl
		  << "edges,    " << counts[1] << " min/max: " << minCounts[1] << "/" << maxCounts[1] << std::endl
		  << "faces,    " << counts[2] << " min/max: " << minCounts[2] << "/" << maxCounts[2] << std::endl
		  << "elements, " << counts[3] << " min/max: " << minCounts[3] << "/" << maxCounts[3] << std::endl;
  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;

}

//--------------------------------------------------------------------------
//-------- delete_edges -----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::delete_edges()
{
  if (debug()) {
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts( *bulkData_ , counts);

    NaluEnv::self().naluOutputP0() << "Realm::delete_edges: before delete_edges, mesh has  "
                    << counts[0] << " nodes, "
                    << counts[1] << " edges, "
                    << counts[2] << " faces, "
                    << counts[3] << " elements" << std::endl;
  }

  stk::mesh::BucketVector const& edge_buckets = bulkData_->get_buckets( stk::topology::EDGE_RANK,  *edgesPart_);
  std::vector<stk::mesh::Entity> edges;
  stk::mesh::get_selected_entities( *edgesPart_ , edge_buckets, edges);

  if (debug()) {
    size_t sz = edges.size(), g_sz=0;
    stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &sz, &g_sz, 1);
    NaluEnv::self().naluOutputP0() << "P[" << bulkData_->parallel_rank() << "] Realm::delete_edges: edge list local size= "
				   << sz << " global size= " << g_sz << std::endl;
  }

  // delete elem -> edge relations
  bulkData_->modification_begin();
  for (unsigned ii=0; ii < edges.size(); ++ii) {
    while (true) {

      if (!bulkData_->is_valid(edges[ii]))
        throw std::runtime_error("bad edge 1");

      unsigned num_elems = bulkData_->num_elements(edges[ii]);

      if (!num_elems)
        break;

      stk::mesh::Entity const * const edge_elems = bulkData_->begin_elements(edges[ii]);
      stk::mesh::ConnectivityOrdinal const* edge_elem_ordinals = bulkData_->begin_element_ordinals(edges[ii]);

      stk::mesh::Entity to_rel = edge_elems[0];
      stk::mesh::RelationIdentifier to_id = edge_elem_ordinals[0];

      bool del = bulkData_->destroy_relation( to_rel, edges[ii], to_id);
      if (!del)
        throw std::runtime_error("delete_edges failed to delete up relation");
    }

    if (3 == metaData_->spatial_dimension()) {
      while (true) {

        if (!bulkData_->is_valid(edges[ii]))
          throw std::runtime_error("bad edge 1");

        unsigned num_faces = bulkData_->num_faces(edges[ii]);

        if (!num_faces)
          break;

        stk::mesh::Entity const * const edge_faces = bulkData_->begin_faces(edges[ii]);
        stk::mesh::ConnectivityOrdinal const* edge_face_ordinals = bulkData_->begin_face_ordinals(edges[ii]);

        stk::mesh::Entity to_rel = edge_faces[0];
        stk::mesh::RelationIdentifier to_id = edge_face_ordinals[0];

        bool del = bulkData_->destroy_relation( to_rel, edges[ii], to_id);
        if (!del)
          throw std::runtime_error("delete_edges failed to delete up relation for face");
      }
    }
  }

  // now delete edges
  for (unsigned ii=0; ii < edges.size(); ++ii) {

    if (bulkData_->is_valid(edges[ii]) && bulkData_->bucket(edges[ii]).owned()) {
      if ( ! bulkData_->destroy_entity( edges[ii] ) ) {
        unsigned num_elems = bulkData_->num_elements(edges[ii]);
        NaluEnv::self().naluOutputP0() << "P[" << bulkData_->parallel_rank() << "] deleting edge num_elems= " << num_elems
				       << std::endl;

        stk::mesh::EntityRank topRank = stk::topology::ELEMENT_RANK;
        if (solutionOptions_->useAdapter_ && solutionOptions_->maxRefinementLevel_ > 0)
          ++topRank;
        for (stk::mesh::EntityRank irank = stk::topology::EDGE_RANK; irank <= topRank; ++irank) {
          unsigned nc = bulkData_->num_connectivity(edges[ii], irank);
          NaluEnv::self().naluOutputP0() << "P[" << bulkData_->parallel_rank() << "] deleting edge nc[" << irank << "]= " << nc
					 << std::endl;
        }

        throw std::runtime_error("delete_edges failed to delete edge");
      }
    }
  }
  bulkData_->modification_end();

  if (debug()) {
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts( *bulkData_ , counts);

    NaluEnv::self().naluOutputP0() << "P[" << bulkData_->parallel_rank() << "] "
                  << "Realm::delete_edges: after delete_edges, mesh has  "
                  << counts[0] << " nodes, "
                  << counts[1] << " edges, "
                  << counts[2] << " faces, "
                  << counts[3] << " elements" << std::endl;
  }
}

//--------------------------------------------------------------------------
//-------- initialize_contact ----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::initialize_contact()
{
  contactManager_->initialize();
}

//--------------------------------------------------------------------------
//-------- initialize_non_conformal ----------------------------------------
//--------------------------------------------------------------------------
void
Realm::initialize_non_conformal()
{
  nonConformalManager_->initialize();
}

//--------------------------------------------------------------------------
//-------- initialize_overset ----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::initialize_overset()
{
  oversetManager_->initialize();
}

//--------------------------------------------------------------------------
//-------- get_coordinates_name ---------------------------------------------
//--------------------------------------------------------------------------
std::string
Realm::get_coordinates_name()
{
  return ( (solutionOptions_->meshMotion_ | solutionOptions_->meshDeformation_ | solutionOptions_->externalMeshDeformation_) 
           ? "current_coordinates" : "coordinates");
}

//--------------------------------------------------------------------------
//-------- has_mesh_motion -------------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::has_mesh_motion()
{
  return solutionOptions_->meshMotion_;
}

//--------------------------------------------------------------------------
//-------- has_mesh_deformation --------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::has_mesh_deformation()
{
  return solutionOptions_->externalMeshDeformation_ | solutionOptions_->meshDeformation_;
}

//--------------------------------------------------------------------------
//-------- does_mesh_move --------------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::does_mesh_move()
{
  return has_mesh_motion() | has_mesh_deformation();
}

//--------------------------------------------------------------------------
//-------- has_non_matching_boundary_face_alg ------------------------------
//--------------------------------------------------------------------------
bool
Realm::has_non_matching_boundary_face_alg()
{
  return hasContact_ | hasNonConformal_ | hasOverset_; 
}

//--------------------------------------------------------------------------
//-------- query_for_overset -----------------------------------------------
//--------------------------------------------------------------------------
bool 
Realm::query_for_overset() 
{
  for (size_t ibc = 0; ibc < boundaryConditions_.size(); ++ibc) {
    BoundaryCondition& bc = *boundaryConditions_[ibc];
    switch(bc.theBcType_) {
    case OVERSET_BC:
      hasOverset_ = true;
      break;
    default:
      hasOverset_ = false;
    }
  }
  return hasOverset_;
}

//--------------------------------------------------------------------------
//-------- process_mesh_motion ---------------------------------------------
//--------------------------------------------------------------------------
void
Realm::process_mesh_motion()
{
  // extract parameters; allows for omega to change...
  std::map<std::string, std::pair<std::vector<std::string>, double> >::const_iterator iter;
  for ( iter = solutionOptions_->meshMotionMap_.begin();
        iter != solutionOptions_->meshMotionMap_.end(); ++iter) {

    std::pair<std::vector<std::string>, double> thePair = iter->second;

    std::vector<std::string> theVector = thePair.first;
    const double theOmega = thePair.second;
    for (size_t k = 0; k < theVector.size(); ++k ) {

      stk::mesh::Part *targetPart = metaData_->get_part(theVector[k]);

      if ( NULL == targetPart ) {
        throw std::runtime_error("Sorry, no part name found" + theVector[k]);
      }
      else {
        set_omega(targetPart, theOmega);
        set_current_displacement(targetPart);
        set_current_coordinates(targetPart);
        set_mesh_velocity(targetPart);
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- set_omega -------------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::set_omega(
  stk::mesh::Part *targetPart,
  double scalarOmega)
{
  ScalarFieldType *omega = metaData_->get_field<ScalarFieldType>(stk::topology::NODE_RANK, "omega");

  stk::mesh::Selector s_all_nodes = stk::mesh::Selector(*targetPart);

  stk::mesh::BucketVector const& node_buckets = bulkData_->get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * bigO = stk::mesh::field_data(*omega, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      bigO[k] = scalarOmega;
    }
  }
}

//--------------------------------------------------------------------------
//-------- set_current_displacement ----------------------------------------
//--------------------------------------------------------------------------
void
Realm::set_current_displacement(
  stk::mesh::Part *targetPart)
{
  const int nDim = metaData_->spatial_dimension();
  const double currentTime = get_current_time();

  VectorFieldType *modelCoords = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates");
  VectorFieldType *displacement = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "mesh_displacement");
  ScalarFieldType *omega = metaData_->get_field<ScalarFieldType>(stk::topology::NODE_RANK, "omega");

  stk::mesh::Selector s_all_nodes = stk::mesh::Selector(*targetPart);

  stk::mesh::BucketVector const& node_buckets = bulkData_->get_buckets( stk::topology::NODE_RANK, s_all_nodes );

  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    const double * bigO = stk::mesh::field_data(*omega, b);
    const double * mCoords = stk::mesh::field_data(*modelCoords, b);
    double * dx = stk::mesh::field_data(*displacement, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // extract omega
      const double theO = bigO[k];

      const int offSet = k*nDim;
      // hacked for 2D
      const double cX = mCoords[offSet];
      const double cY = mCoords[offSet+1];

      dx[offSet] =  ((cos(theO*currentTime) - 1.0 )*cX
                     - sin(theO*currentTime)*cY);
      dx[offSet+1] =  (sin(theO*currentTime)*cX
                       + (cos(theO*currentTime)-1.0)*cY);
    }
  }
}

//--------------------------------------------------------------------------
//-------- set_current_coordinates -----------------------------------------
//--------------------------------------------------------------------------
void
Realm::set_current_coordinates(
  stk::mesh::Part *targetPart)
{
  const int nDim = metaData_->spatial_dimension();

  VectorFieldType *modelCoords = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates");
  VectorFieldType *currentCoords = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "current_coordinates");
  VectorFieldType *displacement = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "mesh_displacement");

  stk::mesh::Selector s_all_nodes = stk::mesh::Selector(*targetPart);

  stk::mesh::BucketVector const& node_buckets = bulkData_->get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    const double * mCoords = stk::mesh::field_data(*modelCoords, b);
    double * cCoords = stk::mesh::field_data(*currentCoords, b);
    const double * dx = stk::mesh::field_data(*displacement, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j )
        cCoords[offSet+j] = mCoords[offSet+j] + dx[offSet+j];
    }
  }
}

//--------------------------------------------------------------------------
//-------- set_mesh_velocity -----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::set_mesh_velocity(
  stk::mesh::Part *targetPart)
{
  const int nDim = metaData_->spatial_dimension();

  VectorFieldType *currentCoords = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "current_coordinates");
  VectorFieldType *meshVelocity = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "mesh_velocity");
  ScalarFieldType *omega = metaData_->get_field<ScalarFieldType>(stk::topology::NODE_RANK, "omega");

  stk::mesh::Selector s_all_nodes = stk::mesh::Selector(*targetPart);

  stk::mesh::BucketVector const& node_buckets = bulkData_->get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    const double * cCoords = stk::mesh::field_data(*currentCoords, b);
    const double * bigO = stk::mesh::field_data(*omega, b);
    double * vnp1 = stk::mesh::field_data(*meshVelocity, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int offSet = k*nDim;

      // hacked for 2D
      const double cX = cCoords[offSet];
      const double cY = cCoords[offSet+1];
      const double theO = bigO[k];
      vnp1[offSet] =  -theO*cY;
      vnp1[offSet+1] = theO*cX;
    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_geometry ------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::compute_geometry()
{
  // perform extrusion of mesh
  if ( hasContact_ )
    extrusionMeshDistanceAlgDriver_->execute();
  computeGeometryAlgDriver_->execute();

  // find total volume if the mesh moves at all
  if ( does_mesh_move() ) {
    double totalVolume = 0.0;
    double maxVolume = -1.0e16;
    double minVolume = 1.0e16;

    ScalarFieldType *dualVolume = metaData_->get_field<ScalarFieldType>(stk::topology::NODE_RANK, "dual_nodal_volume");

    stk::mesh::Selector s_local_nodes
      = metaData_->locally_owned_part() &stk::mesh::selectField(*dualVolume);

    stk::mesh::BucketVector const& node_buckets = bulkData_->get_buckets( stk::topology::NODE_RANK, s_local_nodes );
    for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
	  ib != node_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      const stk::mesh::Bucket::size_type length   = b.size();
      const double * dv = stk::mesh::field_data(*dualVolume, b);
      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        const double theVol = dv[k];
        totalVolume += theVol;
        maxVolume = std::max(theVol, maxVolume);
        minVolume = std::min(theVol, minVolume);
      }
    }

    // get min, max and sum over processes
    double g_totalVolume = 0.0, g_minVolume = 0.0, g_maxVolume = 0.0;
    stk::all_reduce_min(NaluEnv::self().parallel_comm(), &minVolume, &g_minVolume, 1);
    stk::all_reduce_max(NaluEnv::self().parallel_comm(), &maxVolume, &g_maxVolume, 1);
    stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &totalVolume, &g_totalVolume, 1);

    NaluEnv::self().naluOutputP0() << " Volume  " << g_totalVolume
		    << " min: " << g_minVolume
		    << " max: " << g_maxVolume << std::endl;
  }
}

//--------------------------------------------------------------------------
//-------- compute_vrtm ----------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::compute_vrtm()
{
  // compute velocity relative to mesh; must be tied to velocity update...
  if ( solutionOptions_->meshMotion_ || solutionOptions_->externalMeshDeformation_ ) {
    const int nDim = metaData_->spatial_dimension();

    VectorFieldType *velocity = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
    VectorFieldType *meshVelocity = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "mesh_velocity");
    VectorFieldType *velocityRTM = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm");

    stk::mesh::Selector s_all_nodes
       = (metaData_->locally_owned_part() | metaData_->globally_shared_part());

    stk::mesh::BucketVector const& node_buckets = bulkData_->get_buckets( stk::topology::NODE_RANK, s_all_nodes );
    for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
	  ib != node_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      const stk::mesh::Bucket::size_type length   = b.size();
      const double * uNp1 = stk::mesh::field_data(*velocity, b);
      const double * vNp1 = stk::mesh::field_data(*meshVelocity, b);
      double * vrtm = stk::mesh::field_data(*velocityRTM, b);

      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        const int offSet = k*nDim;
        for ( int j=0; j < nDim; ++j ) {
          vrtm[offSet+j] = uNp1[offSet+j] - vNp1[offSet+j];
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- init_current_coordinates -----------------------------------------
//--------------------------------------------------------------------------
void
Realm::init_current_coordinates()
{

  const int nDim = metaData_->spatial_dimension();

  VectorFieldType *modelCoords = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates");
  VectorFieldType *currentCoords = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "current_coordinates");
  VectorFieldType *displacement = metaData_->get_field<VectorFieldType>(stk::topology::NODE_RANK, "mesh_displacement");

  stk::mesh::Selector s_all_nodes
    = (metaData_->locally_owned_part() | metaData_->globally_shared_part());

  stk::mesh::BucketVector const& node_buckets = bulkData_->get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    const double * mCoords = stk::mesh::field_data(*modelCoords, b);
    double * cCoords = stk::mesh::field_data(*currentCoords, b);
    double * dx = stk::mesh::field_data(*displacement, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j ) {
        dx[offSet+j] = 0.0; //RESTART...
        cCoords[offSet+j] = mCoords[offSet+j];
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_l2_scaling ----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::compute_l2_scaling()
{
  // loop over all material propertys  and save off part vector
  stk::mesh::PartVector partVec;
  const std::vector<std::string> targetNames = materialPropertys_.targetNames_;
  for (size_t itarget=0; itarget < targetNames.size(); ++itarget) {
    // target need not be subsetted since nothing below will depend on topo
    stk::mesh::Part *targetPart = metaData_->get_part(targetNames[itarget]);
    partVec.push_back(targetPart);
  }

  size_t totalNodes = 0;

  // selector for all locally owned nodes
  stk::mesh::Selector s_locally_owned_union =
    metaData_->locally_owned_part()
    &stk::mesh::selectUnion(partVec);

  stk::mesh::BucketVector const& node_bucket = bulkData_->get_buckets( stk::topology::NODE_RANK, s_locally_owned_union );

  for ( stk::mesh::BucketVector::const_iterator ib = node_bucket.begin() ;
        ib != node_bucket.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    totalNodes += length;
  }

  // Parallel assembly of total nodes
  size_t g_totalNodes = 0;
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &totalNodes, &g_totalNodes, 1);

  l2Scaling_ = 1.0/std::sqrt(g_totalNodes);

}


//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_nodal_fields(
  stk::mesh::Part *part)
{
  // register high level common fields
  const int nDim = metaData_->spatial_dimension();

  // mesh motion/deformation is high level
  if ( solutionOptions_->meshMotion_ || solutionOptions_->externalMeshDeformation_) {
    VectorFieldType *displacement = &(metaData_->declare_field<VectorFieldType>(stk::topology::NODE_RANK, "mesh_displacement"));
    stk::mesh::put_field(*displacement, *part, nDim);
    VectorFieldType *currentCoords = &(metaData_->declare_field<VectorFieldType>(stk::topology::NODE_RANK, "current_coordinates"));
    stk::mesh::put_field(*currentCoords, *part, nDim);
    VectorFieldType *meshVelocity = &(metaData_->declare_field<VectorFieldType>(stk::topology::NODE_RANK, "mesh_velocity"));
    stk::mesh::put_field(*meshVelocity, *part, nDim);
    VectorFieldType *velocityRTM = &(metaData_->declare_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity_rtm"));
    stk::mesh::put_field(*velocityRTM, *part, nDim);
    // only internal mesh motion requires rotation rate
    if ( solutionOptions_->meshMotion_ ) {
      ScalarFieldType *omega = &(metaData_->declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "omega"));
      stk::mesh::put_field(*omega, *part);
    }
    // only external mesh deformation requires dvi/dxj (for GCL)
    if ( solutionOptions_->externalMeshDeformation_) {
      ScalarFieldType *divV = &(metaData_->declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "div_mesh_velocity"));
      stk::mesh::put_field(*divV, *part);
    }
  }

  if ( averagingInfo_->processAveraging_ )
    register_averaging_variables(part);
}

//--------------------------------------------------------------------------
//-------- register_averaging_variables ------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_averaging_variables(stk::mesh::Part *part)
{
  // register high level common fields
  const int nDim = metaData_->spatial_dimension();

  // always need Reynolds averaged density; add as restart variable
  std::string densReName = "density_ra";
  ScalarFieldType *densityRA =  &(metaData_->declare_field<ScalarFieldType>(stk::topology::NODE_RANK, densReName));
  stk::mesh::put_field(*densityRA, *part);
  augment_restart_variable_list(densReName);

  // favre
  for ( size_t i = 0; i < averagingInfo_->favreFieldNameVec_.size(); ++i ) {
    const std::string varName = averagingInfo_->favreFieldNameVec_[i];
    const std::string favreName = varName + "_fa";
    if ( varName == "velocity" ) {
      VectorFieldType *velocity = &(metaData_->declare_field<VectorFieldType>(stk::topology::NODE_RANK, favreName));
      stk::mesh::put_field(*velocity, *part, nDim);
    }
    else {
      ScalarFieldType *scalarQ = &(metaData_->declare_field<ScalarFieldType>(stk::topology::NODE_RANK, favreName));
      stk::mesh::put_field(*scalarQ, *part);
    }

    // add to restart
    augment_restart_variable_list(favreName);
  }

  // reynolds
  for ( size_t i = 0; i < averagingInfo_->reynoldsFieldNameVec_.size(); ++i ) {
    const std::string varName = averagingInfo_->reynoldsFieldNameVec_[i];
    const std::string reynoldsName = varName + "_ra";
    if ( varName == "velocity" ) {
      VectorFieldType *velocity = &(metaData_->declare_field<VectorFieldType>(stk::topology::NODE_RANK, reynoldsName));
      stk::mesh::put_field(*velocity, *part, nDim);
    }
    else {
      ScalarFieldType *scalarQ = &(metaData_->declare_field<ScalarFieldType>(stk::topology::NODE_RANK, reynoldsName));
      stk::mesh::put_field(*scalarQ, *part);
    }

    // add to restart
    augment_restart_variable_list(reynoldsName);
  }
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_interior_algorithm(
  stk::mesh::Part *part)
{
  //====================================================
  // Register interior algorithms
  //====================================================
  const AlgorithmType algType = INTERIOR;
  std::map<AlgorithmType, Algorithm *>::iterator it
    = computeGeometryAlgDriver_->algMap_.find(algType);
  if ( it == computeGeometryAlgDriver_->algMap_.end() ) {
    ComputeGeometryInteriorAlgorithm *theAlg
      = new ComputeGeometryInteriorAlgorithm(*this, part);
    computeGeometryAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

  // nodal post processing; only averaging thus far
  if ( averagingInfo_->processAveraging_ ) {
    if ( NULL == postConvergedAlgDriver_ )
      postConvergedAlgDriver_ = new AlgorithmDriver(*this);

    std::map<AlgorithmType, Algorithm *>::iterator it_pp
      = postConvergedAlgDriver_->algMap_.find(algType);
    if ( it_pp == postConvergedAlgDriver_->algMap_.end() ) {
      TurbulenceAveragingAlgorithm *theAlg
	= new TurbulenceAveragingAlgorithm(*this, part);
      postConvergedAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      it_pp->second->partVec_.push_back(part);
    }
  }
}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{

  //====================================================
  // Register face (boundary condition) data
  //====================================================

  // push back the part for book keeping and, later, skin mesh
  bcPartVec_.push_back(part);

  const int nDim = metaData_->spatial_dimension();

  // register fields
  MasterElement *meFC = get_surface_master_element(theTopo);
  const int numScsIp = meFC->numIntPoints_;

  GenericFieldType *exposedAreaVec_
    = &(metaData_->declare_field<GenericFieldType>(static_cast<stk::topology::rank_t>(metaData_->side_rank()), "exposed_area_vector"));
  stk::mesh::put_field(*exposedAreaVec_, *part, nDim*numScsIp );

  //====================================================
  // Register wall algorithms
  //====================================================
  const AlgorithmType algType = WALL;
  std::map<AlgorithmType, Algorithm *>::iterator it
    = computeGeometryAlgDriver_->algMap_.find(algType);
  if ( it == computeGeometryAlgDriver_->algMap_.end() ) {
    ComputeGeometryBoundaryAlgorithm *theAlg
      = new ComputeGeometryBoundaryAlgorithm(*this, part);
    computeGeometryAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{

  //====================================================
  // Register face (boundary condition) data
  //====================================================

  // push back the part for book keeping and, later, skin mesh
  bcPartVec_.push_back(part);

  const int nDim = metaData_->spatial_dimension();

  // register fields
  MasterElement *meFC = get_surface_master_element(theTopo);
  const int numScsIp = meFC->numIntPoints_;

  GenericFieldType *exposedAreaVec_
    = &(metaData_->declare_field<GenericFieldType>(static_cast<stk::topology::rank_t>(metaData_->side_rank()), "exposed_area_vector"));
  stk::mesh::put_field(*exposedAreaVec_, *part, nDim*numScsIp );

  //====================================================
  // Register wall algorithms
  //====================================================
  const AlgorithmType algType = INFLOW;
  std::map<AlgorithmType, Algorithm *>::iterator it
    = computeGeometryAlgDriver_->algMap_.find(algType);
  if ( it == computeGeometryAlgDriver_->algMap_.end() ) {
    ComputeGeometryBoundaryAlgorithm *theAlg
      = new ComputeGeometryBoundaryAlgorithm(*this, part);
    computeGeometryAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_open_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{

  //====================================================
  // Register face (boundary condition) data
  //====================================================

  // push back the part for book keeping and, later, skin mesh
  bcPartVec_.push_back(part);

  const int nDim = metaData_->spatial_dimension();

  // register fields
  MasterElement *meFC = get_surface_master_element(theTopo);
  const int numScsIp = meFC->numIntPoints_;

  GenericFieldType *exposedAreaVec_
    = &(metaData_->declare_field<GenericFieldType>(static_cast<stk::topology::rank_t>(metaData_->side_rank()), "exposed_area_vector"));
  stk::mesh::put_field(*exposedAreaVec_, *part, nDim*numScsIp );

  //====================================================
  // Register wall algorithms
  //====================================================
  const AlgorithmType algType = OPEN;
  std::map<AlgorithmType, Algorithm *>::iterator it
    = computeGeometryAlgDriver_->algMap_.find(algType);
  if ( it == computeGeometryAlgDriver_->algMap_.end() ) {
    ComputeGeometryBoundaryAlgorithm *theAlg
      = new ComputeGeometryBoundaryAlgorithm(*this, part);
    computeGeometryAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_contact_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_contact_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo,
  const ContactBoundaryConditionData &contactBCData)
{

  // push back the part for book keeping and, later, skin mesh
  bcPartVec_.push_back(part);

  const AlgorithmType algType = CONTACT;

  hasContact_ = true;

  if ( NULL == extrusionMeshDistanceAlgDriver_ )
    extrusionMeshDistanceAlgDriver_ = new AlgorithmDriver(*this);

  //====================================================
  // Register boundary condition data
  //====================================================

  const int nDim = metaData_->spatial_dimension();

  // register fields
  MasterElement *meFC = get_surface_master_element(theTopo);
  const int numScsIp = meFC->numIntPoints_;

  // exposed area vector
  GenericFieldType *exposedAreaVec_
    = &(metaData_->declare_field<GenericFieldType>(static_cast<stk::topology::rank_t>(metaData_->side_rank()), "exposed_area_vector"));
  stk::mesh::put_field(*exposedAreaVec_, *part, nDim*numScsIp );

  // register nodal field that will hold important information
  VectorFieldType *haloNormal =  &(metaData_->declare_field<VectorFieldType>(stk::topology::NODE_RANK, "halo_normal"));
  stk::mesh::put_field(*haloNormal, *part, nDim);
  VectorFieldType *haloDxj =  &(metaData_->declare_field<VectorFieldType>(stk::topology::NODE_RANK, "halo_dxj"));
  stk::mesh::put_field(*haloDxj, *part, nDim);
  ScalarFieldType *extDistance = &(metaData_->declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "extrusion_distance"));
  stk::mesh::put_field(*extDistance, *part);

  // correction for extrusion distance (for non-planar surfaces)
  ScalarFieldType *extDistanceCorrFac 
    = &(metaData_->declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "extrusion_distance_correct_fac"));
  stk::mesh::put_field(*extDistanceCorrFac, *part);
  ScalarFieldType *extDistanceCorrCount 
    = &(metaData_->declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "extrusion_distance_correct_count"));
  stk::mesh::put_field(*extDistanceCorrCount, *part);

  if ( realmUsesEdges_ ) {
    // need some extra nodal data for edge-based
    ScalarFieldType *haloMdot =  &(metaData_->declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "halo_mdot"));
    stk::mesh::put_field(*haloMdot, *part);

    VectorFieldType *haloAxj =  &(metaData_->declare_field<VectorFieldType>(stk::topology::NODE_RANK, "halo_axj"));
    stk::mesh::put_field(*haloAxj, *part, nDim);
  }
  else {
    throw std::runtime_error("Element-based scheme not supported with contact bc");
  }

  // extract data
  ContactUserData userData = contactBCData.userData_;

  // extract params useful for search
  const double maxSearchRadius = userData.maxSearchRadius_;
  const double minSearchRadius = userData.minSearchRadius_;
  const double extrusionDistance = userData.extrusionDistance_;
  const std::string searchMethodName = userData.searchMethodName_;
  const double expandBoxPercentage = userData.expandBoxPercentage_/100.0;
  const bool clipIsoParametricCoords = userData.clipIsoParametricCoords_;

  // what type of interpolation?
  const bool useHermiteInterp = userData.useHermiteInterpolation_;

  // some output for the user
  NaluEnv::self().naluOutputP0() << "extrusion alg active with distance " << extrusionDistance << std::endl;
  NaluEnv::self().naluOutputP0() << "min/max radius " << minSearchRadius << " " << maxSearchRadius << std::endl;
  NaluEnv::self().naluOutputP0() << "search block name: " << std::endl;
  for ( size_t k = 0; k < userData.searchBlock_.size(); ++k )
    NaluEnv::self().naluOutputP0() << userData.searchBlock_[k] << std::endl;

  // create manager
  if ( NULL == contactManager_ ) {
    contactManager_ = new ContactManager(*this);
  }

  // create contact info for this surface
  ContactInfo *contactInfo
    = new ContactInfo(*this,
                      part->name(),
                      maxSearchRadius,
                      minSearchRadius,
                      userData.searchBlock_,
                      expandBoxPercentage,
                      part,
                      searchMethodName,
                      clipIsoParametricCoords,
                      useHermiteInterp);


  contactManager_->contactInfoVec_.push_back(contactInfo);

  //====================================================
  // Register contact algorithms
  //====================================================

  // handle boundary data; for now this is constant extrusion distance (with non-planar correction)
  std::vector<double> userSpecBc(1);
  userSpecBc[0] = extrusionDistance;
  ConstantAuxFunction *theAuxFuncBc = new ConstantAuxFunction(0, 1, userSpecBc);

  // bc data alg
  AuxFunctionAlgorithm *auxAlgBc
    = new AuxFunctionAlgorithm(*this, part,
                               extDistance, theAuxFuncBc,
                               stk::topology::NODE_RANK);
  bcDataAlg_.push_back(auxAlgBc);

  // ExtrusionMeshDistanceBoundaryAlgorithm handles exposed area vector
  std::map<AlgorithmType, Algorithm *>::iterator it
    = extrusionMeshDistanceAlgDriver_->algMap_.find(algType);
  if ( it == extrusionMeshDistanceAlgDriver_->algMap_.end() ) {
    ExtrusionMeshDistanceBoundaryAlgorithm *theAlg
      = new ExtrusionMeshDistanceBoundaryAlgorithm(*this, part );
    extrusionMeshDistanceAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

  // extruded mesh (volume and area vector)
  if ( realmUsesEdges_ ) {
    std::map<AlgorithmType, Algorithm *>::iterator itext
      = computeGeometryAlgDriver_->algMap_.find(algType);
    if ( itext == computeGeometryAlgDriver_->algMap_.end() ) {
      ComputeGeometryExtrusionBoundaryAlgorithm *theAlg
        = new ComputeGeometryExtrusionBoundaryAlgorithm(*this, part);
      computeGeometryAlgDriver_->algMap_[algType] = theAlg;
    }
    else {
      itext->second->partVec_.push_back(part);
    }
  }

}

//--------------------------------------------------------------------------
//-------- register_symmetry_bc --------------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{

  //====================================================
  // Register face (boundary condition) data
  //====================================================

  // push back the part for book keeping and, later, skin mesh
  bcPartVec_.push_back(part);

  const int nDim = metaData_->spatial_dimension();

  // register fields
  MasterElement *meFC = get_surface_master_element(theTopo);
  const int numScsIp = meFC->numIntPoints_;

  GenericFieldType *exposedAreaVec_
    = &(metaData_->declare_field<GenericFieldType>(static_cast<stk::topology::rank_t>(metaData_->side_rank()), "exposed_area_vector"));
  stk::mesh::put_field(*exposedAreaVec_, *part, nDim*numScsIp );

  //====================================================
  // Register symmetry algorithms
  //====================================================
  const AlgorithmType algType = SYMMETRY;
  std::map<AlgorithmType, Algorithm *>::iterator it
    = computeGeometryAlgDriver_->algMap_.find(algType);
  if ( it == computeGeometryAlgDriver_->algMap_.end() ) {
    ComputeGeometryBoundaryAlgorithm *theAlg
      = new ComputeGeometryBoundaryAlgorithm(*this, part);
    computeGeometryAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_periodic_bc --------------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_periodic_bc(
  stk::mesh::Part *masterMeshPart,
  stk::mesh::Part *slaveMeshPart,
  const double &searchTolerance,
  const std::string &searchMethodName)
{

  // push back the part for book keeping and, later, skin mesh
  bcPartVec_.push_back(masterMeshPart);
  bcPartVec_.push_back(slaveMeshPart);

  if ( NULL == periodicManager_ ) {
    periodicManager_ = new PeriodicManager(*this);
    hasPeriodic_ = true;
  }

  // add the parts to the manager
  periodicManager_->add_periodic_pair(masterMeshPart, slaveMeshPart, searchTolerance, searchMethodName);
}

//--------------------------------------------------------------------------
//-------- setup_non_conformal_bc ------------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_non_conformal_bc(
  stk::mesh::Part *currentPart,
  stk::mesh::Part *opposingPart,
  const NonConformalBoundaryConditionData &nonConformalBCData)
{

  hasNonConformal_ = true;
  
  // extract data
  NonConformalUserData userData = nonConformalBCData.userData_;

  // extract params useful for search
  const std::string searchMethodName = userData.searchMethodName_;
  const double expandBoxPercentage = userData.expandBoxPercentage_/100.0;
  const bool clipIsoParametricCoords = userData.clipIsoParametricCoords_; 
  const double searchTolerance = userData.searchTolerance_;

  // deal with output
  const bool ncAlgDetailedOutput = solutionOptions_->ncAlgDetailedOutput_;

  // create manager
  if ( NULL == nonConformalManager_ ) {
    nonConformalManager_ = new NonConformalManager(*this, ncAlgDetailedOutput);
  }
   
  // create contact info for this surface
  NonConformalInfo *nonConformalInfo
    = new NonConformalInfo(*this,
                           currentPart,
                           opposingPart,
                           expandBoxPercentage,
                           searchMethodName,
                           clipIsoParametricCoords,
                           searchTolerance);
  
  nonConformalManager_->nonConformalInfoVec_.push_back(nonConformalInfo);
}

//--------------------------------------------------------------------------
//-------- register_non_conformal_bc ---------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_non_conformal_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{

  // push back the part for book keeping and, later, skin mesh
  bcPartVec_.push_back(part);

  const AlgorithmType algType = NON_CONFORMAL;

  const int nDim = metaData_->spatial_dimension();
  
  // register fields
  MasterElement *meFC = get_surface_master_element(theTopo);
  const int numScsIp = meFC->numIntPoints_;
  
  // exposed area vector
  GenericFieldType *exposedAreaVec_
    = &(metaData_->declare_field<GenericFieldType>(static_cast<stk::topology::rank_t>(metaData_->side_rank()), "exposed_area_vector"));
  stk::mesh::put_field(*exposedAreaVec_, *part, nDim*numScsIp );
   
  //====================================================
  // Register non-conformal algorithms
  //====================================================
  std::map<AlgorithmType, Algorithm *>::iterator it
    = computeGeometryAlgDriver_->algMap_.find(algType);
  if ( it == computeGeometryAlgDriver_->algMap_.end() ) {
    ComputeGeometryBoundaryAlgorithm *theAlg
      = new ComputeGeometryBoundaryAlgorithm(*this, part);
    computeGeometryAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- setup_overset_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_overset_bc(
  const OversetBoundaryConditionData &oversetBCData)
{
  // setting flag for linear system setup (may have been set via earlier "query")
  hasOverset_ = true;
  
  // create manager while providing overset data
  if ( NULL == oversetManager_ ) {
    oversetManager_ = new OversetManager(*this,oversetBCData.userData_);
  }   
}

//--------------------------------------------------------------------------
//-------- periodic_field_update -------------------------------------------
//--------------------------------------------------------------------------
void
Realm::periodic_field_update(
  stk::mesh::FieldBase *theField,
  const unsigned &sizeOfField,
  const bool &bypassFieldCheck) const
{
  const bool addSlaves = true;
  const bool setSlaves = true;
  periodicManager_->apply_constraints(theField, sizeOfField, bypassFieldCheck, addSlaves, setSlaves);
}

//--------------------------------------------------------------------------
//-------- periodic_delta_solution_update -------------------------------------------
//--------------------------------------------------------------------------
void
Realm::periodic_delta_solution_update(
  stk::mesh::FieldBase *theField,
  const unsigned &sizeOfField) const
{
  const bool bypassFieldCheck = true;
  const bool addSlaves = false;
  const bool setSlaves = true;
  periodicManager_->apply_constraints(theField, sizeOfField, bypassFieldCheck, addSlaves, setSlaves);
}

//--------------------------------------------------------------------------
//-------- periodic_max_field_update ---------------------------------------
//--------------------------------------------------------------------------
void
Realm::periodic_max_field_update(
  stk::mesh::FieldBase *theField,
  const unsigned &sizeOfField) const
{
  periodicManager_->apply_max_field(theField, sizeOfField);
}

//--------------------------------------------------------------------------
//-------- get_slave_part_vector -------------------------------------------
//--------------------------------------------------------------------------
const stk::mesh::PartVector &
Realm::get_slave_part_vector()
{
  if ( hasPeriodic_)
    return periodicManager_->get_slave_part_vector();
  else
    return emptyPartVector_;
}


//--------------------------------------------------------------------------
//-------- overset_field_update -------------------------------------------
//--------------------------------------------------------------------------
void
Realm::overset_orphan_node_field_update(
  stk::mesh::FieldBase *theField,
  const unsigned sizeRow,
  const unsigned sizeCol)
{
  oversetManager_->overset_orphan_node_field_update(theField, sizeRow, sizeCol);
}

//--------------------------------------------------------------------------
//-------- provide_output --------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::provide_output()
{
  stk::diag::TimeBlock mesh_output_timeblock(Simulation::outputTimer());

  if ( outputInfo_->hasOutputBlock_ ) {

    if (outputInfo_->outputFreq_ == 0)
      return;

    const double start_time = stk::cpu_time();

    // process output via io
    const double currentTime = get_current_time();
    const int timeStepCount = get_time_step_count();
    const bool isOutput = (timeStepCount % outputInfo_->outputFreq_) == 0;

    if ( isOutput ) {
      // when adaptivity has occurred, re-create the output mesh file
      if (outputInfo_->meshAdapted_)
        create_output_mesh();

      // not set up for globals
      ioBroker_->process_output_request(resultsFileIndex_, currentTime);
      equationSystems_.provide_output();
    }

    const double stop_time = stk::cpu_time();

    // increment time for output
    timerOutputFields_ += (stop_time - start_time);
  }
}

//--------------------------------------------------------------------------
//-------- provide_restart_output ------------------------------------------
//--------------------------------------------------------------------------
void
Realm::provide_restart_output()
{
  stk::diag::TimeBlock mesh_output_timeblock(Simulation::outputTimer());

  if ( outputInfo_->hasRestartBlock_ ) {

    if (outputInfo_->restartFreq_ == 0)
      return;

    const double start_time = stk::cpu_time();

    // process restart via io
    const double currentTime = get_current_time();
    const int timeStepCount = get_time_step_count();
    const bool isRestartOutputStep = (timeStepCount % outputInfo_->restartFreq_) == 0;

    if ( isRestartOutputStep ) {

      // handle fields
      ioBroker_->begin_output_step(restartFileIndex_, currentTime);
      ioBroker_->write_defined_output_fields(restartFileIndex_);

      // push global variables for time step
      const double timeStepNm1 = timeIntegrator_->get_time_step();
      globalParameters_.set_value("timeStepNm1", timeStepNm1);
      globalParameters_.set_value("timeStepCount", timeStepCount);
      if ( averagingInfo_->processAveraging_ )
        globalParameters_.set_value("currentTimeFilter", averagingInfo_->currentTimeFilter_ );

      stk::util::ParameterMapType::const_iterator i = globalParameters_.begin();
      stk::util::ParameterMapType::const_iterator iend = globalParameters_.end();
      for (; i != iend; ++i)
      {
        std::string parameterName = (*i).first;
        stk::util::Parameter parameter = (*i).second;
        if ( parameter.toRestartFile ) {
          ioBroker_->write_global(restartFileIndex_, parameterName,  parameter.value, parameter.type);
        }
      }

      ioBroker_->end_output_step(restartFileIndex_);
    }

    const double stop_time = stk::cpu_time();

    // increment time for output
    timerOutputFields_ += (stop_time - start_time);
  }

}

//--------------------------------------------------------------------------
//-------- swap_states -----------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::swap_states()
{
  bulkData_->update_field_data_states();
}

//--------------------------------------------------------------------------
//-------- predict_state ---------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::predict_state()
{
  equationSystems_.predict_state();
}

//--------------------------------------------------------------------------
//-------- populate_initial_condition --------------------------------------
//--------------------------------------------------------------------------
void
Realm::populate_initial_condition()
{
  for ( size_t k = 0; k < initCondAlg_.size(); ++k ) {
    initCondAlg_[k]->execute();
  }
}

//--------------------------------------------------------------------------
//-------- boundary_data_to_state_data -------------------------------------
//--------------------------------------------------------------------------
void
Realm::boundary_data_to_state_data()
{
  equationSystems_.boundary_data_to_state_data();
}

//--------------------------------------------------------------------------
//-------- populate_restart ------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::populate_restart(
  double &timeStepNm1, int &timeStepCount)
{
  double foundRestartTime = get_current_time();
  if ( restarted_simulation() ) {
    // allow restart to skip missed required fields
    const double restartTime = outputInfo_->restartTime_;
    std::vector<stk::io::MeshField> missingFields;
    foundRestartTime = ioBroker_->read_defined_input_fields(restartTime, &missingFields);
    if ( missingFields.size() > 0 ){
      for ( size_t k = 0; k < missingFields.size(); ++k) {
        NaluEnv::self().naluOutputP0() << "WARNING: Restart value for Field "
                                       << missingFields[k].field()->name()
                                       << " is missing; may default to IC specification" << std::endl;
      }
      if ( !supportInconsistentRestart_ ) {
        NaluEnv::self().naluOutputP0() << "The user may desire to set the support_inconsistent_multi_state_restart Realm line command" << std::endl;
        NaluEnv::self().naluOutputP0() << "This is applicable for a BDF2 restart run from a previously run Backward Euler simulation" << std::endl;
      }
    }
    
    // extract time parameters; okay if they are missing; no need to let the user know
    const bool abortIfNotFound = false;
    ioBroker_->get_global("timeStepNm1", timeStepNm1, abortIfNotFound);
    ioBroker_->get_global("timeStepCount", timeStepCount, abortIfNotFound);
    if ( averagingInfo_->processAveraging_)
      ioBroker_->get_global("currentTimeFilter", averagingInfo_->currentTimeFilter_, abortIfNotFound);
  }
  return foundRestartTime;
}

//--------------------------------------------------------------------------
//-------- populate_variables_from_input -----------------------------------
//--------------------------------------------------------------------------
void
Realm::populate_variables_from_input()
{
  // no reading fields from mesh if this is a restart
  if ( !restarted_simulation() && solutionOptions_->inputVarFromFileMap_.size() > 0 ) {
    const double timeToRead = 1.0e8;
    ioBroker_->read_defined_input_fields(timeToRead);
  }
}

//--------------------------------------------------------------------------
//-------- populate_derived_quantities -------------------------------------
//--------------------------------------------------------------------------
void
Realm::populate_derived_quantities()
{
  equationSystems_.populate_derived_quantities();
}

//--------------------------------------------------------------------------
//-------- initial_work -----------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::initial_work()
{
  equationSystems_.initial_work();
}

//--------------------------------------------------------------------------
//-------- set_global_id ---------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::set_global_id()
{
  const stk::mesh::Selector s_universal = metaData_->universal_part();
  stk::mesh::BucketVector const& buckets = bulkData_->get_buckets( stk::topology::NODE_RANK, s_universal );

  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin();
        ib != buckets.end(); ++ib ) {
    const stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();
    stk::mesh::EntityId *naluGlobalIds = stk::mesh::field_data(*naluGlobalId_, b);

    for ( stk::mesh::Bucket::size_type k = 0; k < length; ++k ) {
   	  naluGlobalIds[k] = bulkData_->identifier(b[k]);
    }
  }
}

//--------------------------------------------------------------------------
//-------- populate_boundary_data ------------------------------------------
//--------------------------------------------------------------------------
void
Realm::populate_boundary_data()
{
  // realm first
  for ( size_t k = 0; k < bcDataAlg_.size(); ++k ) {
    bcDataAlg_[k]->execute();
  }
  equationSystems_.populate_boundary_data();
}

//--------------------------------------------------------------------------
//-------- output_banner ---------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::output_banner()
{
  // in future, probably want to call through to eq sys..
  NaluEnv::self().naluOutputP0() << " Max Courant: " << maxCourant_ << " Max Reynolds: " << maxReynolds_ << std::endl;
}

//--------------------------------------------------------------------------
//-------- check_job -------------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::check_job(bool get_node_count)
{

  NaluEnv::self().naluOutputP0() << std::endl;
  NaluEnv::self().naluOutputP0() << "Realm memory Review:       " << name_ << std::endl;
  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;

  // set number of nodes, check job run size
  if (get_node_count)
  {
    size_t localNodeCount = ioBroker_->get_input_io_region()->get_property("node_count").get_int();
    stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &localNodeCount, &nodeCount_, 1);
    NaluEnv::self().naluOutputP0() << "Node count from meta data = " << nodeCount_ << std::endl;
  }

  /// estimate memory based on N*bandwidth, N = nodeCount*nDOF,
  ///   bandwidth = NCon(=27 for Hex mesh)*nDOF - we are very conservative here
  const unsigned HexBWFactor = 27;
  const unsigned MatrixStorageFactor = 3;  // for CRS storage, need one A_IJ, and one I and one J, approx
  SizeType memoryEstimate = 0;
  double procGBScale = double(NaluEnv::self().parallel_size())*(1024.*1024.*1024.);
  for (unsigned ieq=0; ieq < equationSystems_.size(); ++ieq)
    {
      if (!equationSystems_[ieq]->linsys_)
        continue;
      SizeType numDof = equationSystems_[ieq]->linsys_->numDof();
      SizeType N = nodeCount_ * numDof;
      SizeType bandwidth = HexBWFactor * numDof;
      memoryEstimate += MatrixStorageFactor * N*bandwidth * sizeof(double);
    }
  NaluEnv::self().naluOutputP0() << "Total memory estimate for Matrix solve (per core)= "
                  << double(memoryEstimate)/procGBScale << " GB." << std::endl;

  SizeType memoryEstimateFields = 0;
  if (metaData_->is_commit())
    {
      const stk::mesh::FieldVector & fields =  metaData_->get_fields();
      unsigned nfields = fields.size();
      for (unsigned ifld = 0; ifld < nfields; ++ifld)
        {
          stk::mesh::FieldBase *field = fields[ifld];
          // FIXME for element, edge, face-based fields
          unsigned fsz = field->max_size(stk::topology::NODE_RANK);
          memoryEstimateFields += nodeCount_ * fsz * sizeof(double);
        }
      NaluEnv::self().naluOutputP0() << "Total memory estimate for Fields (per core)= "
                      << double(memoryEstimateFields)/procGBScale << " GB." << std::endl;
      memoryEstimate += memoryEstimateFields;
    }

  NaluEnv::self().naluOutputP0() << "Total memory estimate (per core) = "
                  << double(memoryEstimate)/procGBScale << " GB." << std::endl;

  if (metaData_->is_commit() && estimateMemoryOnly_)
    {
      throw std::runtime_error("Job requested memory estimate only, shutting down");
    }

  // here's where we can check for estimated memory > given available memory
  if (availableMemoryPerCoreGB_ != 0
      && double(memoryEstimate)/procGBScale > availableMemoryPerCoreGB_)
    {
      NaluEnv::self().naluOutputP0() << "ERROR: property available_memory_per_core_GB is set (= " << availableMemoryPerCoreGB_
                      << ") and estimated memory (= " << double(memoryEstimate)/procGBScale
                      << ") is greater,\n job too large to run, \naborting..." << std::endl;
      throw std::runtime_error("Job shutting down");
    }
}

//--------------------------------------------------------------------------
//-------- dump_simulation_time --------------------------------------------
//--------------------------------------------------------------------------
void
Realm::dump_simulation_time()
{

  NaluEnv::self().naluOutputP0() << "Begin Timer Overview for Realm: " << name_ << std::endl;
  NaluEnv::self().naluOutputP0() << "-------------------------------- " << std::endl;

  // equation system time
  equationSystems_.dump_eq_time();

  const int nprocs = NaluEnv::self().parallel_size();

  // common
  const unsigned ntimers = 4;
  double total_time[ntimers] = {timerReadMesh_,timerOutputFields_, timerInitializeEqs_, timerPropertyEval_};
  double g_min_time[ntimers] = {}, g_max_time[ntimers] = {}, g_total_time[ntimers] = {};

  // get min, max and sum over processes
  stk::all_reduce_min(NaluEnv::self().parallel_comm(), &total_time[0], &g_min_time[0], ntimers);
  stk::all_reduce_max(NaluEnv::self().parallel_comm(), &total_time[0], &g_max_time[0], ntimers);
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &total_time[0], &g_total_time[0], ntimers);

  NaluEnv::self().naluOutputP0() << "Timing for IO: " << std::endl;
  NaluEnv::self().naluOutputP0() << "     io read mesh --  " << " \tavg: " << g_total_time[0]/double(nprocs)
                  << " \tmin: " << g_min_time[0] << " \tmax: " << g_max_time[0] << std::endl;
  NaluEnv::self().naluOutputP0() << " io output fields --  " << " \tavg: " << g_total_time[1]/double(nprocs)
                  << " \tmin: " << g_min_time[1] << " \tmax: " << g_max_time[1] << std::endl;

  NaluEnv::self().naluOutputP0() << "Timing for connectivity/finalize lysys: " << std::endl;
  NaluEnv::self().naluOutputP0() << "         eqs init --  " << " \tavg: " << g_total_time[2]/double(nprocs)
                  << " \tmin: " << g_min_time[2] << " \tmax: " << g_max_time[2] << std::endl;

  NaluEnv::self().naluOutputP0() << "Timing for property evaluation:         " << std::endl;
  NaluEnv::self().naluOutputP0() << "         props    --  " << " \tavg: " << g_total_time[3]/double(nprocs)
                  << " \tmin: " << g_min_time[3] << " \tmax: " << g_max_time[3] << std::endl;

  if (solutionOptions_->useAdapter_ && solutionOptions_->maxRefinementLevel_) {
    double g_total_adapt = 0.0, g_min_adapt = 0.0, g_max_adapt = 0.0;
    stk::all_reduce_min(NaluEnv::self().parallel_comm(), &timerAdapt_, &g_min_adapt, 1);
    stk::all_reduce_max(NaluEnv::self().parallel_comm(), &timerAdapt_, &g_max_adapt, 1);
    stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &timerAdapt_, &g_total_adapt, 1);

    NaluEnv::self().naluOutputP0() << "Timing for adaptivity:         " << std::endl;
    NaluEnv::self().naluOutputP0() << "         adapt    --  " << " \tavg: " << g_total_adapt/double(nprocs)
                    << " \tmin: " << g_min_adapt << " \tmax: " << g_max_adapt << std::endl;
  }

  // now edge creation; if applicable
  if ( realmUsesEdges_ ) {
    double g_total_edge = 0.0, g_min_edge = 0.0, g_max_edge = 0.0;
    stk::all_reduce_min(NaluEnv::self().parallel_comm(), &timerCreateEdges_, &g_min_edge, 1);
    stk::all_reduce_max(NaluEnv::self().parallel_comm(), &timerCreateEdges_, &g_max_edge, 1);
    stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &timerCreateEdges_, &g_total_edge, 1);

    NaluEnv::self().naluOutputP0() << "Timing for Edge: " << std::endl;

    NaluEnv::self().naluOutputP0() << "    edge creation --  " << " \tavg: " << g_total_edge/double(nprocs)
                    << " \tmin: " << g_min_edge << " \tmax: " << g_max_edge << std::endl;
  }

  // periodic
  if ( hasPeriodic_ ){
    double periodicSearchTime = periodicManager_->get_search_time();
    double g_minPeriodicSearchTime = 0.0, g_maxPeriodicSearchTime = 0.0, g_periodicSearchTime = 0.0;
    stk::all_reduce_min(NaluEnv::self().parallel_comm(), &periodicSearchTime, &g_minPeriodicSearchTime, 1);
    stk::all_reduce_max(NaluEnv::self().parallel_comm(), &periodicSearchTime, &g_maxPeriodicSearchTime, 1);
    stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &periodicSearchTime, &g_periodicSearchTime, 1);

    NaluEnv::self().naluOutputP0() << "Timing for Periodic: " << std::endl;
    NaluEnv::self().naluOutputP0() << "    search --         " << " \tavg: " << g_periodicSearchTime/double(nprocs)
                     << " \tmin: " << g_minPeriodicSearchTime << " \tmax: " << g_maxPeriodicSearchTime << std::endl;
  }

  // contact
  if ( hasContact_ ) {
    double g_totalContact = 0.0, g_minContact= 0.0, g_maxContact = 0.0;
    stk::all_reduce_min(NaluEnv::self().parallel_comm(), &timerContact_, &g_minContact, 1);
    stk::all_reduce_max(NaluEnv::self().parallel_comm(), &timerContact_, &g_maxContact, 1);
    stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &timerContact_, &g_totalContact, 1);

    NaluEnv::self().naluOutputP0() << "Timing for Contact: " << std::endl;

    NaluEnv::self().naluOutputP0() << "    contact bc    --  " << " \tavg: " << g_totalContact/double(nprocs)
                    << " \tmin: " << g_minContact << " \tmax: " << g_maxContact << std::endl;
  }

  // transfer
  if ( hasTransfer_ ) {
    double g_totalXfer = 0.0, g_minXfer= 0.0, g_maxXfer = 0.0;
    stk::all_reduce_min(NaluEnv::self().parallel_comm(), &timerTransferSearch_, &g_minXfer, 1);
    stk::all_reduce_max(NaluEnv::self().parallel_comm(), &timerTransferSearch_, &g_maxXfer, 1);
    stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &timerTransferSearch_, &g_totalXfer, 1);

    NaluEnv::self().naluOutputP0() << "Timing for Tranfer (fromRealm):    " << std::endl;

    NaluEnv::self().naluOutputP0() << "    search --         " << " \tavg: " << g_totalXfer/double(nprocs)
                    << " \tmin: " << g_minXfer << " \tmax: " << g_maxXfer << std::endl;
  }

}

//--------------------------------------------------------------------------
//-------- provide_mean_norm -----------------------------------------------
//--------------------------------------------------------------------------
double
Realm::provide_mean_norm()
{
  return equationSystems_.provide_mean_system_norm();
}

//--------------------------------------------------------------------------
//-------- get_volume_master_element ---------------------------------------
//--------------------------------------------------------------------------
MasterElement *
Realm::get_volume_master_element(
  const stk::topology & theTopo)
{

  MasterElement *theElem = NULL;

  std::map<stk::topology, MasterElement *>::iterator it =
    volumeMeMap_.find(theTopo);
  if ( it == volumeMeMap_.end() ) {
    // not found; will need to create it and add it
    switch ( theTopo.value() ) {

      case stk::topology::HEX_8:
        theElem = new HexSCV();
        break;

      case stk::topology::HEX_27:
        theElem = new Hex27SCV();
        break;

      case stk::topology::TET_4:
        theElem = new TetSCV();
        break;

      case stk::topology::PYRAMID_5:
        theElem = new PyrSCV();
        break;

      case stk::topology::WEDGE_6:
        theElem = new WedSCV();
        break;

      case stk::topology::QUAD_4_2D:
        theElem = new Quad2DSCV();
        break;

      case stk::topology::QUAD_9_2D:
        theElem = new Quad92DSCV();
        break;

      case stk::topology::TRI_3_2D:
        theElem = new Tri2DSCV();
        break;

      default:
        NaluEnv::self().naluOutputP0() << "sorry, we only support hex8, tet4, wed6, pyr5, quad4, and tri3 volume elements" << std::endl;
        break;
    }

    volumeMeMap_[theTopo] = theElem;

  }
  else {
    // found it
    theElem = it->second;
  }

  return theElem;
}

//--------------------------------------------------------------------------
//-------- get_surface_master_element ---------------------------------------
//--------------------------------------------------------------------------
MasterElement *
Realm::get_surface_master_element(
  const stk::topology & theTopo)
{

  MasterElement *theElem = NULL;

  std::map<stk::topology, MasterElement *>::iterator it =
    surfaceMeMap_.find(theTopo);
  if ( it == surfaceMeMap_.end() ) {
    // not found; will need to create it and add it
    switch ( theTopo.value() ) {

      case stk::topology::HEX_8:
        theElem = new HexSCS();
        break;

      case stk::topology::HEX_27:
        theElem = new Hex27SCS();
        break;

      case stk::topology::TET_4:
        theElem = new TetSCS();
        break;

      case stk::topology::PYRAMID_5:
        theElem = new PyrSCS();
        break;

      case stk::topology::WEDGE_6:
        theElem = new WedSCS();
        break;

      case stk::topology::QUAD_4:
        theElem =  new Quad3DSCS();
        break;

      case stk::topology::QUAD_9:
        theElem =  new Quad93DSCS();
        break;

      case stk::topology::TRI_3:
        theElem = new Tri3DSCS();
        break;

      case stk::topology::QUAD_4_2D:
        theElem =  new Quad2DSCS();
        break;

      case stk::topology::QUAD_9_2D:
        theElem =  new Quad92DSCS();
        break;

      case stk::topology::TRI_3_2D:
        theElem = new Tri2DSCS();
        break;

      case stk::topology::LINE_2:
        theElem = new Edge2DSCS();
        break;

      case stk::topology::LINE_3:
        theElem = new Edge32DSCS();
        break;

      default:
        NaluEnv::self().naluOutputP0() << "sorry, we only support hex8, tet4, pyr5, wed6, quad2d, quad3d, tri2d, tri3d and edge2d surface elements" << std::endl;
        NaluEnv::self().naluOutputP0() << "you're type is " << theTopo.value() << std::endl;
        break;

    }

    surfaceMeMap_[theTopo] = theElem;

  }

  else {
    // found it!
    theElem = it->second;
  }

  return theElem;

}

//--------------------------------------------------------------------------
//-------- get_hybrid_factor -----------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_hybrid_factor(
      const std::string dofName )
{
  double factor = solutionOptions_->hybridDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->hybridMap_.find(dofName);
  if (iter != solutionOptions_->hybridMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_alpha_factor ------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_alpha_factor(
      const std::string dofName )
{
  double factor = solutionOptions_->alphaDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->alphaMap_.find(dofName);
  if (iter != solutionOptions_->alphaMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_alpha_upw_factor --------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_alpha_upw_factor(
      const std::string dofName )
{
  double factor = solutionOptions_->alphaUpwDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->alphaUpwMap_.find(dofName);
  if (iter != solutionOptions_->alphaUpwMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_upw_factor ------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_upw_factor(
      const std::string dofName )
{
  double factor = solutionOptions_->upwDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->upwMap_.find(dofName);
  if (iter != solutionOptions_->upwMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- primitive_uses_limiter ------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::primitive_uses_limiter(
  const std::string dofName )
{
  bool usesIt = false;
  std::map<std::string, bool>::const_iterator iter
    = solutionOptions_->limiterMap_.find(dofName);
  if (iter != solutionOptions_->limiterMap_.end()) {
    usesIt = (*iter).second;
  }
  return usesIt;
}

//--------------------------------------------------------------------------
//-------- get_lam_schmidt -------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_lam_schmidt(
  const std::string dofName )
{
  double factor = solutionOptions_->lamScDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->lamScMap_.find(dofName);
  if (iter != solutionOptions_->lamScMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_lam_prandtl -------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_lam_prandtl(
  const std::string dofName, bool &prProvided )
{
  double factor = 1.0;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->lamPrMap_.find(dofName);
  if (iter != solutionOptions_->lamPrMap_.end()) {
    factor = (*iter).second;
    prProvided = true;
  }
  else {
    prProvided = false;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_turb_schmidt ------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_turb_schmidt(
  const std::string dofName )
{
  double factor = solutionOptions_->turbScDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->turbScMap_.find(dofName);
  if (iter != solutionOptions_->turbScMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_turb_prandtl ------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_turb_prandtl(
  const std::string dofName )
{
  double factor = solutionOptions_->turbPrDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->turbPrMap_.find(dofName);
  if (iter != solutionOptions_->turbPrMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_noc_usage ---------------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_noc_usage(
  const std::string dofName )
{
  bool factor = solutionOptions_->nocDefault_;
  std::map<std::string, bool>::const_iterator iter
    = solutionOptions_->nocMap_.find(dofName);
  if (iter != solutionOptions_->nocMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_divU --------------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_divU()
{
  return solutionOptions_->includeDivU_;
}

//--------------------------------------------------------------------------
//-------- get_mdot_interp -------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_mdot_interp()
{
  return solutionOptions_->mdotInterpRhoUTogether_ ? 1.0 : 0.0;
}

//--------------------------------------------------------------------------
//-------- get_cvfem_shifted_mdot ------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_cvfem_shifted_mdot()
{
  return solutionOptions_->cvfemShiftMdot_;
}

//--------------------------------------------------------------------------
//-------- get_cvfem_shifted_poisson ---------------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_cvfem_shifted_poisson()
{
  return solutionOptions_->cvfemShiftPoisson_;
}

//--------------------------------------------------------------------------
//-------- get_cvfem_reduced_sens_poisson ---------------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_cvfem_reduced_sens_poisson()
{
  return solutionOptions_->cvfemReducedSensPoisson_;
}

//--------------------------------------------------------------------------
//-------- has_nc_gauss_labatto_quadrature ---------------------------------
//--------------------------------------------------------------------------
bool
Realm::has_nc_gauss_labatto_quadrature()
{
  return solutionOptions_->ncAlgGaussLabatto_;
}

//--------------------------------------------------------------------------
//-------- get_nc_alg_type -------------------------------------------------
//--------------------------------------------------------------------------
NonConformalAlgType
Realm::get_nc_alg_type()
{
  return solutionOptions_->ncAlgType_;
}

//--------------------------------------------------------------------------
//-------- get_nc_alg_upwind_advection -------------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_nc_alg_upwind_advection()
{
  return solutionOptions_->ncAlgUpwindAdvection_;
}

//--------------------------------------------------------------------------
//-------- get_nc_alg_include_pstab ----------------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_nc_alg_include_pstab()
{
  return solutionOptions_->ncAlgIncludePstab_;
}

//--------------------------------------------------------------------------
//-------- get_material_prop_eval ------------------------------------------
//--------------------------------------------------------------------------
PropertyEvaluator *
Realm::get_material_prop_eval(
  const PropertyIdentifier thePropID)
{
  PropertyEvaluator *thePropEval = NULL;
  std::map<PropertyIdentifier, PropertyEvaluator*>::const_iterator iter
    = materialPropertys_.propertyEvalMap_.find(thePropID);
  if (iter != materialPropertys_.propertyEvalMap_.end()) {
    thePropEval = (*iter).second;
  }
  return thePropEval;
}

//--------------------------------------------------------------------------
//-------- is_turbulent ----------------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::is_turbulent()
{
  return solutionOptions_->isTurbulent_;
}

//--------------------------------------------------------------------------
//-------- is_turbulent ----------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::is_turbulent( bool isIt )
{
  isTurbulent_ = isIt;
  solutionOptions_->isTurbulent_ = isIt;
}

//--------------------------------------------------------------------------
//-------- needs_enthalpy --------------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::needs_enthalpy()
{
  return needsEnthalpy_;
}

//--------------------------------------------------------------------------
//-------- needs_enthalpy --------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::needs_enthalpy( bool needsEnthalpy )
{
  needsEnthalpy_ = needsEnthalpy;
}

//--------------------------------------------------------------------------
//-------- number_of_states ------------------------------------------------
//--------------------------------------------------------------------------
int
Realm::number_of_states()
{
  const int numStates = (timeIntegrator_->secondOrderTimeAccurate_) ? 3 : 2;
  return numStates;
}

//--------------------------------------------------------------------------
//-------- augment_transfer_list -------------------------------------------
//--------------------------------------------------------------------------
void
Realm::augment_transfer_vector(Transfer *transfer)
{
  transferVec_.push_back(transfer);
  hasTransfer_ = true;
}

//--------------------------------------------------------------------------
//-------- process_transfer ------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::process_transfer()
{
  if ( !hasTransfer_ )
    return;

  std::vector<Transfer *>::iterator ii;
  for( ii=transferVec_.begin(); ii!=transferVec_.end(); ++ii )
    (*ii)->execute();
}

//--------------------------------------------------------------------------
//-------- post_converged_work ---------------------------------------------
//--------------------------------------------------------------------------
void
Realm::post_converged_work()
{
  if ( NULL != postConvergedAlgDriver_ )
    postConvergedAlgDriver_->execute();

  for ( size_t k = 0; k < postConvergedAlg_.size(); ++k)
    postConvergedAlg_[k]->execute();

  if ( NULL != solutionNormPostProcessing_ )
    solutionNormPostProcessing_->execute();

  equationSystems_.post_converged_work();
}

//--------------------------------------------------------------------------
//-------- get_current_time() ----------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_current_time()
{
  return timeIntegrator_->get_current_time();
}

//--------------------------------------------------------------------------
//-------- get_time_step() ----------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_time_step()
{
  return timeIntegrator_->get_time_step();
}

//--------------------------------------------------------------------------
//-------- get_gamma1() ----------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_gamma1()
{
  return timeIntegrator_->get_gamma1();
}

//--------------------------------------------------------------------------
//-------- get_gamma2() ----------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_gamma2()
{
  return timeIntegrator_->get_gamma2();
}

//--------------------------------------------------------------------------
//-------- get_gamma3() ----------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_gamma3()
{
  return timeIntegrator_->get_gamma3();
}

//--------------------------------------------------------------------------
//-------- get_time_step_count() ----------------------------------------------
//--------------------------------------------------------------------------
int
Realm::get_time_step_count() const
{
  return timeIntegrator_->get_time_step_count();
}

//--------------------------------------------------------------------------
//-------- restarted_simulation() ------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::restarted_simulation()
{
  return outputInfo_->activateRestart_ ;
}

//--------------------------------------------------------------------------
//-------- support_inconsistent_restart() ----------------------------------
//--------------------------------------------------------------------------
bool
Realm::support_inconsistent_restart()
{
  return supportInconsistentRestart_ ;
}

//--------------------------------------------------------------------------
//-------- get_stefan_boltzmann() ------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_stefan_boltzmann()
{
  return solutionOptions_->stefanBoltzmann_;
}

//--------------------------------------------------------------------------
//-------- get_turb_model_constant() ------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_turb_model_constant(
   const TurbulenceModelConstant turbModelEnum)
{
  std::map<TurbulenceModelConstant, double>::iterator it
    = solutionOptions_->turbModelConstantMap_.find(turbModelEnum);
  if ( it != solutionOptions_->turbModelConstantMap_.end() ) {
    return it->second;
  }
  else {
    throw std::runtime_error("unknown (not found) turbulence model constant");
  }
}

//--------------------------------------------------------------------------
//-------- process_adaptivity() ----------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::process_adaptivity()
{
  const int timeStepCount = get_time_step_count();
  const bool processAdaptivity = (timeStepCount % solutionOptions_->adaptivityFrequency_) == 0 ? true : false;
  return processAdaptivity;
}

//--------------------------------------------------------------------------
//-------- get_buckets() ----------------------------------------------
//--------------------------------------------------------------------------
stk::mesh::BucketVector const& Realm::get_buckets( stk::mesh::EntityRank rank,
                                                   const stk::mesh::Selector & selector ,
                                                   bool get_all) const
{
  if (metaData_->spatial_dimension() == 3 && rank == stk::topology::EDGE_RANK)
    return bulkData_->get_buckets(rank, selector);

  if (!get_all && solutionOptions_->useAdapter_ && solutionOptions_->maxRefinementLevel_ > 0)
    {
      stk::mesh::Selector new_selector = selector;
      if (rank != stk::topology::NODE_RANK)
        {
          // adapterSelector_ avoids parent elements
          new_selector = selector & adapterSelector_[rank];
        }
      return bulkData_->get_buckets(rank, new_selector);
    }
  else
    {
      return bulkData_->get_buckets(rank, selector);
    }
}

//--------------------------------------------------------------------------
//-------- bulk_data() -----------------------------------------------------
//--------------------------------------------------------------------------
stk::mesh::BulkData &
Realm::bulk_data()
{
  return *bulkData_;
}

//--------------------------------------------------------------------------
//-------- meta_data() -----------------------------------------------------
//--------------------------------------------------------------------------
stk::mesh::MetaData &
Realm::meta_data()
{
  return *metaData_;
}

//--------------------------------------------------------------------------
//-------- get_activate_aura() -----------------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_activate_aura()
{
  return activateAura_;
}

//--------------------------------------------------------------------------
//-------- get_inactive_selector() -----------------------------------------
//--------------------------------------------------------------------------
stk::mesh::Selector
Realm::get_inactive_selector()
{
  // provide inactive part that excludes background surface
  return (hasOverset_) ? (stk::mesh::Selector(*oversetManager_->inActivePart_) 
                          &!(stk::mesh::selectUnion(oversetManager_->orphanPointSurfaceVecBackground_)))
    : stk::mesh::Selector();
}

//--------------------------------------------------------------------------
//-------- push_equation_to_systems() --------------------------------------
//--------------------------------------------------------------------------
void
Realm::push_equation_to_systems(
  EquationSystem *eqSystem)
{
  equationSystems_.equationSystemVector_.push_back(eqSystem);
}

} // namespace nalu
} // namespace Sierra
