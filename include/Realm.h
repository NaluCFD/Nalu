/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Realm_h
#define Realm_h

#include <Enums.h>
#include <FieldTypeDef.h>

// yaml for parsing..
#include <yaml-cpp/yaml.h>

#include <BoundaryConditions.h>
#include <InitialConditions.h>
#include <MaterialPropertys.h>
#include <EquationSystems.h>
#include <Teuchos_RCP.hpp>
#include <overset/OversetManager.h>

#include <stk_util/util/ParameterList.hpp>

// standard c++
#include <map>
#include <string>
#include <vector>
#include <stdint.h>

namespace stk {
namespace mesh {
class Part;
}
namespace io {
  class StkMeshIoBroker;
}
}

namespace YAML {
class Node;
}

namespace sierra{
namespace nalu{

class Algorithm;
class AlgorithmDriver;
class AuxFunctionAlgorithm;
class ComputeGeometryAlgorithmDriver;
// class OversetManager;
class NonConformalManager;
class EquationSystems;
class OutputInfo;
class PostProcessingInfo;
class PeriodicManager;
class Realms;
class Simulation;
class SolutionOptions;
class TimeIntegrator;
class MasterElement;
class PropertyEvaluator;
class HDF5FilePtr;
class Transfer;

class SolutionNormPostProcessing;
class TurbulenceAveragingPostProcessing;
class DataProbePostProcessing;
class Actuator;

/** Representation of a computational domain and physics equations solved on
 * this domain.
 *
 */
class Realm {
 public:

  Realm(Realms&, const YAML::Node & node);
  virtual ~Realm();
  
  typedef size_t SizeType;

  virtual void load(const YAML::Node & node);
  void look_ahead_and_creation(const YAML::Node & node);

  virtual void breadboard();

  virtual void initialize();
 
  Simulation *root() const;
  Simulation *root();
  Realms *parent() const;
  Realms *parent();

  bool debug() const;
  bool get_activate_memory_diagnostic();
  void provide_memory_summary();
  std::string convert_bytes(double bytes);

  void create_mesh();

  void setup_nodal_fields();
  void setup_edge_fields();
  void setup_element_fields();

  void setup_interior_algorithms();
  void setup_post_processing_algorithms();
  void setup_bc();
  void enforce_bc_on_exposed_faces();
  void setup_initial_conditions();
  void setup_property();
  void augment_property_map(
    PropertyIdentifier propID,
    ScalarFieldType *theField);

  void makeSureNodesHaveValidTopology();

  void initialize_global_variables();

  void balance_nodes();

  void create_output_mesh();
  void create_restart_mesh();
  void input_variables_from_mesh();

  void augment_output_variable_list(
      const std::string fieldName);
  
  void augment_restart_variable_list(
      std::string restartFieldName);

  void create_edges();
  void provide_entity_count();
  void delete_edges();
  void commit();

  void process_mesh_motion();
  void compute_centroid_on_parts(
    std::vector<std::string> partNames,
    std::vector<double> &centroid);

  void init_current_coordinates();

  std::string get_coordinates_name();
  bool has_mesh_motion() const;
  bool has_mesh_deformation() const;
  bool does_mesh_move() const;
  bool has_non_matching_boundary_face_alg() const;

  // overset boundary condition requires elemental field registration
  bool query_for_overset();
  
  void set_omega(
    stk::mesh::Part *targetPart,
    double omega);
  void set_current_displacement(
    stk::mesh::Part *targetPart,
    const std::vector<double> &centroidCoords,
    const std::vector<double> &unitVec);
  void set_current_coordinates(
    stk::mesh::Part *targetPart);
  void set_mesh_velocity(
    stk::mesh::Part *targetPart,
    const std::vector<double> &centroidCoords,
    const std::vector<double> &unitVec);
  void mesh_velocity_cross_product(double *o, double *c, double *u);

  // allow for a single mesh displacement (coords change, however, mesh motion fixed)
  void process_initial_displacement();
  void set_initial_displacement(
    stk::mesh::Part *targetPart,
    const std::vector<double> &centroidCoords,
    const std::vector<double> &unitVec,
    const double theAngle);

  // non-conformal-like algorithm suppoer
  void initialize_non_conformal();
  void initialize_overset();
  void initialize_post_processing_algorithms();

  void compute_geometry();
  void compute_vrtm();
  void compute_l2_scaling();
  void output_converged_results();
  void provide_output();
  void provide_restart_output();

  void register_interior_algorithm(
    stk::mesh::Part *part);

  void register_nodal_fields(
    stk::mesh::Part *part);

  void register_wall_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

  void register_inflow_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

  void register_open_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

  void register_symmetry_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

  void register_periodic_bc(
    stk::mesh::Part *masterMeshPart,
    stk::mesh::Part *slaveMeshPart,
    const double &searchTolerance,
    const std::string &searchMethodName);

  void setup_non_conformal_bc(
    stk::mesh::PartVector currentPartVec,
    stk::mesh::PartVector opposingPartVec,
    const NonConformalBoundaryConditionData &nonConformalBCData);

  void register_non_conformal_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

  void setup_overset_bc(
    const OversetBoundaryConditionData &oversetBCData);

  void periodic_field_update(
    stk::mesh::FieldBase *theField,
    const unsigned &sizeOfTheField,
    const bool &bypassFieldCheck = true) const;

  void periodic_delta_solution_update(
     stk::mesh::FieldBase *theField,
     const unsigned &sizeOfField) const;

  void periodic_max_field_update(
     stk::mesh::FieldBase *theField,
     const unsigned &sizeOfField) const;

  const stk::mesh::PartVector &get_slave_part_vector();

  void overset_orphan_node_field_update(
    stk::mesh::FieldBase *theField,
    const unsigned sizeRow,
    const unsigned sizeCol);

  virtual void populate_initial_condition();
  virtual void populate_boundary_data();
  virtual void boundary_data_to_state_data();
  virtual double populate_variables_from_input(const double currentTime);
  virtual void populate_external_variables_from_input(const double currentTime) {}
  virtual double populate_restart( double &timeStepNm1, int &timeStepCount);
  virtual void populate_derived_quantities();
  virtual void evaluate_properties();
  virtual double compute_adaptive_time_step();
  virtual void swap_states();
  virtual void predict_state();
  virtual void pre_timestep_work();
  virtual void output_banner();
  virtual void advance_time_step();
  virtual bool active_time_step();
  virtual void initial_work();
  
  void set_global_id();

  /** Initialize the HYPRE global row IDs
   *
   *  \sa Realm::hypreGlobalId_
   */
  void set_hypre_global_id();
 
  /// check job for fitting in memory
  void check_job(bool get_node_count);

  void dump_simulation_time();
  double provide_mean_norm();

  double get_hybrid_factor(
    const std::string dofname);
  double get_alpha_factor(
    const std::string dofname);
  double get_alpha_upw_factor(
    const std::string dofname);
  double get_upw_factor(
    const std::string dofname);
  bool primitive_uses_limiter(
    const std::string dofname);
  double get_lam_schmidt(
    const std::string dofname);
  double get_lam_prandtl(
    const std::string dofname, bool &prProvided);
  double get_turb_schmidt(
    const std::string dofname);
  double get_turb_prandtl(
    const std::string dofname);
  bool get_noc_usage(
    const std::string dofname);
  bool get_shifted_grad_op(
    const std::string dofname);
  bool get_skew_symmetric(
    const std::string dofname);
  double get_divU();

  // tanh factor specifics
  std::string get_tanh_functional_form(
    const std::string dofname);
  double get_tanh_trans(
    const std::string dofname);
  double get_tanh_width(
    const std::string dofname);

  // consistent mass matrix for projected nodal gradient
  bool get_consistent_mass_matrix_png(
    const std::string dofname);

  // pressure poisson nuance
  bool get_cvfem_shifted_mdot();
  bool get_cvfem_reduced_sens_poisson();
  
  bool has_nc_gauss_labatto_quadrature();
  bool get_nc_alg_upwind_advection();
  bool get_nc_alg_include_pstab();
  bool get_nc_alg_current_normal();

  void get_material_prop_eval(
    const PropertyIdentifier thePropID,
    std::vector<PropertyEvaluator*> &propEvalVec);

  bool is_turbulent();
  void is_turbulent(
    bool isIt);

  bool needs_enthalpy();
  void needs_enthalpy(bool needsEnthalpy);

  int number_of_states();

  std::string name();

  // redirection of stk::mesh::get_buckets to allow global selector
  stk::mesh::BucketVector const& get_buckets( 
    stk::mesh::EntityRank rank,
    const stk::mesh::Selector & selector) const;

  // get aura, bulk and meta data
  bool get_activate_aura();
  stk::mesh::BulkData & bulk_data();
  const stk::mesh::BulkData & bulk_data() const;
  stk::mesh::MetaData & meta_data();
  const stk::mesh::MetaData & meta_data() const;

  // inactive part
  stk::mesh::Selector get_inactive_selector();

  // push back equation to equation systems vector
  void push_equation_to_systems(
    EquationSystem *eqSystem);

  // provide all of the physics target names
  const std::vector<std::string> &get_physics_target_names();
  double get_tanh_blending(const std::string dofName);

  Realms& realms_;

  std::string name_;
  std::string type_;
  std::string inputDBName_;
  unsigned spatialDimension_;

  bool realmUsesEdges_;
  int solveFrequency_;
  bool isTurbulent_;
  bool needsEnthalpy_;

  double l2Scaling_;

  // ioBroker, meta and bulk data
  stk::mesh::MetaData *metaData_;
  stk::mesh::BulkData *bulkData_;
  stk::io::StkMeshIoBroker *ioBroker_;

  size_t resultsFileIndex_;
  size_t restartFileIndex_;

  // nalu field data
  GlobalIdFieldType *naluGlobalId_;

  // algorithm drivers managed by region
  ComputeGeometryAlgorithmDriver *computeGeometryAlgDriver_;
  AlgorithmDriver *errorIndicatorAlgDriver_;
  unsigned numInitialElements_;
  // for element, side, edge, node rank (node not used)
  stk::mesh::Selector adapterSelector_[4];
  Teuchos::RCP<stk::mesh::Selector> activePartForIO_;

  TimeIntegrator *timeIntegrator_;

  BoundaryConditions boundaryConditions_;
  InitialConditions initialConditions_;
  MaterialPropertys materialPropertys_;
  
  EquationSystems equationSystems_;

  double maxCourant_;
  double maxReynolds_;
  double targetCourant_;
  double timeStepChangeFactor_;
  int currentNonlinearIteration_;

  SolutionOptions *solutionOptions_;
  OutputInfo *outputInfo_;
  PostProcessingInfo *postProcessingInfo_;
  SolutionNormPostProcessing *solutionNormPostProcessing_;
  TurbulenceAveragingPostProcessing *turbulenceAveragingPostProcessing_;
  DataProbePostProcessing *dataProbePostProcessing_;
  Actuator *actuator_;

  std::vector<Algorithm *> propertyAlg_;
  std::map<PropertyIdentifier, ScalarFieldType *> propertyMap_;
  std::vector<Algorithm *> initCondAlg_;

  SizeType nodeCount_;
  bool estimateMemoryOnly_;
  double availableMemoryPerCoreGB_;
  double timerCreateMesh_;
  double timerPopulateMesh_;
  double timerPopulateFieldData_;
  double timerOutputFields_;
  double timerCreateEdges_;
  double timerNonconformal_;
  double timerInitializeEqs_;
  double timerPropertyEval_;
  double timerAdapt_;
  double timerTransferSearch_;
  double timerTransferExecute_;
  double timerSkinMesh_;
  double timerSortExposedFace_;

  NonConformalManager *nonConformalManager_;
  OversetManager *oversetManager_;
  bool hasNonConformal_;
  bool hasOverset_;

  // three type of transfer operations
  bool hasMultiPhysicsTransfer_;
  bool hasInitializationTransfer_;
  bool hasIoTransfer_;
  bool hasExternalDataTransfer_;

  PeriodicManager *periodicManager_;
  bool hasPeriodic_;
  bool hasFluids_;

  // global parameter list
  stk::util::ParameterList globalParameters_;

  // part for all exposed surfaces in the mesh
  stk::mesh::Part *exposedBoundaryPart_;

  // part for new edges
  stk::mesh::Part *edgesPart_;

  // cheack that all exposed surfaces have a bc applied
  bool checkForMissingBcs_;

  // check if there are negative Jacobians
  bool checkJacobians_;
  
  // types of physics
  bool isothermal_;
  bool uniform_;

  // some post processing of entity counts
  bool provideEntityCount_;

  // pointer to HDF5 file structure holding table
  HDF5FilePtr *HDF5ptr_;

  // automatic mesh decomposition; None, rib, rcb, multikl, etc.
  std::string autoDecompType_;

  // allow aura to be optional
  bool activateAura_;

  // allow detailed output (memory) to be provided
  bool activateMemoryDiagnostic_;

  // sometimes restarts can be missing states or dofs
  bool supportInconsistentRestart_;

  bool doBalanceNodes_;
  struct BalanceNodeOptions
  {
    BalanceNodeOptions() :
      target(1.0),
      numIters(5)
    {};

    double target;
    int numIters;
  };
  BalanceNodeOptions balanceNodeOptions_;

  // beginning wall time
  double wallTimeStart_;

  // mesh parts for all interior domains
  stk::mesh::PartVector interiorPartVec_;

  /** Vector holding side sets that have been registered with the boundary
   * conditions in the input file.
   *
   * The member is intended to for use in Realm::enforce_bc_on_exposed_faces to
   * check for "exposed surfaces" that might have not been assigned BCs in the
   * input file.
   *
   */
  stk::mesh::PartVector bcPartVec_;

  // empty part vector should it be required
  stk::mesh::PartVector emptyPartVector_;

  // base mesh parts
  stk::mesh::PartVector basePartVector_;

  std::vector<AuxFunctionAlgorithm *> bcDataAlg_;

  // transfer information; three types
  std::vector<Transfer *> multiPhysicsTransferVec_;
  std::vector<Transfer *> initializationTransferVec_;
  std::vector<Transfer *> ioTransferVec_;
  std::vector<Transfer *> externalDataTransferVec_;
  void augment_transfer_vector(Transfer *transfer, const std::string transferObjective, Realm *toRealm);
  void process_multi_physics_transfer(bool forcedXfer = false);
  void process_initialization_transfer();
  void process_io_transfer();
  void process_external_data_transfer();
  
  // process end of time step converged work
  void post_converged_work();

  // time information; calls through timeIntegrator
  double get_current_time();
  double get_time_step();
  double get_gamma1();
  double get_gamma2();
  double get_gamma3();
  int get_time_step_count() const;
  double get_time_step_from_file();
  bool get_is_fixed_time_step();
  bool get_is_terminate_based_on_time();
  double get_total_sim_time();
  int get_max_time_step_count();

  // restart
  bool restarted_simulation();
  bool support_inconsistent_restart();

  double get_stefan_boltzmann();
  double get_turb_model_constant(
    const TurbulenceModelConstant turbModelEnum);

  // id for the input mesh
  size_t inputMeshIdx_;

  // save off the node
  const YAML::Node & node_;

  // flag for CVFEM usage
  bool usesCVFEM_;

  std::string physics_part_name(std::string) const;
  std::vector<std::string> physics_part_names(std::vector<std::string>) const;
  std::string get_quad_type() const;

  // check for mesh changing
  bool mesh_changed() const;

  stk::mesh::PartVector allPeriodicInteractingParts_;
  stk::mesh::PartVector allNonConformalInteractingParts_;

  bool isFinalOuterIter_{false};

  /** The starting index (global) of the HYPRE linear system in this MPI rank
   *
   *  Note that this is actually the offset into the linear system. This index
   *  must be adjusted accordingly to account for multiple degrees of freedom on
   *  a particular node. This is performed in sierra::nalu::HypreLinearSystem.
   */
  stk::mesh::EntityId hypreILower_;

  /** The ending index (global) of the HYPRE linear system in this MPI rank
   *
   *  Note that this is actually the offset into the linear system. This index
   *  must be adjusted accordingly to account for multiple degrees of freedom on
   *  a particular node. This is performed in sierra::nalu::HypreLinearSystem.
   */
  stk::mesh::EntityId hypreIUpper_;

  /** The total number of HYPRE nodes in the linear system
   *
   *  Note that this is not an MPI rank local quantity
   */
  stk::mesh::EntityId hypreNumNodes_;

  /** Global Row IDs for the HYPRE linear system
   *
   *  The HYPRE IDs are different from STK IDs and Realm::naluGlobalId_ because
   *  HYPRE expects contiguous IDs for matrix rows and further requires that the
   *  IDs be ordered across MPI ranks; i.e., startIdx (MPI_rank + 1) =
   *  endIdx(MPI_rank) + 1.
   */
  HypreIDFieldType* hypreGlobalId_{nullptr};

  /** Flag indicating whether Hypre solver is being used for any of the equation
   * systems.
   */
  bool hypreIsActive_{false};

};

} // namespace nalu
} // namespace Sierra

#endif
