/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "gas_dynamics/GasDynamicsEquationSystem.h"
#include "AlgorithmDriver.h"
#include "AuxFunctionAlgorithm.h"
#include "ConstantAuxFunction.h"
#include "CopyFieldAlgorithm.h"
#include "EquationSystem.h"
#include "EquationSystems.h"
#include "Enums.h"
#include "FieldFunctions.h"
#include "NaluEnv.h"
#include "NaluParsing.h"
#include "Realm.h"
#include "Realms.h"
#include "Simulation.h"
#include "TimeIntegrator.h"

#include <property_evaluator/TemperaturePropAlgorithm.h>

// gas dynamics specifically
#include "gas_dynamics/AssembleGasDynamicsAlgorithmDriver.h"
#include "gas_dynamics/AssembleGasDynamicsFluxAlgorithm.h"
#include "gas_dynamics/AssembleGasDynamicsOpenAlgorithm.h"
#include "gas_dynamics/AssembleGasDynamicsSymmetryAlgorithm.h"
#include "gas_dynamics/AssembleGasDynamicsNonConformalAlgorithm.h"

// suplemental
#include "gas_dynamics/AssembleGasDynamicsCourantReynoldsElemAlgorithm.h"

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
// GasDynamicsEquationSystem - manages compressible (laminar) system
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
GasDynamicsEquationSystem::GasDynamicsEquationSystem(
  EquationSystems& eqSystems,
  bool debugOutput)
  : EquationSystem(eqSystems, "GasDynamicsEQS", "mixture_fraction"),
    density_(NULL),        // rho
    momentum_(NULL),       // rhoUj
    totalEnergy_(NULL),    // rhoE^t
    velocity_(NULL),       // uj
    totalEnthalpy_(NULL),  // rhoH^t = rho*(h + 1/2*uSq)
    staticEnthalpy_(NULL), // h = cpT
    pressure_(NULL),
    temperature_(NULL),
    machNumber_(NULL),
    speedOfSound_(NULL),
    cp_(NULL),
    cv_(NULL),
    viscosity_(NULL),
    thermalCond_(NULL),
    gamma_(NULL),
    dualNodalVolume_(NULL),
    rhsGasDyn_(NULL),
    assembleGasDynAlgDriver_(new AssembleGasDynamicsAlgorithmDriver(realm_)),
    cflReyAlgDriver_(new AlgorithmDriver(realm_)),
    isInit_(true),
    debugOutput_(debugOutput),
    fakeNorm_(0.0)
{
  // must be edge-based; AUSM+ based; no gradient extrapolation and limiting
  if ( !realm_.realmUsesEdges_ )
    throw std::runtime_error("GasDynamicsEquationSystem MUST activate edges");
  
  // push back EQ to manager
  realm_.push_equation_to_systems(this);
  
  // inform realm
  realm_.hasFluids_ = true;

  // advertise as gas dynamics; manages enthalpy (static) evaluation
  realm_.gasDynamics_ = true;

  // advertise as non-isothermal (still uniform)
  realm_.isothermal_ = false;
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
GasDynamicsEquationSystem::~GasDynamicsEquationSystem()
{
  delete assembleGasDynAlgDriver_;
  delete cflReyAlgDriver_;

  std::vector<Algorithm *>::iterator ii;
  for( ii=gasDynBcDataMapAlg_.begin(); ii!=gasDynBcDataMapAlg_.end(); ++ii )
    delete *ii;
  for( ii=densityAlg_.begin(); ii!=densityAlg_.end(); ++ii )
    delete *ii;
  for( ii=enthalpyAlg_.begin(); ii!=enthalpyAlg_.end(); ++ii )
    delete *ii;
}

//--------------------------------------------------------------------------
//-------- initial_work ----------------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::initial_work()
{
  // call base class method
  EquationSystem::initial_work();

  // proceed with a bunch of initial work; wrap in timer
  const double timeA = NaluEnv::self().nalu_time();

  // populate bc values for velocity and temperature to the nodal fields
  for ( size_t k = 0; k < gasDynBcDataMapAlg_.size(); ++k )
    gasDynBcDataMapAlg_[k]->execute();

  compute_density(); 
  compute_momentum();
  compute_static_enthalpy();
  compute_total_enthalpy();
  compute_total_energy();
  compute_speed_of_sound();
  compute_mach_number();

  // process CFL/Reynolds
  cflReyAlgDriver_->execute();

  // deal with state... just populated state Np1 above; copy state np1 to state n
  VectorFieldType &momN = momentum_->field_of_state(stk::mesh::StateN);
  VectorFieldType &momNp1 = momentum_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), momNp1, momN, realm_.get_activate_aura());
  ScalarFieldType &rhoN = density_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &rhoNp1 = density_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), rhoNp1, rhoN, realm_.get_activate_aura());
  ScalarFieldType &eN = totalEnergy_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &eNp1 = totalEnergy_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), eNp1, eN, realm_.get_activate_aura());

  if ( debugOutput_ )
    dump_state("GasDynamicsEquationSystem::initial_work(): post");
  
  const double timeB = NaluEnv::self().nalu_time();
  timerMisc_ += (timeB-timeA);
}

//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::register_nodal_fields(
  stk::mesh::Part *part)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();
  const int numStates = realm_.number_of_states();

  // dofs first
  density_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "density", numStates));
  stk::mesh::put_field_on_mesh(*density_, *part, nullptr);
  momentum_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "momentum", numStates));
  stk::mesh::put_field_on_mesh(*momentum_, *part, nDim, nullptr);
  stk::io::set_field_output_type(*momentum_, stk::io::FieldOutputType::VECTOR_3D);
  totalEnergy_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "total_energy", numStates));
  stk::mesh::put_field_on_mesh(*totalEnergy_, *part, nullptr);

  // aux
  velocity_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "velocity"));
  stk::mesh::put_field_on_mesh(*velocity_, *part, nDim, nullptr);
  stk::io::set_field_output_type(*velocity_, stk::io::FieldOutputType::VECTOR_3D);
  totalEnthalpy_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "total_enthalpy"));
  stk::mesh::put_field_on_mesh(*totalEnthalpy_, *part, nullptr);
  staticEnthalpy_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "enthalpy"));
  stk::mesh::put_field_on_mesh(*staticEnthalpy_, *part, nullptr);
  pressure_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "pressure"));
  stk::mesh::put_field_on_mesh(*pressure_, *part, nullptr);
  temperature_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "temperature"));
  stk::mesh::put_field_on_mesh(*temperature_, *part, nullptr);
  machNumber_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "mach_number"));
  stk::mesh::put_field_on_mesh(*machNumber_, *part, nullptr);
  speedOfSound_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "speed_of_sound"));
  stk::mesh::put_field_on_mesh(*speedOfSound_, *part, nullptr);

  // properties: user specifies: mu, R, MW, gamma, and kappa
  //             property evaluator populates nodal fields for 
  //             mu, gamma, kappa, cp and cv (each based on gamma and R/MW)
  viscosity_ =  &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "viscosity"));
  stk::mesh::put_field_on_mesh(*viscosity_, *part, nullptr);
  cp_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "specific_heat"));
  stk::mesh::put_field_on_mesh(*cp_, *part, nullptr);
  cv_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "specific_heat_v"));
  stk::mesh::put_field_on_mesh(*cv_, *part, nullptr);
  thermalCond_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "thermal_conductivity"));
  stk::mesh::put_field_on_mesh(*thermalCond_, *part, nullptr);
  gamma_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "gamma"));
  stk::mesh::put_field_on_mesh(*gamma_, *part, nullptr);

  // push to property list: density needs to be populated for IC
  realm_.augment_property_map(DENSITY_ID, density_);
  realm_.augment_property_map(VISCOSITY_ID, viscosity_);
  realm_.augment_property_map(THERMAL_COND_ID, thermalCond_);
  // Cv and Cp derived from gamma and R/MW
  realm_.augment_property_map(GAMMA_ID, gamma_);

  // dual nodal volume (should push up...)
  dualNodalVolume_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "dual_nodal_volume"));
  stk::mesh::put_field_on_mesh(*dualNodalVolume_, *part, nullptr);

  // residual; special ordering... Pick momentum first, then continuity, then total energy
  rhsGasDyn_ = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "rhs_gas_dynamics"));
  stk::mesh::put_field_on_mesh(*rhsGasDyn_, *part, nDim+2, nullptr);

  // fileds that require restart
  realm_.augment_restart_variable_list("density");
  realm_.augment_restart_variable_list("momentum");
  realm_.augment_restart_variable_list("total_energy");
  realm_.augment_restart_variable_list("pressure");
  
  // make sure all states are properly populated (restart can handle this)
  if ( numStates > 2 ) {
    throw std::runtime_error("BDF2 implicit integrator not supported");
  }
}

//--------------------------------------------------------------------------
//-------- register_edge_fields --------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::register_edge_fields(
  stk::mesh::Part *part)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();
  const int nDim = meta_data.spatial_dimension();
  VectorFieldType *edgeAreaVec 
    = &(meta_data.declare_field<double>(stk::topology::EDGE_RANK, "edge_area_vector"));
  stk::mesh::put_field_on_mesh(*edgeAreaVec, *part, nDim, nullptr);
  stk::io::set_field_output_type(*edgeAreaVec, stk::io::FieldOutputType::VECTOR_3D);
}

//--------------------------------------------------------------------------
//-------- register_element_fields -----------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::register_element_fields(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{
  stk::mesh::MetaData &meta_data = realm_.meta_data();

  // register the intersected elemental field
  if ( realm_.query_for_overset() ) {
    const int sizeOfElemField = 1;
    GenericFieldType *intersectedElement
      = &(meta_data.declare_field<double>(stk::topology::ELEMENT_RANK, "intersected_element"));
    stk::mesh::put_field_on_mesh(*intersectedElement, *part, sizeOfElemField, nullptr);
  }

  // provide mean element Peclet and Courant fields; always...
  GenericFieldType *elemReynolds
    = &(meta_data.declare_field<double>(stk::topology::ELEMENT_RANK, "element_reynolds"));
  stk::mesh::put_field_on_mesh(*elemReynolds, *part, nullptr);
  GenericFieldType *elemCourant
    = &(meta_data.declare_field<double>(stk::topology::ELEMENT_RANK, "element_courant"));
  stk::mesh::put_field_on_mesh(*elemCourant, *part, nullptr);
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::register_interior_algorithm(
  stk::mesh::Part *part)
{
  // algorithm type supprted
  const AlgorithmType algI = INTERIOR;

  // edge-based (only) to compute RHS residual; no mass term, no source terms
  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleGasDynAlgDriver_->algMap_.find(algI);
  if ( it == assembleGasDynAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg = new AssembleGasDynamicsFluxAlgorithm(realm_, part, density_, 
                                                             momentum_, velocity_, totalEnthalpy_, 
                                                             pressure_, temperature_, speedOfSound_, viscosity_,
                                                             thermalCond_, rhsGasDyn_);
    assembleGasDynAlgDriver_->algMap_[algI] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

  // gas dynamics manages population of DOFs and certain properties

  // extract material prop evaluation for density and create alg to compute it
  std::vector<PropertyEvaluator *> rhoEvalVec;
  realm_.get_material_prop_eval(DENSITY_ID, rhoEvalVec);  
  TemperaturePropAlgorithm *rhoAlg
    = new TemperaturePropAlgorithm( realm_, part, density_, rhoEvalVec[0], temperature_->name());
  densityAlg_.push_back(rhoAlg);

  // extract material prop evaluation for enthalpy and create alg to compute it
  std::vector<PropertyEvaluator *> enthEvalVec;
  realm_.get_material_prop_eval(ENTHALPY_ID, enthEvalVec);  
  TemperaturePropAlgorithm *enthAlg
    = new TemperaturePropAlgorithm( realm_, part, staticEnthalpy_, enthEvalVec[0], temperature_->name());
  enthalpyAlg_.push_back(enthAlg);

  // non-solver Courant number alg
  it = cflReyAlgDriver_->algMap_.find(algI);
  if ( it == cflReyAlgDriver_->algMap_.end() ) {
    AssembleGasDynamicsCourantReynoldsElemAlgorithm *theAlg
      = new AssembleGasDynamicsCourantReynoldsElemAlgorithm(realm_, part);
    cflReyAlgDriver_->algMap_[algI] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }
  // GDFIXME: Add speed of sound of Courant
}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const InflowBoundaryConditionData &inflowBCData)
{
  // dirichlet mesh part
  dirichletPart_.push_back(part);

  stk::mesh::MetaData &meta_data = realm_.meta_data();
  const unsigned nDim = meta_data.spatial_dimension();

  // general variables that will be used
  InflowUserData userData = inflowBCData.userData_;
  UserDataType theDataType = CONSTANT_UD;
  std::string primitiveName = "na";

  // register boundary data; velocity
  primitiveName = "velocity";
  VectorFieldType *uBcField = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "velocity_bc"));
  stk::mesh::put_field_on_mesh(*uBcField, *part, nDim, nullptr);
  stk::io::set_field_output_type(*uBcField, stk::io::FieldOutputType::VECTOR_3D);

  theDataType = get_bc_data_type(userData, primitiveName);
  if ( CONSTANT_UD != theDataType )
    throw std::runtime_error("only constant");

  Velocity ux = userData.u_;
  std::vector<double> uSpec(nDim);
  uSpec[0] = ux.ux_;
  uSpec[1] = ux.uy_;
  if ( nDim > 2)
    uSpec[2] = ux.uz_;
    
  // new it
  AuxFunction *uAuxFunc = new ConstantAuxFunction(0, nDim, uSpec);

  // bc data alg
  AuxFunctionAlgorithm *uauxAlg
    = new AuxFunctionAlgorithm(realm_, part,
			       uBcField, uAuxFunc,
			       stk::topology::NODE_RANK);
  
  realm_.initCondAlg_.push_back(uauxAlg);

  // copy bc 
  CopyFieldAlgorithm *theCopyAlgU
    = new CopyFieldAlgorithm(realm_, part,
                             uBcField, velocity_,
                             0, nDim,
                             stk::topology::NODE_RANK);
  gasDynBcDataMapAlg_.push_back(theCopyAlgU);

  // register boundary data; temperature
  primitiveName = "temperature";
  ScalarFieldType *tBcField = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "temperature_bc"));
  stk::mesh::put_field_on_mesh(*tBcField, *part, nullptr);
  
  theDataType = get_bc_data_type(userData, primitiveName);
  if ( CONSTANT_UD != theDataType )
    throw std::runtime_error("only constant");

  Temperature theTemp = userData.temperature_;
  std::vector<double> tSpec(1);
  tSpec[0] = theTemp.temperature_;
 
  // new it
  AuxFunction *tAuxFunc = new ConstantAuxFunction(0, 1, tSpec);

  // bc data alg
  AuxFunctionAlgorithm *tauxAlg
    = new AuxFunctionAlgorithm(realm_, part,
			       tBcField, tAuxFunc,
			       stk::topology::NODE_RANK);
  
  realm_.initCondAlg_.push_back(tauxAlg);

  // copy bc 
  CopyFieldAlgorithm *theCopyAlgT
    = new CopyFieldAlgorithm(realm_, part,
                             tBcField, temperature_,
                             0, 1,
                             stk::topology::NODE_RANK);
  gasDynBcDataMapAlg_.push_back(theCopyAlgT);
}

//--------------------------------------------------------------------------
//-------- register_open_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &partTopo,
  const OpenBoundaryConditionData &openBCData)
{
  // algorithm type
  const AlgorithmType algType = OPEN;
  
  stk::mesh::MetaData &meta_data = realm_.meta_data();
  const unsigned nDim = meta_data.spatial_dimension();

  // general variables that will be used
  OpenUserData userData = openBCData.userData_;
  UserDataType theDataType = CONSTANT_UD;
  std::string primitiveName = "na";

  // register boundary data; velocity
  primitiveName = "velocity";
  VectorFieldType *uBcField = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "velocity_bc"));
  stk::mesh::put_field_on_mesh(*uBcField, *part, nDim, nullptr);
  stk::io::set_field_output_type(*uBcField, stk::io::FieldOutputType::VECTOR_3D);

  theDataType = get_bc_data_type(userData, primitiveName);
  if ( CONSTANT_UD != theDataType )
    throw std::runtime_error("only constant");

  Velocity ux = userData.u_;
  std::vector<double> uSpec(nDim);
  uSpec[0] = ux.ux_;
  uSpec[1] = ux.uy_;
  if ( nDim > 2)
    uSpec[2] = ux.uz_;
    
  // new it
  AuxFunction *uAuxFunc = new ConstantAuxFunction(0, nDim, uSpec);

  // bc data alg
  AuxFunctionAlgorithm *uauxAlg
    = new AuxFunctionAlgorithm(realm_, part,
			       uBcField, uAuxFunc,
			       stk::topology::NODE_RANK);
  
  realm_.initCondAlg_.push_back(uauxAlg);

  // register boundary data; temperature
  primitiveName = "temperature";
  ScalarFieldType *tBcField = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "temperature_bc"));
  stk::mesh::put_field_on_mesh(*tBcField, *part, nullptr);
  
  theDataType = get_bc_data_type(userData, primitiveName);
  if ( CONSTANT_UD != theDataType )
    throw std::runtime_error("only constant");

  Temperature theTemp = userData.temperature_;
  std::vector<double> tSpec(1);
  tSpec[0] = theTemp.temperature_;
 
  // new it
  AuxFunction *tAuxFunc = new ConstantAuxFunction(0, 1, tSpec);

  // bc data alg
  AuxFunctionAlgorithm *tauxAlg
    = new AuxFunctionAlgorithm(realm_, part,
			       tBcField, tAuxFunc,
			       stk::topology::NODE_RANK);
  
  realm_.initCondAlg_.push_back(tauxAlg);

  // register boundary data; pressure
  primitiveName = "pressure";
  ScalarFieldType *pBcField = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "pressure_bc"));
  stk::mesh::put_field_on_mesh(*pBcField, *part, nullptr);
  
  theDataType = get_bc_data_type(userData, primitiveName);
  if ( CONSTANT_UD != theDataType )
    throw std::runtime_error("only constant");

  Pressure thePress = userData.p_;
  std::vector<double> pSpec(1);
  pSpec[0] = thePress.pressure_;
 
  // new it
  AuxFunction *pAuxFunc = new ConstantAuxFunction(0, 1, pSpec);

  // bc data alg
  AuxFunctionAlgorithm *pauxAlg
    = new AuxFunctionAlgorithm(realm_, part,
			       pBcField, pAuxFunc,
			       stk::topology::NODE_RANK);
  
  realm_.initCondAlg_.push_back(pauxAlg);

  // outflow (convected out only)
  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleGasDynAlgDriver_->algMap_.find(algType);
  if ( it == assembleGasDynAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg = new AssembleGasDynamicsOpenAlgorithm(realm_, part, 
                                                             density_, momentum_, totalEnthalpy_, 
                                                             pBcField, tBcField, cp_, gamma_, 
                                                             rhsGasDyn_);
    assembleGasDynAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }  
}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const WallBoundaryConditionData &wallBCData)
{
  // dirichlet mesh part; no wall function
  dirichletPart_.push_back(part);

  stk::mesh::MetaData &meta_data = realm_.meta_data();
  const unsigned nDim = meta_data.spatial_dimension();

  // general variables that will be used
  WallUserData userData = wallBCData.userData_;
  UserDataType theDataType = CONSTANT_UD;
  std::string primitiveName = "na";

  // register boundary data; velocity
  primitiveName = "velocity";
  VectorFieldType *uBcField = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "velocity_bc"));
  stk::mesh::put_field_on_mesh(*uBcField, *part, nDim, nullptr);
  stk::io::set_field_output_type(*uBcField, stk::io::FieldOutputType::VECTOR_3D);

  theDataType = get_bc_data_type(userData, primitiveName);
  if ( CONSTANT_UD != theDataType )
    throw std::runtime_error("only constant");

  Velocity ux = userData.u_;
  std::vector<double> uSpec(nDim);
  uSpec[0] = ux.ux_;
  uSpec[1] = ux.uy_;
  if ( nDim > 2)
    uSpec[2] = ux.uz_;
  
  // new it
  AuxFunction *uAuxFunc = new ConstantAuxFunction(0, nDim, uSpec);

  // bc data alg
  AuxFunctionAlgorithm *uauxAlg
    = new AuxFunctionAlgorithm(realm_, part,
			       uBcField, uAuxFunc,
			       stk::topology::NODE_RANK);
  
  realm_.initCondAlg_.push_back(uauxAlg);

  // copy bc 
  CopyFieldAlgorithm *theCopyAlgU
    = new CopyFieldAlgorithm(realm_, part,
                             uBcField, velocity_,
                             0, nDim,
                             stk::topology::NODE_RANK);
  gasDynBcDataMapAlg_.push_back(theCopyAlgU);

  // register boundary data; temperature; either temperature specified or adiabatic...
  primitiveName = "temperature";
  if ( bc_data_specified(userData, primitiveName) ) {
    
    ScalarFieldType *tBcField = &(meta_data.declare_field<double>(stk::topology::NODE_RANK, "temperature_bc"));
    stk::mesh::put_field_on_mesh(*tBcField, *part, nullptr);
    
    theDataType = get_bc_data_type(userData, primitiveName);
    if ( CONSTANT_UD != theDataType )
      throw std::runtime_error("only constant");
    
    Temperature theTemp = userData.temperature_;
    std::vector<double> tSpec(1);
    tSpec[0] = theTemp.temperature_;
    
    // new it
    AuxFunction *tAuxFunc = new ConstantAuxFunction(0, 1, tSpec);
    
    // bc data alg
    AuxFunctionAlgorithm *tauxAlg
      = new AuxFunctionAlgorithm(realm_, part,
                                 tBcField, tAuxFunc,
                                 stk::topology::NODE_RANK);
    
    realm_.initCondAlg_.push_back(tauxAlg);

    // copy bc 
    CopyFieldAlgorithm *theCopyAlgT
      = new CopyFieldAlgorithm(realm_, part,
                               tBcField, temperature_,
                               0, 1,
                               stk::topology::NODE_RANK);
    gasDynBcDataMapAlg_.push_back(theCopyAlgT);
  }
  else {
    NaluEnv::self().naluOutputP0() << "adiabatic wall" << std::endl;
  }
}

//--------------------------------------------------------------------------
//-------- register_symmetry_bc --------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/,
  const SymmetryBoundaryConditionData &/*wallBCData*/)
{
  // algorithm type
  const AlgorithmType algType = SYMMETRY;

  // symmetry (normal stress)
  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleGasDynAlgDriver_->algMap_.find(algType);
  if ( it == assembleGasDynAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg = new AssembleGasDynamicsSymmetryAlgorithm(realm_, part, pressure_, 
                                                                 rhsGasDyn_);
    assembleGasDynAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_non_conformal_bc ---------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::register_non_conformal_bc(
  stk::mesh::Part *part,
  const stk::topology &/*theTopo*/)
{
  const AlgorithmType algType = NON_CONFORMAL;
  
  // edge-based (only) to compute RHS residual; no mass term, no source terms
  std::map<AlgorithmType, Algorithm *>::iterator it
    = assembleGasDynAlgDriver_->algMap_.find(algType);
  if ( it == assembleGasDynAlgDriver_->algMap_.end() ) {
    Algorithm *theAlg = new AssembleGasDynamicsNonConformalAlgorithm(
     realm_, part, density_, 
     momentum_, velocity_, totalEnthalpy_, 
     pressure_, temperature_, speedOfSound_, viscosity_,
     thermalCond_, rhsGasDyn_);
    assembleGasDynAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }
}

//--------------------------------------------------------------------------
//-------- register_overset_bc ---------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::register_overset_bc()
{
  const int nDim = realm_.meta_data().spatial_dimension();

  UpdateOversetFringeAlgorithmDriver* theAlg = new UpdateOversetFringeAlgorithmDriver(realm_);
  // Perform fringe updates before all equation system solves
  equationSystems_.preIterAlgDriver_.push_back(theAlg);
  theAlg->fields_.push_back(
    std::unique_ptr<OversetFieldData>(new OversetFieldData(density_,1,1)));
  theAlg->fields_.push_back(
    std::unique_ptr<OversetFieldData>(new OversetFieldData(momentum_,1,nDim)));
  theAlg->fields_.push_back(
    std::unique_ptr<OversetFieldData>(new OversetFieldData(totalEnergy_,1,1)));

  if ( realm_.has_mesh_motion() ) {
    UpdateOversetFringeAlgorithmDriver* theAlgPost = new UpdateOversetFringeAlgorithmDriver(realm_,false);
    // Perform fringe updates after all equation system solves (ideally on the post_time_step)
    equationSystems_.postIterAlgDriver_.push_back(theAlgPost);
    theAlgPost->fields_.push_back(
      std::unique_ptr<OversetFieldData>(new OversetFieldData(density_,1,1)));
    theAlgPost->fields_.push_back(
      std::unique_ptr<OversetFieldData>(new OversetFieldData(momentum_,1,nDim)));
    theAlgPost->fields_.push_back(
      std::unique_ptr<OversetFieldData>(new OversetFieldData(totalEnergy_,1,1)));
  }
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::initialize()
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- reinitialize_linear_system --------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::reinitialize_linear_system()
{
  // nothing to do
}


//--------------------------------------------------------------------------
//-------- predict_state ---------------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::predict_state()
{
  // copy old state to new; this allows all of the explict algs to operate on n+1 state..
  VectorFieldType &momN = momentum_->field_of_state(stk::mesh::StateN);
  VectorFieldType &momNp1 = momentum_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), momN, momNp1, realm_.get_activate_aura());

  ScalarFieldType &rhoN = density_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &rhoNp1 = density_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), rhoN, rhoNp1, realm_.get_activate_aura());

  ScalarFieldType &eN = totalEnergy_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &eNp1 = totalEnergy_->field_of_state(stk::mesh::StateNP1);
  field_copy(realm_.meta_data(), realm_.bulk_data(), eN, eNp1, realm_.get_activate_aura());
}

//--------------------------------------------------------------------------
//-------- solve_and_update ------------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::solve_and_update()
{
  // initialization 
  if ( isInit_ ) {
    isInit_ = false;
  }

  // explicit approach
  double timeA = NaluEnv::self().nalu_time();
  assemble_gas_dynamics();
  update_gas_dynamics();

  // process CFL/Reynolds
  cflReyAlgDriver_->execute();

  double timeB = NaluEnv::self().nalu_time();
  timerAssemble_ += (timeB-timeA);  
}

//--------------------------------------------------------------------------
//-------- assemble_gas_dynamics -------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::assemble_gas_dynamics()
{
  if ( debugOutput_ )
    dump_state("GasDynamicsEquationSystem::assemble_gas_dynamics(): pre");
  
  // execute the explicit alg; pre and post manage zeroing and specialized assemblies
  assembleGasDynAlgDriver_->execute();
}

//--------------------------------------------------------------------------
//-------- update_gas_dynamics ---------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::update_gas_dynamics()
{
  NaluEnv::self().naluOutputP0() << "GasDynamicsEquationSystem::update_gas_dynamics" << std::endl;
 
  const int nDim = realm_.meta_data().spatial_dimension();
  const int cOffset = nDim;
  const int eOffset = nDim + 1;
  const int totalSize = nDim + 2;

  const double dt = realm_.get_time_step();

  // toggle for debug
  const double updateFac = 1.0;

  // fake norm for regression testing
  double l_fakeNorm = 0.0;

  // extract fields of state
  VectorFieldType &momentumN = momentum_->field_of_state(stk::mesh::StateN);
  VectorFieldType &momentumNp1 = momentum_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &densityN = density_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);
  ScalarFieldType &totalEnergyN = totalEnergy_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &totalEnergyNp1 = totalEnergy_->field_of_state(stk::mesh::StateNP1);

  // select everything other than dirichlet parts
  stk::mesh::Selector s_nodes = (!stk::mesh::selectUnion(dirichletPart_)) 
    & stk::mesh::selectField(*rhsGasDyn_)
    & !(realm_.get_inactive_selector());

  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    const double * dualNodalVolume = stk::mesh::field_data(*dualNodalVolume_, b);
    const double * rhsGasDyn = stk::mesh::field_data(*rhsGasDyn_, b);
    double * momN = stk::mesh::field_data(momentumN, b);
    double * momNp1 = stk::mesh::field_data(momentumNp1, b);
    double * rhoN = stk::mesh::field_data(densityN, b);
    double * rhoNp1 = stk::mesh::field_data(densityNp1, b);
    double * teN = stk::mesh::field_data(totalEnergyN, b);
    double * teNp1 = stk::mesh::field_data(totalEnergyNp1, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      
      const int kTotalS = k*totalSize;
      const int kNdim = k*nDim;

      const double fac = updateFac*dt/dualNodalVolume[k];

      // momentum first
      for ( int j = 0; j < nDim; ++j ) {
        momNp1[kNdim+j] = momN[kNdim+j] + fac*rhsGasDyn[kTotalS+j];
      }
      
      // continuity and energy
      rhoNp1[k] = rhoN[k] + fac*rhsGasDyn[kTotalS+cOffset];
      teNp1[k] = teN[k] + fac*rhsGasDyn[kTotalS+eOffset];

      l_fakeNorm += rhsGasDyn[kTotalS+cOffset]*rhsGasDyn[kTotalS+cOffset];
    }
  }
  
  // overset update; reconstruct the dofs and let aux variables go for the ride below
  if ( realm_.hasOverset_ ) {
    realm_.overset_constraint_node_field_update(density_, 1, nDim);
    realm_.overset_constraint_node_field_update(momentum_, 1, 1);
    realm_.overset_constraint_node_field_update(totalEnergy_, 1, 1);
  }
  
  // update aux variables
  compute_velocity();
  compute_pressure();
  compute_temperature();
  compute_static_enthalpy();
  compute_total_enthalpy();
  compute_speed_of_sound();
  compute_mach_number();
  
  if ( debugOutput_ )
    dump_state("GasDynamicsEquationSystem::assemble_gas_dynamics(): post");

  // compute fake norm sum
  double g_fakeNorm = 0.0;
  stk::ParallelMachine comm = NaluEnv::self().parallel_comm();
  stk::all_reduce_sum(comm, &l_fakeNorm, &g_fakeNorm, 1);

  // advertise
  fakeNorm_ = std::sqrt(g_fakeNorm);
}

//--------------------------------------------------------------------------
//-------- provide_scaled_norm ---------------------------------------------
//--------------------------------------------------------------------------
double
GasDynamicsEquationSystem::provide_scaled_norm()
{
  return fakeNorm_;
}

//--------------------------------------------------------------------------
//-------- provide_norm ----------------------------------------------------
//--------------------------------------------------------------------------
double
GasDynamicsEquationSystem::provide_norm()
{
  return fakeNorm_;
}

//--------------------------------------------------------------------------
//-------- dump_state ------------------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::dump_state(
  const std::string indicator)
{
  NaluEnv::self().naluOutputP0() << indicator << std::endl;

  stk::mesh::MetaData &metaData = realm_.meta_data();
  const int nDim = metaData.spatial_dimension();

  // extract coords
  VectorFieldType *coords_ 
    = metaData.get_field<double>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
  
  // deal with state
  ScalarFieldType &densityN = density_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  VectorFieldType &momentumN = momentum_->field_of_state(stk::mesh::StateN);
  VectorFieldType &momentumNp1 = momentum_->field_of_state(stk::mesh::StateNP1);

  ScalarFieldType &totalEnergyN = totalEnergy_->field_of_state(stk::mesh::StateN);
  ScalarFieldType &totalEnergyNp1 = totalEnergy_->field_of_state(stk::mesh::StateNP1);

  stk::mesh::Selector s_nodes = stk::mesh::selectField(*momentum_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    const double * coords = stk::mesh::field_data(*coords_, b);

    // state, old and new
    const double * rN = stk::mesh::field_data(densityN, b);
    const double * mN = stk::mesh::field_data(momentumN, b);
    const double * teN = stk::mesh::field_data(totalEnergyN, b);
    const double * rNp1 = stk::mesh::field_data(densityNp1, b);
    const double * mNp1 = stk::mesh::field_data(momentumNp1, b);
    const double * teNp1 = stk::mesh::field_data(totalEnergyNp1, b);

    // aux
    const double * velocity = stk::mesh::field_data(*velocity_, b);
    const double * totalEnthalpy = stk::mesh::field_data(*totalEnthalpy_, b);
    const double * staticEnthalpy = stk::mesh::field_data(*staticEnthalpy_, b);
    const double * pressure = stk::mesh::field_data(*pressure_, b);
    const double * temperature = stk::mesh::field_data(*temperature_, b);
    const double * machNumber = stk::mesh::field_data(*machNumber_, b);
    const double * speedOfSound = stk::mesh::field_data(*speedOfSound_, b);
    const double * cp = stk::mesh::field_data(*cp_, b);
    const double * cv = stk::mesh::field_data(*cv_, b);
    const double * viscosity = stk::mesh::field_data(*viscosity_, b);
    const double * thermalCond = stk::mesh::field_data(*thermalCond_, b);
    const double * gamma = stk::mesh::field_data(*gamma_, b);
    const double * dualNodalVolume = stk::mesh::field_data(*dualNodalVolume_, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int kNdim = k*nDim;
      NaluEnv::self().naluOutput() << "==================================" << std::endl;
      for ( int i = 0; i < nDim; ++i ) {
        NaluEnv::self().naluOutput() << "Coords:       " << coords[kNdim+i] << std::endl;
        NaluEnv::self().naluOutput() << "MomentumN:    " << mN[kNdim+i] << std::endl;
        NaluEnv::self().naluOutput() << "MomentumNp1:  " << mNp1[kNdim+i] << std::endl;
        NaluEnv::self().naluOutput() << "Velocity:     " << velocity[kNdim+i] << std::endl;
      }
   
       NaluEnv::self().naluOutput() << "DensityN:     " << rN[k] << std::endl;
       NaluEnv::self().naluOutput() << "DensityNp1:   " << rNp1[k] << std::endl;
       NaluEnv::self().naluOutput() << "totalEnerN:   " << teN[k] << std::endl;
       NaluEnv::self().naluOutput() << "totalEnerNp1: " << teNp1[k] << std::endl;
       NaluEnv::self().naluOutput() << "totalEnth:    " << totalEnthalpy[k] << std::endl;
       NaluEnv::self().naluOutput() << "staticEnth:   " << staticEnthalpy[k] << std::endl;
       NaluEnv::self().naluOutput() << "pressure:     " << pressure[k] << std::endl;
       NaluEnv::self().naluOutput() << "temperature:  " << temperature[k] << std::endl;
       NaluEnv::self().naluOutput() << "machNumber:   " << machNumber[k] << std::endl;
       NaluEnv::self().naluOutput() << "speedOfSound: " << speedOfSound[k] << std::endl;
       NaluEnv::self().naluOutput() << "cp:           " << cp[k] << std::endl;
       NaluEnv::self().naluOutput() << "cv:           " << cv[k] << std::endl;
       NaluEnv::self().naluOutput() << "viscosity:    " << viscosity[k] << std::endl;
       NaluEnv::self().naluOutput() << "thermalCond:  " << thermalCond[k] << std::endl;
       NaluEnv::self().naluOutput() << "gamma:        " << gamma[k] << std::endl;
       NaluEnv::self().naluOutput() << "dualNodalVol: " << dualNodalVolume[k] << std::endl;
    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_density -------------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::compute_density()
{
  // ideal gas law; rho = P*M/R/T
  for ( size_t k = 0; k < densityAlg_.size(); ++k )
    densityAlg_[k]->execute();
}

//--------------------------------------------------------------------------
//-------- compute_static_enthalpy -----------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::compute_static_enthalpy()
{
  // should be cp*dT; with Tref of zero
  for ( size_t k = 0; k < enthalpyAlg_.size(); ++k )
    enthalpyAlg_[k]->execute();  
}

//--------------------------------------------------------------------------
//-------- compute_momentum ------------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::compute_momentum()
{
  // rhoUj = rho*uj
  const int nDim = realm_.meta_data().spatial_dimension();

  stk::mesh::Selector s_nodes = stk::mesh::selectField(*momentum_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    const double * velocity = stk::mesh::field_data(*velocity_, b);
    const double * rho = stk::mesh::field_data(*density_, b);
    double * momentum = stk::mesh::field_data(*momentum_, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const double rhoK = rho[k];
      const int kNdim = k*nDim;
      for ( int j = 0; j < nDim; ++j ) {
        momentum[kNdim+j] = rhoK*velocity[kNdim+j];
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_total_energy --------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::compute_total_energy()
{
  // rhoE^t = p/(gamma - 1) + 1/2*rho*uSq = rho*(h + 1/2*rho*uSq) - p 
  const int nDim = realm_.meta_data().spatial_dimension();

  stk::mesh::Selector s_nodes = stk::mesh::selectField(*totalEnergy_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * totalEnergy = stk::mesh::field_data(*totalEnergy_, b);
    const double * pressure = stk::mesh::field_data(*pressure_, b);
    const double * rho = stk::mesh::field_data(*density_, b);
    const double * velocity = stk::mesh::field_data(*velocity_, b);
    const double * gamma = stk::mesh::field_data(*gamma_, b);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int kNdim = k*nDim;
      double uSq = 0.0;
      for ( int j = 0; j < nDim; ++j ) {
        uSq += velocity[kNdim+j]*velocity[kNdim+j];
      }
      totalEnergy[k] = pressure[k]/(gamma[k]-1.0) + 0.5*rho[k]*uSq;
    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_total_enthalpy ------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::compute_total_enthalpy()
{
  // rhoH^t = rho*(h + 1/2 rho*uSq)
  const int nDim = realm_.meta_data().spatial_dimension();

  stk::mesh::Selector s_nodes = stk::mesh::selectField(*totalEnthalpy_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    double * totalEnthalpy = stk::mesh::field_data(*totalEnthalpy_, b);
    const double * rho = stk::mesh::field_data(*density_, b);
    const double * velocity = stk::mesh::field_data(*velocity_, b);
    const double * staticEnthalpy = stk::mesh::field_data(*staticEnthalpy_, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int kNdim = k*nDim;
      double uSq = 0.0;
      for ( int j = 0; j < nDim; ++j ) {
        uSq += velocity[kNdim+j]*velocity[kNdim+j];
      }
      totalEnthalpy[k] = rho[k]*(staticEnthalpy[k] + 0.5*uSq);
    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_velocity ------------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::compute_velocity()
{
  // uj = rhoUj/rho
  const int nDim = realm_.meta_data().spatial_dimension();

  stk::mesh::Selector s_nodes = stk::mesh::selectField(*velocity_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    const double * momentum = stk::mesh::field_data(*momentum_, b);
    const double * rho = stk::mesh::field_data(*density_, b);
    double * velocity = stk::mesh::field_data(*velocity_, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const double inv_rhoK = 1.0/rho[k];
      const int kNdim = k*nDim;
      for ( int j = 0; j < nDim; ++j ) {
        velocity[kNdim+j] = momentum[kNdim+j]*inv_rhoK;
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_pressure ------------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::compute_pressure()
{
  // p = (gamma - 1)*[rhoE^t - 1/2 rho*uSq ]
  const int nDim = realm_.meta_data().spatial_dimension();

  stk::mesh::Selector s_nodes = stk::mesh::selectField(*pressure_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    const double * totalEnergy = stk::mesh::field_data(*totalEnergy_, b);
    const double * rho = stk::mesh::field_data(*density_, b);
    const double * velocity = stk::mesh::field_data(*velocity_, b);
    const double * gamma = stk::mesh::field_data(*gamma_, b);
    double * pressure = stk::mesh::field_data(*pressure_, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int kNdim = k*nDim;
      double uSq = 0.0;
      for ( int j = 0; j < nDim; ++j ) {
        uSq += velocity[kNdim+j]*velocity[kNdim+j];
      } 
      pressure[k] = (gamma[k] - 1.0)*(totalEnergy[k] - 0.5*rho[k]*uSq);
    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_temperature ---------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::compute_temperature()
{
  // T = internal_energy / cv = (rhoE^t/rho - 1/2 uSq)/Cv
  const int nDim = realm_.meta_data().spatial_dimension();

  stk::mesh::Selector s_nodes = stk::mesh::selectField(*temperature_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    const double * totalEnergy = stk::mesh::field_data(*totalEnergy_, b);
    const double * rho = stk::mesh::field_data(*density_, b);
    const double * velocity = stk::mesh::field_data(*velocity_, b);
    const double * cv = stk::mesh::field_data(*cv_, b);
    double * temperature = stk::mesh::field_data(*temperature_, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int kNdim = k*nDim;
      double uSq = 0.0;
      for ( int j = 0; j < nDim; ++j ) {
        uSq += velocity[kNdim+j]*velocity[kNdim+j];
      } 
      const double internalE = totalEnergy[k]/rho[k] - 0.5*uSq;
      temperature[k] = internalE/cv[k];
    }
  }
}
  
//--------------------------------------------------------------------------
//-------- compute_speed_of_sound ------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::compute_speed_of_sound()
{
  // c = sqrt(gamma*P/rho)
  stk::mesh::Selector s_nodes = stk::mesh::selectField(*speedOfSound_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    const double * rho = stk::mesh::field_data(*density_, b);
    const double * gamma = stk::mesh::field_data(*gamma_, b);
    const double * pressure = stk::mesh::field_data(*pressure_, b);
    double * speedOfSound = stk::mesh::field_data(*speedOfSound_, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      speedOfSound[k] = std::sqrt(gamma[k]*pressure[k]/rho[k]);
    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_mach_number ---------------------------------------------
//--------------------------------------------------------------------------
void
GasDynamicsEquationSystem::compute_mach_number()
{
  // M = sqrt(uSq)/lilc
  const int nDim = realm_.meta_data().spatial_dimension();

  stk::mesh::Selector s_nodes = stk::mesh::selectField(*machNumber_);
  stk::mesh::BucketVector const& node_buckets =
    realm_.get_buckets( stk::topology::NODE_RANK, s_nodes );
  
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    const double * rho = stk::mesh::field_data(*density_, b);
    const double * gamma = stk::mesh::field_data(*gamma_, b);
    const double * pressure = stk::mesh::field_data(*pressure_, b);
    const double * velocity = stk::mesh::field_data(*velocity_, b);
    double * machNumber = stk::mesh::field_data(*machNumber_, b);
    
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int kNdim = k*nDim;
      double uSq = 0.0;
      for ( int j = 0; j < nDim; ++j ) {
        uSq += velocity[kNdim+j]*velocity[kNdim+j];
      } 
      const double lilc = std::sqrt(gamma[k]*pressure[k]/rho[k]);
      machNumber[k] = std::sqrt(uSq)/lilc;
    }
  }
}
  
} // namespace nalu
} // namespace Sierra
