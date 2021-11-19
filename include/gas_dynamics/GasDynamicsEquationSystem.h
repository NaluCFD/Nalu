/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef GasDynamicsEquationSystem_h
#define GasDynamicsEquationSystem_h

#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>

namespace stk{
struct topology;
}

namespace sierra{
namespace nalu{

class AlgorithmDriver;
class AssembleGasDynamicsAlgorithmDriver;
class Realm;
class EquationSystems;

class GasDynamicsEquationSystem : public EquationSystem {

public:

  GasDynamicsEquationSystem(
    EquationSystems& equationSystems,
    bool debugOutput);
  virtual ~GasDynamicsEquationSystem();
  
  virtual void initial_work();

  virtual void register_nodal_fields(
    stk::mesh::Part *part);

  virtual void register_edge_fields(
    stk::mesh::Part *part);

  virtual void register_element_fields(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

  void register_interior_algorithm(
    stk::mesh::Part *part);
  
  void register_inflow_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const InflowBoundaryConditionData &inflowBCData);
  
  void register_open_bc(
    stk::mesh::Part *part,
    const stk::topology &partTopo,
    const OpenBoundaryConditionData &openBCData);

  void register_wall_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const WallBoundaryConditionData &wallBCData);

  virtual void register_symmetry_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const SymmetryBoundaryConditionData &symmetryBCData);

  virtual void register_overset_bc();

  void initialize();
  void reinitialize_linear_system();
  
  void predict_state();
  
  void solve_and_update();

  void assemble_gas_dynamics();
  void update_gas_dynamics();  
  void dump_state(const std::string indicator);

  // keep it real
  double provide_scaled_norm() {return 1.0;}
  double provide_norm() {return 1.0;}
  double provide_norm_increment() {return 1.0;}

  void compute_density();
  void compute_momentum();
  void compute_static_enthalpy();
  void compute_total_enthalpy();
  void compute_total_energy();

  void compute_velocity();
  void compute_pressure();
  void compute_temperature();
  void compute_speed_of_sound();
  void compute_mach_number();

  ScalarFieldType *density_;
  VectorFieldType *momentum_;
  ScalarFieldType *totalEnergy_;
  VectorFieldType *velocity_;
  ScalarFieldType *totalEnthalpy_;
  ScalarFieldType *staticEnthalpy_;
  ScalarFieldType *pressure_;
  ScalarFieldType *temperature_;
  ScalarFieldType *machNumber_;
  ScalarFieldType *speedOfSound_;
  ScalarFieldType *cp_;
  ScalarFieldType *cv_;
  ScalarFieldType *viscosity_;
  ScalarFieldType *thermalCond_;
  ScalarFieldType *gamma_;
  ScalarFieldType *dualNodalVolume_;
  GenericFieldType *rhsGasDyn_;

  AssembleGasDynamicsAlgorithmDriver *assembleGasDynAlgDriver_;
  AlgorithmDriver *cflReyAlgDriver_;

  bool isInit_;
  const bool debugOutput_;

  // boundary condition mapping
  std::vector<Algorithm *> gasDynBcDataMapAlg_;
  
  // density and enthalpy algorithm
  std::vector<Algorithm *> densityAlg_;
  std::vector<Algorithm *> enthalpyAlg_;
  
  // saved of mesh parts that are not to be touched
  std::vector<stk::mesh::Part *> dirichletPart_;
};


} // namespace nalu
} // namespace Sierra

#endif
