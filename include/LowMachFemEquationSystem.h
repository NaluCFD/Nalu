/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef LowMachFemEquationSystem_h
#define LowMachFemEquationSystem_h

#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>

namespace stk{
struct topology;
}

namespace sierra{
namespace nalu{

class AlgorithmDriver;
class Realm;
class MomentumFemEquationSystem;
class ContinuityFemEquationSystem;
class LinearSystem;
class ProjectedNodalGradientEquationSystem;

/** Low-Mach formulation of the Navier-Stokes Equations (FEM)
 *
 *  This class is a thin-wrapper around sierra::nalu::ContinuityFemEquationSystem
 *  and sierra::nalu::MomentumFemEquationSystem that orchestrates the interactions
 *  between the velocity and the pressure Possion solves in the
 *  LowMachFemEquationSystem::solve_and_update method.
 */
class LowMachFemEquationSystem : public EquationSystem {

public:

  LowMachFemEquationSystem (
    EquationSystems& equationSystems);
  virtual ~LowMachFemEquationSystem();
  
  virtual void initialize();

  virtual void register_nodal_fields(
    stk::mesh::Part *part);

  virtual void register_element_fields(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

  virtual void register_open_bc(
    stk::mesh::Part *part,
    const stk::topology &partTopo,
    const OpenBoundaryConditionData &openBCData);

  virtual void pre_iter_work();
  virtual void solve_and_update();

  void copy_lagged();
  
  void project_nodal_velocity();

  MomentumFemEquationSystem *momentumEqSys_;
  ContinuityFemEquationSystem *continuityEqSys_;

  ScalarFieldType *density_;
  ScalarFieldType *viscosity_;

  VectorFieldType *dpdxL_;
  VectorFieldType *vrtmL_;

  bool isInit_;
};

/** Representation of the Momentum conservation equations in 2-D and 3-D
 *
 */
class MomentumFemEquationSystem : public EquationSystem {

public:

  MomentumFemEquationSystem(
    EquationSystems& equationSystems);
  virtual ~MomentumFemEquationSystem();

  virtual void initial_work();

  virtual void register_nodal_fields(
    stk::mesh::Part *part);

  virtual void register_interior_algorithm(
    stk::mesh::Part *part);

  virtual void register_inflow_bc(
    stk::mesh::Part *part,
    const stk::topology &partTopo,
    const InflowBoundaryConditionData &inflowBCData);

  virtual void register_open_bc(
    stk::mesh::Part *part,
    const stk::topology &partTopo,
    const OpenBoundaryConditionData &openBCData);

  virtual void register_wall_bc(
    stk::mesh::Part *part,
    const stk::topology &partTopo,
    const WallBoundaryConditionData &wallBCData);

  virtual void register_overset_bc();

  virtual void initialize();
  virtual void reinitialize_linear_system();
  
  virtual void predict_state();

  virtual void register_initial_condition_fcn(
    stk::mesh::Part *part,
    const std::map<std::string, std::string> &theNames,
    const std::map<std::string, std::vector<double> > &theParams);
  
  VectorFieldType *velocity_;
  VectorFieldType *uTmp_;
  ScalarFieldType *density_;
  ScalarFieldType *viscosity_;
  AlgorithmDriver *cflReyAlgDriver_;

  // saved of mesh parts that are not to be projected
  std::vector<stk::mesh::Part *> notProjectedPart_;
};

class ContinuityFemEquationSystem : public EquationSystem {

public:

  ContinuityFemEquationSystem(
    EquationSystems& equationSystems);
  virtual ~ContinuityFemEquationSystem();

  virtual void register_nodal_fields(
    stk::mesh::Part *part);

  virtual void register_interior_algorithm(
    stk::mesh::Part *part);

  virtual void register_overset_bc();

  virtual void initialize();
  virtual void reinitialize_linear_system();    
  
  virtual void register_initial_condition_fcn(
      stk::mesh::Part *part,
      const std::map<std::string, std::string> &theNames,
      const std::map<std::string, std::vector<double> > &theParams);

  virtual void register_inflow_bc(
    stk::mesh::Part *part,
    const stk::topology &partTopo,
    const InflowBoundaryConditionData &inflowBCData);

  virtual void register_open_bc(
    stk::mesh::Part *part,
    const stk::topology &partTopo,
    const OpenBoundaryConditionData &openBCData);

  virtual void manage_projected_nodal_gradient(
    EquationSystems& eqSystems);
  virtual void compute_projected_nodal_gradient();
  
  const bool managePNG_;
  ScalarFieldType *pressure_;
  VectorFieldType *dpdx_;
  ScalarFieldType *density_;

  ScalarFieldType *pTmp_;

  ProjectedNodalGradientEquationSystem *projectedNodalGradEqs_;
};

} // namespace nalu
} // namespace Sierra

#endif
