/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef VolumeOfFluidEquationSystem_h
#define VolumeOfFluidEquationSystem_h

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
class AssembleNodalGradAlgorithmDriver;
class LinearSystem;
class EquationSystems;
class ProjectedNodalGradientEquationSystem;

class VolumeOfFluidEquationSystem : public EquationSystem {

public:

  VolumeOfFluidEquationSystem(
    EquationSystems& equationSystems,
    const bool outputClippingDiag,
    const double deltaZClip);
  virtual ~VolumeOfFluidEquationSystem();

  void populate_derived_quantities();
  
  virtual void register_nodal_fields(
    stk::mesh::Part *part);

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

  virtual void register_initial_condition_fcn(
      stk::mesh::Part *part,
      const std::map<std::string, std::string> &theNames,
      const std::map<std::string, std::vector<double> > &theParams);

  void initialize();
  void reinitialize_linear_system();
  
  void predict_state();
  
  void solve_and_update();
  void update_and_clip();

  void manage_projected_nodal_gradient(
    EquationSystems& eqSystems);
  void compute_projected_nodal_gradient();
  
  void sharpen_interface_explicit();
  void smooth_vof();
  void compute_interface_normal();

  const bool managePNG_;
  const bool outputClippingDiag_;
  const double deltaVofClip_;

  ScalarFieldType *vof_;
  ScalarFieldType *vofSmoothed_;
  ScalarFieldType *interfaceNormal_;
  VectorFieldType *dvofdx_;
  ScalarFieldType *vofTmp_;
  
  AssembleNodalGradAlgorithmDriver *assembleNodalGradAlgDriver_;
  
  ProjectedNodalGradientEquationSystem *projectedNodalGradEqs_;
  
  bool isInit_;
};


} // namespace nalu
} // namespace Sierra

#endif
