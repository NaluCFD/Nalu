/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef LowMachMonolithicEquationSystem_h
#define LowMachMonolithicEquationSystem_h

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
class AssembleNodalGradUAlgorithmDriver;
class LinearSystem;

/** Low-Mach formulation of the Navier-Stokes Equations
 *
 *  This class supports the monolithic uvwp system.
 */
class LowMachMonolithicEquationSystem : public EquationSystem {

public:

  LowMachMonolithicEquationSystem (
    EquationSystems& equationSystems);
  virtual ~LowMachMonolithicEquationSystem();
  
  virtual void initialize();

  virtual void initial_work();
  
  virtual void register_nodal_fields(
    stk::mesh::Part *part);
 
  virtual void register_element_fields(
    stk::mesh::Part *part,
    const stk::topology &theTopo);
  
  virtual void register_edge_fields(
    stk::mesh::Part *part);

  virtual void register_interior_algorithm(
    stk::mesh::Part *part);
  
  virtual void register_inflow_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const InflowBoundaryConditionData &inflowBCData);

  virtual void register_open_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const OpenBoundaryConditionData &openBCData);

  virtual void register_wall_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const WallBoundaryConditionData &wallBCData);
    
  virtual void register_symmetry_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const SymmetryBoundaryConditionData &symmetryBCData);

  virtual void register_initial_condition_fcn(
      stk::mesh::Part *part,
      const std::map<std::string, std::string> &theNames,
      const std::map<std::string, std::vector<double> > &theParams);

  virtual void pre_timestep_work();
  virtual void pre_iter_work();
  virtual void solve_and_update();
  virtual void predict_state();
  virtual void post_converged_work();

  void copy_uvwp_to_uvw_p();
  void copy_uvw_p_to_uvwp();

  VectorFieldType *velocity_;
  ScalarFieldType *pressure_;
  GenericFieldType *dudx_;
  VectorFieldType *dpdx_;
  VectorFieldType *dpdxOld_;
  GenericFieldType *uvwp_;
  GenericFieldType *uvwpTmp_;
  VectorFieldType *coordinates_;
  ScalarFieldType *density_;
  ScalarFieldType *viscosity_;

  AssembleNodalGradUAlgorithmDriver *assembleNodalGradUAlgDriver_;
  AssembleNodalGradAlgorithmDriver *assembleNodalGradPAlgDriver_;
  AlgorithmDriver *cflReyAlgDriver_;

  const int nDim_;
  const int sizeOfSystem_;
  bool edgeNodalGradientU_;
  bool edgeNodalGradientP_;
  bool isInit_;    
};
  
} // namespace nalu
} // namespace Sierra

#endif
