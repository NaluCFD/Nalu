/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MixtureFractionFemEquationSystem_h
#define MixtureFractionFemEquationSystem_h

#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>

namespace stk{
struct topology;
}

namespace sierra{
namespace nalu{

class AlgorithmDriver;
class EquationSystems;
class LinearSystem;
class ProjectedNodalGradientEquationSystem;
class Realm;

class MixtureFractionFemEquationSystem : public EquationSystem {

public:

  MixtureFractionFemEquationSystem(
    EquationSystems& equationSystems,
    const bool outputClippingDiag,
    const double deltaZClip,
    const bool computePng);
  virtual ~MixtureFractionFemEquationSystem();

  void populate_derived_quantities();
  
  virtual void register_nodal_fields(
    stk::mesh::Part *part);

  void register_element_fields(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

  void register_interior_algorithm(
    stk::mesh::Part *part);

  virtual void register_inflow_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const InflowBoundaryConditionData &inflowBCData);

  virtual void register_open_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const OpenBoundaryConditionData &OpenBCData);

  virtual void register_wall_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const WallBoundaryConditionData &symmetryBCData);

  virtual void register_symmetry_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const SymmetryBoundaryConditionData &symmetryBCData);

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

  const bool outputClippingDiag_;
  const double deltaZClip_;
  const bool computePNG_;
  const bool managePNG_;

  ScalarFieldType *mixFrac_;
  ScalarFieldType *mixFracUF_;
  ScalarFieldType *zTmp_;
  VectorFieldType *Gjz_;
  VectorFieldType *velocity_;
  ScalarFieldType *density_;
  ScalarFieldType *visc_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *evisc_;
  AlgorithmDriver *diffFluxCoeffAlgDriver_; 
  AlgorithmDriver *cflReyAlgDriver_; 

  ProjectedNodalGradientEquationSystem *projectedNodalGradEqs_;
  
  bool isInit_;
};


} // namespace nalu
} // namespace Sierra

#endif
