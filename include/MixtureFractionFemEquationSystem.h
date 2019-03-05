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
class Realm;
class LinearSystem;
class EquationSystems;

class MixtureFractionFemEquationSystem : public EquationSystem {

public:

  MixtureFractionFemEquationSystem(
    EquationSystems& equationSystems,
    const bool outputClippingDiag,
    const double deltaZClip);
  virtual ~MixtureFractionFemEquationSystem();

  void populate_derived_quantities();
  
  virtual void register_nodal_fields(
    stk::mesh::Part *part);

  void register_element_fields(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

  void register_interior_algorithm(
    stk::mesh::Part *part);

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

  const bool outputClippingDiag_;
  const double deltaZClip_;

  ScalarFieldType *mixFrac_;
  ScalarFieldType *mixFracUF_;
  ScalarFieldType *zTmp_;
  VectorFieldType *velocity_;
  ScalarFieldType *density_;
  ScalarFieldType *visc_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *evisc_;
  AlgorithmDriver *diffFluxCoeffAlgDriver_; 
  AlgorithmDriver *cflReyAlgDriver_; 
  bool isInit_;
};


} // namespace nalu
} // namespace Sierra

#endif
