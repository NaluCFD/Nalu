/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef EnthalpyEquationSystem_h
#define EnthalpyEquationSystem_h

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
class AssembleWallHeatTransferAlgorithmDriver;
class LinearSystem;
class EquationSystems;
class TemperaturePropAlgorithm;

class EnthalpyEquationSystem : public EquationSystem {

public:

  EnthalpyEquationSystem(
    EquationSystems& equationSystems,
    const double minT,
    const double maxT);
  virtual ~EnthalpyEquationSystem();
  
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
    const stk::topology &theTopo,
    const OpenBoundaryConditionData &openBCData);

  void register_wall_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const WallBoundaryConditionData &wallBCData);

  void register_contact_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const ContactBoundaryConditionData &contactBCData);

  virtual void register_symmetry_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const SymmetryBoundaryConditionData &symmetryBCData);

  virtual void register_non_conformal_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

  virtual void register_overset_bc();

  void initialize();
  void reinitialize_linear_system();
  
  void predict_state();
  
  void solve_and_update();
  void post_iter_work();
  void post_adapt_work();
  void extract_temperature();
  void post_converged_work();
  void initial_work();
  
  void temperature_bc_setup(
    std::vector<double> userSpecData,
    stk::mesh::Part *part,
    ScalarFieldType *temperatureBc,
    ScalarFieldType *enthalpyBc,
    const bool isInterface = false,
    const bool copyBcVal = true);
  
  const double minimumT_;
  const double maximumT_;
  
  ScalarFieldType *enthalpy_;
  ScalarFieldType *temperature_;
  VectorFieldType *dhdx_;
  ScalarFieldType *hTmp_;
  ScalarFieldType *visc_;
  ScalarFieldType *tvisc_;
  ScalarFieldType *evisc_;
  ScalarFieldType *thermalCond_;
  ScalarFieldType *specHeat_;
  ScalarFieldType *divQ_;
  ScalarFieldType *pOld_;
  
  AssembleNodalGradAlgorithmDriver *assembleNodalGradAlgDriver_;
  AlgorithmDriver *diffFluxCoeffAlgDriver_;
  AssembleWallHeatTransferAlgorithmDriver *assembleWallHeatTransferAlgDriver_;
  
  bool pmrCouplingActive_;
  bool lowSpeedCompressActive_;
  bool isInit_;

  std::vector<TemperaturePropAlgorithm *> enthalpyFromTemperatureAlg_;
  std::vector<Algorithm *> bdf2CopyStateAlg_;

  // bc enthalpy
  std::vector<TemperaturePropAlgorithm *> bcEnthalpyFromTemperatureAlg_;
  std::vector<Algorithm *> bcCopyStateAlg_;
  
};


} // namespace nalu
} // namespace Sierra

#endif
