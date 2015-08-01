/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MeshDisplacementEquationSystem_h
#define MeshDisplacementEquationSystem_h

#include <EquationSystem.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>

namespace stk{
struct topology;
}

namespace sierra{
namespace nalu{

class AlgorithmDriver;
class AssembleNodalGradUAlgorithmDriver;
class Realm;
class LinearSystem;

class MeshDisplacementEquationSystem : public EquationSystem {

public:

  MeshDisplacementEquationSystem(
    EquationSystems& equationSystems,
    const bool activateMass,
    const bool deformWrtModelCoords);
  virtual ~MeshDisplacementEquationSystem();

  void initial_work();

  void register_nodal_fields(
    stk::mesh::Part *part);

  void register_element_fields(
    stk::mesh::Part *part,
    const stk::topology &theTopo);

  void register_interior_algorithm(
    stk::mesh::Part *part);

  void register_wall_bc(
    stk::mesh::Part *part,
    const stk::topology &theTopo,
    const WallBoundaryConditionData &wallBCData);

  void register_overset_bc();

  void initialize();
  void reinitialize_linear_system();
  
  void predict_state();
  void solve_and_update();
  void compute_current_coordinates();
  void compute_div_mesh_velocity();
  
  const bool activateMass_;
  const bool deformWrtModelCoords_;
  bool isInit_;
  VectorFieldType *meshDisplacement_;
  VectorFieldType *meshVelocity_;
  GenericFieldType *dvdx_;
  ScalarFieldType *divV_;
  VectorFieldType *coordinates_;
  VectorFieldType *currentCoordinates_;
  ScalarFieldType *dualNodalVolume_;
  ScalarFieldType *density_;
  ScalarFieldType *lameMu_;
  ScalarFieldType *lameLambda_;
  VectorFieldType *dxTmp_;

  AssembleNodalGradUAlgorithmDriver *assembleNodalGradAlgDriver_;

};

} // namespace nalu
} // namespace Sierra

#endif
