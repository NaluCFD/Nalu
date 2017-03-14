/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef RadTransFemElemSuppAlg_h
#define RadTransFemElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class RadiativeTransportEquationSystem;
class Realm;
class MasterElement;
class Hex8FEM;

class RadTransFemElemSuppAlg : public SupplementalAlgorithm
{
public:

  RadTransFemElemSuppAlg(
    Realm &realm,
    RadiativeTransportEquationSystem *radEqSystem);

  virtual ~RadTransFemElemSuppAlg();

  virtual void setup();

  virtual void elem_resize(
    MasterElement *meSCS,
    MasterElement *meSCV);

  virtual void elem_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    MasterElement *meSCS,
    MasterElement *meSCV);

  const RadiativeTransportEquationSystem *radEqSystem_;
  
  const stk::mesh::BulkData *bulkData_;

  ScalarFieldType *intensity_;
  ScalarFieldType *absorption_;
  ScalarFieldType *scattering_;
  ScalarFieldType *scalarFlux_;
  ScalarFieldType *radiationSource_;
  VectorFieldType *coordinates_;

  // master element
  Hex8FEM * meFEM_;
  double *ipWeight_;
  const double invPi_;
  const int nDim_;
  const double rowSumLump_;
  const double consistentMass_;
  const bool linearNorm_;
  const bool useUpper_;

  // fixed space
  std::vector<double> ws_Sk_;

  // scratch space; geometry
  std::vector<double> ws_dndx_;
  std::vector<double> ws_deriv_;
  std::vector<double> ws_det_j_;
  std::vector<double> ws_shape_function_;
  std::vector<double> ws_gUpper_;
  std::vector<double> ws_gLower_;

  // scratch space; fields
  std::vector<double> ws_intensity_;
  std::vector<double> ws_absorption_;
  std::vector<double> ws_scattering_;
  std::vector<double> ws_scalarFlux_;
  std::vector<double> ws_radiationSource_;
  std::vector<double> ws_coordinates_;
};

} // namespace nalu
} // namespace Sierra

#endif
