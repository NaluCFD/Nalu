/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumAdvDiffElemSuppAlg_h
#define MomentumAdvDiffElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class MomentumAdvDiffElemSuppAlg : public SupplementalAlgorithm
{
public:

  MomentumAdvDiffElemSuppAlg(
    Realm &realm,
    VectorFieldType *velocity,
    ScalarFieldType *viscosity);

  virtual ~MomentumAdvDiffElemSuppAlg() {}

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
  
  const stk::mesh::BulkData *bulkData_;

  VectorFieldType *velocityNp1_;
  VectorFieldType *coordinates_;
  ScalarFieldType *viscosity_;
  GenericFieldType *massFlowRate_;

  const int nDim_;
  const double includeDivU_;

  // fixed space
  std::vector<double> ws_uIp_;

  // scratch space; geometry
  std::vector<double> ws_scs_areav_;
  std::vector<double> ws_dndx_;
  std::vector<double> ws_deriv_;
  std::vector<double> ws_det_j_;
  std::vector<double> ws_shape_function_;

  // scratch space; fields
  std::vector<double> ws_uNp1_;
  std::vector<double> ws_coordinates_;
  std::vector<double> ws_viscosity_;
};

} // namespace nalu
} // namespace Sierra

#endif
