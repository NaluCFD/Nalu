/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef HeatCondFemElemSuppAlg_h
#define HeatCondFemElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;
class Hex8FEM;

class HeatCondFemElemSuppAlg : public SupplementalAlgorithm
{
public:

  HeatCondFemElemSuppAlg(
    Realm &realm);

  virtual ~HeatCondFemElemSuppAlg() {}

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

  ScalarFieldType *temperature_;
  ScalarFieldType *thermalCond_;
  VectorFieldType *coordinates_;

  // master element
  Hex8FEM * meFEM_;
  double *ipWeight_;
  const int nDim_;

  // scratch space; geometry
  std::vector<double> ws_dndx_;
  std::vector<double> ws_deriv_;
  std::vector<double> ws_det_j_;
  std::vector<double> ws_shape_function_;

  // scratch space; fields
  std::vector<double> ws_temperature_;
  std::vector<double> ws_thermalCond_;
  std::vector<double> ws_coordinates_;
};

} // namespace nalu
} // namespace Sierra

#endif
