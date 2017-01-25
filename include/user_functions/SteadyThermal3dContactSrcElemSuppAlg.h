/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef SteadyThermal3dContactSrcElemSuppAlg_h
#define SteadyThermal3dContactSrcElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class SteadyThermal3dContactSrcElemSuppAlg : public SupplementalAlgorithm
{
public:

  SteadyThermal3dContactSrcElemSuppAlg(
    Realm &realm);

  virtual ~SteadyThermal3dContactSrcElemSuppAlg() {}

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

  VectorFieldType *coordinates_;

  const double a_;
  const double k_;
  const double pi_;
  const bool useShifted_;
  const int nDim_;
  const bool evalAtIps_;

  // scratch space
  std::vector<double> scvCoords_;
  std::vector<double> ws_shape_function_;
  std::vector<double> ws_coordinates_;
  std::vector<double> ws_scv_volume_;
  std::vector<double> ws_nodalSrc_;
};

} // namespace nalu
} // namespace Sierra

#endif
