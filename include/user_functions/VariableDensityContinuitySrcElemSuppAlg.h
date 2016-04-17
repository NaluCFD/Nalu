/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef VariableDensityContinuitySrcElemSuppAlg_h
#define VariableDensityContinuitySrcElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class VariableDensityContinuitySrcElemSuppAlg : public SupplementalAlgorithm
{
public:

  VariableDensityContinuitySrcElemSuppAlg(
    Realm &realm);

  virtual ~VariableDensityContinuitySrcElemSuppAlg() {}

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

  const int nDim_;
  const double unot_;
  const double vnot_;
  const double wnot_;
  const double znot_;
  const double rhoP_;
  const double rhoS_;
  const double a_;
  const double amf_;
  const double pi_;
  double projTimeScale_;
  const bool useShifted_;

  // scratch space (at constructor)
  std::vector<double> scvCoords_;
  // at elem_resize
  std::vector<double> ws_shape_function_;
  std::vector<double> ws_coordinates_;
  std::vector<double> ws_scv_volume_;
};

} // namespace nalu
} // namespace Sierra

#endif
