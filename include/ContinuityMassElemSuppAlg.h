/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ContinuityMassElemSuppAlg_h
#define ContinuityMassElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class ContinuityMassElemSuppAlg : public SupplementalAlgorithm
{
public:

  ContinuityMassElemSuppAlg(
    Realm &realm);

  virtual ~ContinuityMassElemSuppAlg() {}

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

  ScalarFieldType *densityNm1_;
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  VectorFieldType *coordinates_;

  double dt_;
  double gamma1_;
  double gamma2_;
  double gamma3_;
  const int nDim_;
  const bool useShifted_;

  // scratch space
  std::vector<double> ws_shape_function_;
  std::vector<double> ws_qm1_;
  std::vector<double> ws_qN_;
  std::vector<double> ws_qNp1_;
  std::vector<double> ws_rhoNm1_;
  std::vector<double> ws_rhoN_;
  std::vector<double> ws_rhoNp1_;
  std::vector<double> ws_coordinates_;
  std::vector<double> ws_scv_volume_;
};

} // namespace nalu
} // namespace Sierra

#endif
