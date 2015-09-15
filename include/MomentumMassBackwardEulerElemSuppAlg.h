/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumMassBackwardEulerElemSuppAlg_h
#define MomentumMassBackwardEulerElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class MomentumMassBackwardEulerElemSuppAlg : public SupplementalAlgorithm
{
public:

  MomentumMassBackwardEulerElemSuppAlg(
    Realm &realm);

  virtual ~MomentumMassBackwardEulerElemSuppAlg() {}

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

  VectorFieldType *velocityN_;
  VectorFieldType *velocityNp1_;
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  VectorFieldType *Gjp_;
  VectorFieldType *coordinates_;

  double dt_;
  const int nDim_;
  const bool useShifted_;

  // scratch space
  std::vector<double> uNScv_;
  std::vector<double> uNp1Scv_;
  std::vector<double> GjpScv_;

  std::vector<double> ws_shape_function_;
  std::vector<double> ws_uN_;
  std::vector<double> ws_uNp1_;
  std::vector<double> ws_Gjp_;
  std::vector<double> ws_rhoN_;
  std::vector<double> ws_rhoNp1_;
  std::vector<double> ws_coordinates_;
  std::vector<double> ws_scv_volume_;
};

} // namespace nalu
} // namespace Sierra

#endif
