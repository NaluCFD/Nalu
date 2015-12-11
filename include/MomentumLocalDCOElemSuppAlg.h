/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumLocalDCOElemSuppAlg_h
#define MomentumLocalDCOElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class MomentumLocalDCOElemSuppAlg : public SupplementalAlgorithm
{
public:

  MomentumLocalDCOElemSuppAlg(
    Realm &realm,
    VectorFieldType *velocity,
    GenericFieldType *Gju,
    ScalarFieldType *diffFluxCoeff,
    const double fourthFac);

  virtual ~MomentumLocalDCOElemSuppAlg() {}

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

  VectorFieldType *velocityNm1_;
  VectorFieldType *velocityN_;
  VectorFieldType *velocityNp1_;
  ScalarFieldType *densityNm1_;
  ScalarFieldType *densityN_;
  ScalarFieldType *densityNp1_;
  ScalarFieldType *pressure_;
  VectorFieldType *velocityRTM_;
  VectorFieldType *coordinates_;
  ScalarFieldType *diffFluxCoeff_;
  GenericFieldType *Gju_;
  double dt_;
  const int nDim_;
  double gamma1_;
  double gamma2_;
  double gamma3_;
  const double Cupw_;
  const double small_;
  const double fourthFac_;

  // fixed space
  std::vector<double> ws_dukdxScs_;
  std::vector<double> ws_rhovScs_;
  std::vector<double> ws_vrtmScs_;
  std::vector<double> ws_dpdxScs_;
  
  // scratch space; geometry
  std::vector<double> ws_scs_areav_;
  std::vector<double> ws_dndx_;
  std::vector<double> ws_deriv_;
  std::vector<double> ws_det_j_;
  std::vector<double> ws_shape_function_;
  std::vector<double> ws_gUpper_;
  std::vector<double> ws_gLower_;

  // scratch space; fields
  std::vector<double> ws_uNm1_;
  std::vector<double> ws_uN_;
  std::vector<double> ws_uNp1_;
  std::vector<double> ws_rhoNm1_;
  std::vector<double> ws_rhoN_;
  std::vector<double> ws_rhoNp1_;
  std::vector<double> ws_pressure_;
  std::vector<double> ws_velocityRTM_;
  std::vector<double> ws_coordinates_;
  //std::vector<double> ws_diffFluxCoeff_;
  //std::vector<double> ws_Gju_;
};

} // namespace nalu
} // namespace Sierra

#endif
