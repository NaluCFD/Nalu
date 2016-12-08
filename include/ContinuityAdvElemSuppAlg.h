/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ContinuityAdvElemSuppAlg_h
#define ContinuityAdvElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class ContinuityAdvElemSuppAlg : public SupplementalAlgorithm
{
public:

  ContinuityAdvElemSuppAlg(
    Realm &realm);

  virtual ~ContinuityAdvElemSuppAlg() {}

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


  // extract fields; nodal
  VectorFieldType *velocityRTM_;
  VectorFieldType *Gpdx_;
  ScalarFieldType *pressure_;
  ScalarFieldType *densityNp1_;
  VectorFieldType *coordinates_;

  double projTimeScale_;
  const int nDim_;

  const bool meshMotion_;
  const bool shiftMdot_;
  const bool shiftPoisson_;
  const bool reducedSensitivities_;
  const double interpTogether_;
  const double om_interpTogether_;

  // fixed size
  std::vector<double> ws_uIp_;
  std::vector<double> ws_rho_uIp_;
  std::vector<double> ws_Gpdx_Ip_;
  std::vector<double> ws_dpdxIp_;

  // scratch space
  std::vector<double> ws_velocityRTM_;
  std::vector<double> ws_Gpdx_;
  std::vector<double> ws_pressure_;
  std::vector<double> ws_densityNp1_;
  std::vector<double> ws_coordinates_;
  std::vector<double> ws_scs_areav_;
  std::vector<double> ws_dndx_;
  std::vector<double> ws_dndx_lhs_;
  std::vector<double> ws_deriv_;
  std::vector<double> ws_det_j_;
  std::vector<double> ws_shape_function_;
};

} // namespace nalu
} // namespace Sierra

#endif
