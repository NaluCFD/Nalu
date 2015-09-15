/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef MomentumBuoyancySrcElemSuppAlg_h
#define MomentumBuoyancySrcElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class MomentumBuoyancySrcElemSuppAlg : public SupplementalAlgorithm
{
public:

  MomentumBuoyancySrcElemSuppAlg(
    Realm &realm);

  virtual ~MomentumBuoyancySrcElemSuppAlg() {}

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

  ScalarFieldType *densityNp1_;
  VectorFieldType *coordinates_;

  const int nDim_;
  double rhoRef_;
  const bool useShifted_;
  std::vector<double> gravity_;

  // scratch space
  std::vector<double> ws_shape_function_;
  std::vector<double> ws_rhoNp1_;
  std::vector<double> ws_coordinates_;
  std::vector<double> ws_scv_volume_;
};

} // namespace nalu
} // namespace Sierra

#endif
