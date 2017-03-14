/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ScalarDiffFemSuppAlg_h
#define ScalarDiffFemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <AlgTraits.h>
#include <ElemDataRequests.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;
class Hex8FEM;

template<class AlgTraits>
class ScalarDiffFemSuppAlg : public SupplementalAlgorithm
{
public:
  ScalarDiffFemSuppAlg(
    Realm &realm,
    ScalarFieldType *temperature,
    ScalarFieldType *thermalCond);

  virtual ~ScalarDiffFemSuppAlg();

  virtual void element_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element,
    ScratchViews &scratchViews);
  
  const stk::mesh::BulkData *bulkData_;

  ScalarFieldType *temperature_;
  ScalarFieldType *thermalCond_;
  VectorFieldType *coordinates_;

  // master element
  Hex8FEM * meFEM_;
  double *ipWeight_;

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
