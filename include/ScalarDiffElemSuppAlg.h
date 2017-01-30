/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ScalarDiffElemSuppAlg_h
#define ScalarDiffElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

class ScalarDiffElemSuppAlg : public SupplementalAlgorithm
{
public:

  ScalarDiffElemSuppAlg(
    Realm &realm,
    ScalarFieldType *scalarQ,
    ScalarFieldType *diffFluxCoeff,
    const stk::topology &theTopo);

  virtual ~ScalarDiffElemSuppAlg() {}

  virtual void setup();

  virtual void elem_resize() {}

  virtual void element_execute(
    double *lhs,
    double *rhs,
    stk::mesh::Entity element);
  
  const stk::mesh::BulkData *bulkData_;

  ScalarFieldType *scalarQ_;
  ScalarFieldType *diffFluxCoeff_;
  VectorFieldType *coordinates_;

  // master element
  MasterElement  *meSCS_;
  const int *lrscv_;
  const int nodesPerElement_;
  const int numScsIp_;
  const int nDim_;

  // scratch space; geometry
  std::vector<double> ws_scs_areav_;
  std::vector<double> ws_dndx_;
  std::vector<double> ws_deriv_;
  std::vector<double> ws_det_j_;
  std::vector<double> ws_shape_function_;

  // scratch space; fields
  std::vector<double> ws_scalarQ_;
  std::vector<double> ws_diffFluxCoeff_;
  std::vector<double> ws_coordinates_;
};

} // namespace nalu
} // namespace Sierra

#endif
