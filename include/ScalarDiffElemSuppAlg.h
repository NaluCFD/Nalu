/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ScalarDiffElemSuppAlg_h
#define ScalarDiffElemSuppAlg_h

#include <SupplementalAlgorithm.h>
#include <AlgTraits.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_topology/topology.hpp>

// Kokkos
#include <Kokkos_Core.hpp>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

template<class AlgTraits>
class ScalarDiffElemSuppAlg : public SupplementalAlgorithm
{
public:

  ScalarDiffElemSuppAlg(
    Realm &realm,
    ScalarFieldType *scalarQ,
    ScalarFieldType *diffFluxCoeff,
    const stk::topology &theTopo);

  virtual ~ScalarDiffElemSuppAlg() {}

  void element_execute(
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

  // scratch space; geometry
  Kokkos::View<double**> ws_scs_areav_;
  Kokkos::View<double***> ws_dndx_;
  Kokkos::View<double*> ws_deriv_;
  Kokkos::View<double*> ws_det_j_;
  Kokkos::View<double**> ws_shape_function_;

  // scratch space; fields
  Kokkos::View<double*> ws_scalarQ_;
  Kokkos::View<double*> ws_diffFluxCoeff_;
  Kokkos::View<double**> ws_coordinates_;
};

} // namespace nalu
} // namespace Sierra

#endif
