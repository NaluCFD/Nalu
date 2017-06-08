/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SCALARDIFFFEMKERNEL_H
#define SCALARDIFFFEMKERNEL_H

#include "Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>
#include <vector>

namespace sierra {
namespace nalu {

class MasterElement;
class Hex8FEM;
class ElemDataRequests;
class ScratchViews;

/** CVFEM scalar advection/diffusion kernel
 */
template<typename AlgTraits>
class ScalarDiffFemKernel: public Kernel
{
public:
  ScalarDiffFemKernel(
    const stk::mesh::BulkData&,
    ScalarFieldType*,
    ScalarFieldType*,
    bool,
    ElemDataRequests&);

  virtual ~ScalarDiffFemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<double**>&,
    SharedMemView<double*>&,
    ScratchViews&);

private:
  ScalarDiffFemKernel() = delete;

  const stk::mesh::BulkData* bulkData_;
  ScalarFieldType *scalarQ_{nullptr};
  ScalarFieldType *diffFluxCoeff_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  // master element
  Hex8FEM * meFEM_;
  double *ipWeight_;

  /// Shape functions
  Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "v_shape_func" };
};

}  // nalu
}  // sierra

#endif /* SCALARDIFFFEMKERNEL_H */
