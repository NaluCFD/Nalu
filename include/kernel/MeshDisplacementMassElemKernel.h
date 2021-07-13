/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MeshDisplacementMassElemKernel_H
#define MeshDisplacementMassElemKernel_H

#include "kernel/Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class TimeIntegrator;
class SolutionOptions;
class MasterElement;
class ElemDataRequests;

template<typename AlgTraits>
class MeshDisplacementMassElemKernel: public Kernel
{
public:
  MeshDisplacementMassElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    VectorFieldType*,
    ElemDataRequests&,
    const bool);

  virtual ~MeshDisplacementMassElemKernel();

  /** Perform pre-timestep work for the computational kernel
   */
  virtual void setup(const TimeIntegrator&);

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  MeshDisplacementMassElemKernel() = delete;

  VectorFieldType *meshDisplacement_{nullptr};
  VectorFieldType *meshDisplacementNp1_{nullptr};
  VectorFieldType *meshDisplacementNm1_{nullptr};

  VectorFieldType *coordinates_{nullptr};
  ScalarFieldType *dualNodalVolume_{nullptr};
  ScalarFieldType *density_{nullptr};

  double dt_{0.0};
  double gamma_1_{0.0};
  double gamma_2_{0.0};
  double gamma_3_{0.0};

  /// Integration point to node mapping
  const int* ipNodeMap_;

  const bool lumpedMass_;

  /// Shape functions
  AlignedViewType<DoubleType[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_ {"view_shape_func"};
};

}  // nalu
}  // sierra

#endif /* MeshDisplacementMassElemKernel_H */
