/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MeshDisplacementElasticElemKernel_H
#define MeshDisplacementElasticElemKernel_H

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
class MeshDisplacementElasticElemKernel: public Kernel
{
public:
  MeshDisplacementElasticElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    VectorFieldType*,
    ElemDataRequests&);

  virtual ~MeshDisplacementElasticElemKernel();

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
  MeshDisplacementElasticElemKernel() = delete;

  VectorFieldType *meshDisplacement_{nullptr};

  VectorFieldType *coordinates_{nullptr};
  ScalarFieldType *dualNodalVolume_{nullptr};

  ScalarFieldType *mu_{nullptr};
  ScalarFieldType *lambda_{nullptr};

  /// Integration point to node maping
  const int* lrscv_;

  /// Shape functions
  AlignedViewType<DoubleType[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_ {"view_shape_func"};
};

}  // nalu
}  // sierra

#endif /* MeshDisplacementElasticElemKernel_H */
