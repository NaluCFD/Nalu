/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef VolumeOfFluidSucvNsoElemKernel_H
#define VolumeOfFluidSucvNsoElemKernel_H

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

/** CMM (BDF2/BE) for scalar equation
 */
template<typename AlgTraits>
class VolumeOfFluidSucvNsoElemKernel: public Kernel
{
public:
  VolumeOfFluidSucvNsoElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ScalarFieldType*,
    const double,
    const double,
    ElemDataRequests&);

  virtual ~VolumeOfFluidSucvNsoElemKernel();

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
  VolumeOfFluidSucvNsoElemKernel() = delete;

  ScalarFieldType *vofNm1_{nullptr};
  ScalarFieldType *vofN_{nullptr};
  ScalarFieldType *vofNp1_{nullptr};
  VectorFieldType *velocityRTM_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  const double sucvFac_;
  const double nsoFac_;

  const int* lrscv_;

  double dt_{0.0};
  double gamma1_{0.0};
  double gamma2_{0.0};
  double gamma3_{0.0};

  const double Cupw_{0.1};
  const double small_{1.0e-16};

  // correction for gij for elements that do not span -1:1
  double gijFac_{1.0};

  /// Integration point to node mapping
  const int* ipNodeMap_;

  /// Shape functions
  AlignedViewType<DoubleType[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_ {"view_shape_func"};
};

}  // nalu
}  // sierra

#endif /* VolumeOfFluidSucvNsoElemKernel_H */
