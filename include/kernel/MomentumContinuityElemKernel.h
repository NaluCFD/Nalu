/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MomentumContinuityElemKernel_h
#define MomentumContinuityElemKernel_h

#include "Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class ElemDataRequests;
class MasterElement;
class SolutionOptions;
class TimeIntegrator;

/** Advection diffusion term for momentum equation (velocity DOF)
 */
template<typename AlgTraits>
class MomentumContinuityElemKernel: public Kernel
{
public:
  MomentumContinuityElemKernel(
    const stk::mesh::BulkData &bulkData,
    const SolutionOptions &solnOpts,
    VectorFieldType *velocity,
    ScalarFieldType *density,
    ScalarFieldType *viscosity,
    const bool lumpedMass,
    ElemDataRequests &dataPreReqs);

  virtual ~MomentumContinuityElemKernel();

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
  MomentumContinuityElemKernel() = delete;

  ScalarFieldType *densityNp1_{nullptr};
  ScalarFieldType *densityN_{nullptr};
  ScalarFieldType *densityNm1_{nullptr};
  VectorFieldType *velocityNp1_{nullptr};
  VectorFieldType *velocityN_{nullptr};
  VectorFieldType *velocityNm1_{nullptr};
  VectorFieldType *Gjp_{nullptr};
  VectorFieldType *GjpOld_{nullptr};
  VectorFieldType *coordinates_{nullptr};
  ScalarFieldType *pressure_{nullptr};
  ScalarFieldType *viscosity_{nullptr};

  double includeDivU_{0.0};
  double dt_{0.0};
  double gamma1_{0.0};
  double gamma2_{0.0};
  double gamma3_{0.0};
  double projTimeScale_{1.0};

  /// ip node mappings
  const int* lrscv_;
  const int* ipNodeMap_;

  const int dofSize_;

  // fixed scratch space
  AlignedViewType<DoubleType[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_scs_{"v_shape_function_scs"};
  AlignedViewType<DoubleType[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_scv_{"v_shape_function_scv"};  
};

}  // nalu
}  // sierra

#endif /* MomentumContinuityElemKernel_H */
