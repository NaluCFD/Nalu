/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MOMENTUMUPWADVDIFFELEMKERNEL_H
#define MOMENTUMUPWADVDIFFELEMKERNEL_H

#include "kernel/Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class ElemDataRequests;
class EquationSystem;
class MasterElement;
template <typename T> class PecletFunction;
class SolutionOptions;

/** Advection (with upwind) diffusion term for momentum equation (velocity DOF)
 */
template<typename AlgTraits>
class MomentumUpwAdvDiffElemKernel: public Kernel
{
public:
  MomentumUpwAdvDiffElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    EquationSystem*,
    VectorFieldType*,
    ScalarFieldType*,
    GenericFieldType*,
    ElemDataRequests&);

  virtual ~MomentumUpwAdvDiffElemKernel();

  virtual void setup(const TimeIntegrator&);

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

  virtual DoubleType van_leer(
    const DoubleType &dqm,
    const DoubleType &dqp);

private:
  MomentumUpwAdvDiffElemKernel() = delete;

  const SolutionOptions &solnOpts_;
  
  VectorFieldType *velocityNp1_{nullptr};
  VectorFieldType *coordinates_{nullptr};
  ScalarFieldType *density_{nullptr};
  ScalarFieldType *viscosity_{nullptr};
  GenericFieldType *Gju_{nullptr};
  GenericFieldType *massFlowRate_{nullptr};
  VectorFieldType *velocityRTM_{nullptr};

  const int* lrscv_;

  /// Name of the primitive variable (for upwind options lookup in solution options)
  const std::string dofName_;

  double alpha_;
  double alphaUpw_;
  double hoUpwind_;
  bool useLimiter_;
  double om_alpha_;
  double om_alphaUpw_;
  const double includeDivU_;
  const bool shiftedGradOp_;
  const double small_{1.0e-16};
  
  // Peclet function
  PecletFunction<DoubleType>* pecletFunction_{nullptr};

  // fixed scratch space
  AlignedViewType<DoubleType[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
  AlignedViewType<DoubleType[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_adv_shape_function_{"v_adv_shape_function"};
};

}  // nalu
}  // sierra

#endif /* MOMENTUMUPWADVDIFFELEMKERNEL_H */
