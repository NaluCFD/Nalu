/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SCALARUPWADVDIFFELEMKERNEL_H
#define SCALARUPWADVDIFFELEMKERNEL_H

#include "Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class SolutionOptions;
class MasterElement;
class ElemDataRequests;
class ScratchViews;
class PecletFunction;
class EquationSystem;

/** CVFEM scalar upwind advection/diffusion kernel
 */
template<typename AlgTraits>
class ScalarUpwAdvDiffElemKernel: public Kernel
{
public:
  ScalarUpwAdvDiffElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    EquationSystem*,
    ScalarFieldType*,
    VectorFieldType*,
    ScalarFieldType*,
    ElemDataRequests&);

  virtual ~ScalarUpwAdvDiffElemKernel();

  virtual void setup(const TimeIntegrator&);

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<double**>&,
    SharedMemView<double*>&,
    ScratchViews&);

  virtual double van_leer(
    const double &dqm,
    const double &dqp);

private:
  ScalarUpwAdvDiffElemKernel() = delete;

  const SolutionOptions& solnOpts_;

  ScalarFieldType *scalarQ_{nullptr};
  VectorFieldType *Gjq_{nullptr};
  ScalarFieldType *diffFluxCoeff_{nullptr};
  VectorFieldType *velocityRTM_{nullptr};
  ScalarFieldType *density_{nullptr};
  VectorFieldType *coordinates_{nullptr};
  GenericFieldType *massFlowRate_{nullptr};

  /// Left right node indicators
  const int* lrscv_;

  /// Name of the primitive variable (for upwind options lookup in solution options)
  const std::string dofName_;

  double alpha_;
  double alphaUpw_;
  double hoUpwind_;
  bool useLimiter_;
  double om_alpha_;
  double om_alphaUpw_;
  const bool shiftedGradOp_;
  const double small_{1.0e-16};

  /// Peclet function
  PecletFunction* pecletFunction_{nullptr};

  /// Shape functions
  Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "view_shape_func" };

  /// Scratch space to hold coordinates at the integration point
  Kokkos::View<double[AlgTraits::nDim_]> v_coordIp_{"v_coordIp"};
};

}  // nalu
}  // sierra

#endif /* SCALARUPWADVDIFFELEMKERNEL_H */
