/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef CONTINUITYADVELEMKERNEL_H
#define CONTINUITYADVELEMKERNEL_H

#include "Kernel.h"
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
class ScratchViews;

/** CMM (BDF2) for continuity equation (pressure DOF)
 */
template<typename AlgTraits>
class ContinuityAdvElemKernel: public Kernel
{
public:
  ContinuityAdvElemKernel(
    const stk::mesh::BulkData&,
    SolutionOptions&,
    ElemDataRequests&);

  virtual ~ContinuityAdvElemKernel();

  /** Perform pre-timestep work for the computational kernel
   */
  virtual void setup(const TimeIntegrator&);

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<double**>&,
    SharedMemView<double*>&,
    stk::mesh::Entity,
    ScratchViews&);

private:
  ContinuityAdvElemKernel() = delete;

  // extract fields; nodal
  VectorFieldType *velocityRTM_{nullptr};
  VectorFieldType *Gpdx_{nullptr};
  ScalarFieldType *pressure_{nullptr};
  ScalarFieldType *densityNp1_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  double projTimeScale_{1.0};

  const bool meshMotion_;
  const bool shiftMdot_;
  const bool shiftPoisson_;
  const bool reducedSensitivities_;
  const double interpTogether_;
  const double om_interpTogether_;

  // fixed size
  Kokkos::View<double[AlgTraits::nDim_]> v_uIp_{"view_uIp"};
  Kokkos::View<double[AlgTraits::nDim_]> v_rho_uIp_{"view_rhoUIp"};
  Kokkos::View<double[AlgTraits::nDim_]> v_Gpdx_Ip_{"view_GpdxIp"};
  Kokkos::View<double[AlgTraits::nDim_]> v_dpdxIp_{"view_dpdxIp"};

  // scratch space
  Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "view_shape_func" };

  const int* lrscv_;
};

}  // nalu
}  // sierra

#endif /* CONTINUITYADVELEMKERNEL_H */
