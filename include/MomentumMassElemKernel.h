/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MOMENTUMMASSELEMKERNEL_H
#define MOMENTUMMASSELEMKERNEL_H

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

/** CMM (BDF2/BE) for momentum equation (velocity DOF)
 */
template<typename AlgTraits>
class MomentumMassElemKernel: public Kernel
{
public:
  MomentumMassElemKernel(
    const stk::mesh::BulkData&,
    SolutionOptions&,
    ElemDataRequests&,
    const bool);

  virtual ~MomentumMassElemKernel();

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
  MomentumMassElemKernel() = delete;

  VectorFieldType *velocityNm1_{nullptr};
  VectorFieldType *velocityN_{nullptr};
  VectorFieldType *velocityNp1_{nullptr};
  ScalarFieldType *densityNm1_{nullptr};
  ScalarFieldType *densityN_{nullptr};
  ScalarFieldType *densityNp1_{nullptr};
  VectorFieldType *Gjp_{nullptr};
  VectorFieldType *coordinates_{nullptr};

  double dt_{0.0};
  double gamma1_{0.0};
  double gamma2_{0.0};
  double gamma3_{0.0};
  const bool lumpedMass_;

  /// Integration point to node mapping
  const int* ipNodeMap_;

  Kokkos::View<double[AlgTraits::nDim_]> v_uNm1_ {"v_uNm1"};
  Kokkos::View<double[AlgTraits::nDim_]> v_uN_   {"v_uN"};
  Kokkos::View<double[AlgTraits::nDim_]> v_uNp1_ {"v_uNp1"};
  Kokkos::View<double[AlgTraits::nDim_]> v_Gjp_  {"v_Gjp"};

  /// Shape functions
  Kokkos::View<double[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_ {"view_shape_func"};
};

}  // nalu
}  // sierra

#endif /* MOMENTUMMASSELEMKERNEL_H */
