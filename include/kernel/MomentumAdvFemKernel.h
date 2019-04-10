/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MomentumAdvFemKernel_h
#define MomentumAdvFemKernel_h

#include "kernel/Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class ElemDataRequests;
class SolutionOptions;
class TimeIntegrator;

/** FEM advection kernel
 */
template<typename AlgTraits>
class MomentumAdvFemKernel: public Kernel
{
public:
  MomentumAdvFemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    VectorFieldType*,
    ElemDataRequests&);

  virtual ~MomentumAdvFemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

 private:
  MomentumAdvFemKernel() = delete;
  
  /** Perform pre-timestep work for the computational kernel
   */
  virtual void setup(const TimeIntegrator&);
  
  VectorFieldType *velocityNp1_{nullptr};
  ScalarFieldType *density_{nullptr};
  VectorFieldType *coordinates_{nullptr};
  
  VectorFieldType *vrtmL_{nullptr};
  VectorFieldType *GjpL_{nullptr};  
  ScalarFieldType *pressure_{nullptr};

  const bool shiftedGradOp_;
  const double includePstab_;
  const double skewFac_;
  const double om_skewFac_;
  double projTimeScale_;

  /// Shape functions
  AlignedViewType<DoubleType[AlgTraits::numGp_]> v_ip_weight_{ "v_ip_weight" };
  AlignedViewType<DoubleType[AlgTraits::numGp_][AlgTraits::nodesPerElement_]> v_shape_function_ { "v_shape_func" };
};

}  // nalu
}  // sierra

#endif /* MomentumAdvFemKernel_h */
