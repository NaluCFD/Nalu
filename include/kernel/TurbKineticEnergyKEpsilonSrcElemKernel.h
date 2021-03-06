/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef TurbKineticEnergyKEpsilonSrcElemKernel_H
#define TurbKineticEnergyKEpsilonSrcElemKernel_H

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

/** Add Ksgs source term for kernel-based algorithm approach
 */
template<typename AlgTraits>
class TurbKineticEnergyKEpsilonSrcElemKernel: public Kernel
{
public:
  TurbKineticEnergyKEpsilonSrcElemKernel(
    const stk::mesh::BulkData&,
    const SolutionOptions&,
    ElemDataRequests&,
    const bool);

  virtual ~TurbKineticEnergyKEpsilonSrcElemKernel();
  
  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<DoubleType**>&,
    SharedMemView<DoubleType*>&,
    ScratchViews<DoubleType>&);

private:
  TurbKineticEnergyKEpsilonSrcElemKernel() = delete;
  
  VectorFieldType *coordinates_{nullptr};
  VectorFieldType *velocityNp1_{nullptr};
  ScalarFieldType *tkeNp1_{nullptr};
  ScalarFieldType *epsNp1_{nullptr};
  ScalarFieldType *densityNp1_{nullptr};
  ScalarFieldType *tvisc_{nullptr};

  const bool lumpedMass_;
  const bool shiftedGradOp_;

  const double tkeProdLimitRatio_;
  const double includeDivU_;
  const double twoThirds_;

  // Integration point to node mapping
  const int* ipNodeMap_;
  
  // fixed scratch space
  AlignedViewType<DoubleType[AlgTraits::numScvIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
};
 
}  // nalu
}  // sierra

#endif /* TurbKineticEnergyKEpsilonSrcElemKernel_H */
