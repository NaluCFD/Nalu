/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef SCALARDIFFFEMKERNEL_H
#define SCALARDIFFFEMKERNEL_H

#include "Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>
#include <vector>

namespace sierra {
namespace nalu {

class MasterElement;
class Hex8FEM;
class ScratchViews;

/** CVFEM scalar advection/diffusion kernel
 */
template<typename AlgTraits>
class ScalarDiffFemKernel: public Kernel
{
public:
  ScalarDiffFemKernel(
    const stk::mesh::BulkData&,
    ScalarFieldType*,
    ScalarFieldType*);

  virtual ~ScalarDiffFemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<double**>&,
    SharedMemView<double*>&,
    stk::mesh::Entity,
    ScratchViews&);

private:
  ScalarDiffFemKernel() = delete;

  const stk::mesh::BulkData* bulkData_;

  ScalarFieldType *temperature_;
  ScalarFieldType *thermalCond_;
  VectorFieldType *coordinates_;

  // master element
  Hex8FEM * meFEM_;
  double *ipWeight_;

  // scratch space; geometry
  std::vector<double> ws_dndx_;
  std::vector<double> ws_deriv_;
  std::vector<double> ws_det_j_;
  std::vector<double> ws_shape_function_;

  // scratch space; fields
  std::vector<double> ws_temperature_;
  std::vector<double> ws_thermalCond_;
  std::vector<double> ws_coordinates_;
};

}  // nalu
}  // sierra

#endif /* SCALARDIFFFEMKERNEL_H */
