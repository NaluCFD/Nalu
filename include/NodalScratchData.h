/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef NodalScratchData_h
#define NodalScratchData_h

#include <SolverAlgorithm.h>
#include <ElemDataRequests.h>
#include <FieldTypeDef.h>
#include <nalu_make_unique.h>

#include <KokkosInterface.h>
#include <SimdInterface.h>

namespace sierra{
namespace nalu{

class Realm;
class MasterElement;

template <typename T> class NodalScratchData
{
  /** Data holder for nodal loops.  Allows one to request and get gathered local data structures for a node */

public:
  typedef T value_type;
  NodalScratchData(const TeamHandleType& team, const stk::mesh::BulkData& bulk, ElemDataRequests& dataRequests)
  {
    const int num_fields = bulk.mesh_meta_data().get_fields().size();
    data0D_.resize(num_fields);
    data1D_.resize(num_fields);
    data2D_.resize(num_fields);

    const auto& neededFields = dataRequests.get_fields();
    for(const FieldInfo& fieldInfo : neededFields) {
      const unsigned scalarsDim1 = fieldInfo.scalarsDim1;
      const unsigned scalarsDim2 = fieldInfo.scalarsDim2;
      if (scalarsDim1 > 1u && scalarsDim2 <= 1u) {
        data1D_[fieldInfo.field->mesh_meta_data_ordinal()] = get_shmem_view_1D<T>(team, scalarsDim1);
      }
      if (scalarsDim1 > 1u && scalarsDim2 > 1u) {
        data2D_[fieldInfo.field->mesh_meta_data_ordinal()] = get_shmem_view_2D<T>(team, scalarsDim1, scalarsDim2);
      }
    }
  }

  T& get_value(const stk::mesh::FieldBase& field) { return data0D_[field.mesh_meta_data_ordinal()]; };
  SharedMemView<T*> get_1D(const stk::mesh::FieldBase& field)  { return data1D_[field.mesh_meta_data_ordinal()]; };
  SharedMemView<T**> get_2D(const stk::mesh::FieldBase& field) { return data2D_[field.mesh_meta_data_ordinal()]; };

  const T& get_value(const stk::mesh::FieldBase& field) const { return data0D_[field.mesh_meta_data_ordinal()]; };
  SharedMemView<const T*> get_1D(const stk::mesh::FieldBase& field) const { return data1D_[field.mesh_meta_data_ordinal()]; };
  SharedMemView<const T**> get_2D(const stk::mesh::FieldBase& field) const { return data2D_[field.mesh_meta_data_ordinal()]; };

private:
  ScalarAlignedVector data0D_;
  std::vector<SharedMemView<T*>> data1D_;
  std::vector<SharedMemView<T**>> data2D_;
};


template <typename T>
const int bytes_needed_for_nodal_kernels(ElemDataRequests& neededNodalData)
{
  int numScalars = 0;
  const auto& neededFields = neededNodalData.get_fields();
  for(const FieldInfo& fieldInfo : neededFields) {
    unsigned scalarsDim1 = fieldInfo.scalarsDim1;
    unsigned scalarsDim2 = fieldInfo.scalarsDim2;
    if (scalarsDim1 <= 1u && scalarsDim2 <= 1u) {
      numScalars += 1;
    } else if (scalarsDim1 > 1u && scalarsDim2 <= 1u) {
      numScalars += scalarsDim1;
    }
    else if (scalarsDim1 > 1u && scalarsDim2 > 1u) {
      numScalars += scalarsDim1 * scalarsDim2;
    }
  }
  return numScalars * sizeof(T);
}


KOKKOS_INLINE_FUNCTION int bytes_needed_for_nodal_utility_kernels(ElemDataRequests& neededInputs, ElemDataRequests& neededOutputs)
{
  int bytes_per_thread = bytes_needed_for_nodal_kernels<double>(neededInputs)
                      +  bytes_needed_for_nodal_kernels<double>(neededOutputs);
  bytes_per_thread *= 2 * stk::simd::ndoubles;
  return bytes_per_thread;
}



} // namespace nalu
} // namespace Sierra

#endif

