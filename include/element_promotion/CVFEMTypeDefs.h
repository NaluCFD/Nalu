/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef CVFEMTypeDefs_h
#define CVFEMTypeDefs_h

#include <SupplementalAlgorithm.h>
#include <AlgTraits.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

// Kokkos
#include <Kokkos_Core.hpp>

namespace sierra{
namespace nalu{
  
  template <typename T, typename real_type = double>
  using vector_view = Kokkos::View<real_type[T::nDim_]>;

  template <typename T, typename real_type = double>
  using tensor_view = Kokkos::View<real_type[T::nDim_][T::nDim_]>;

  template <typename T, typename real_type = double>
  using nodal_scalar_view = Kokkos::View<real_type[T::nodes1D_][T::nodes1D_]>;

  template <typename T, typename real_type = double>
  using nodal_vector_view = Kokkos::View<real_type[T::nDim_][T::nodes1D_][T::nodes1D_]>;

  template <typename T, typename real_type = double>
  using nodal_tensor_view = Kokkos::View<real_type[T::nDim_][T::nDim_][T::nodes1D_][T::nodes1D_]>;

  template <typename T, typename real_type = double>
  using scs_scalar_view = Kokkos::View<real_type[T::nscs_][T::nodes1D_]>;

  template <typename T, typename real_type = double>
  using scs_vector_view = Kokkos::View<real_type[T::nDim_][T::nscs_][T::nodes1D_]>;

  template <typename T, typename real_type = double>
  using scs_tensor_view = Kokkos::View<real_type[T::nDim_][T::nDim_][T::nscs_][T::nodes1D_]>;

  template <typename T, typename real_type = double>
  using matrix_view = Kokkos::View<real_type[T::nodesPerElement_][T::nodesPerElement_]>;

  template <int p, typename real_type = double>
  using nodal_matrix_view = Kokkos::View<real_type[p+1][p+1]>;

  template <int p, typename real_type = double>
  using scs_matrix_view = Kokkos::View<real_type[p+1][p+1]>;

  template <int p, typename real_type = double>
  using linear_nodal_matrix_view = Kokkos::View<real_type[2][p+1]>;

  template <int p, typename real_type = double>
  using linear_scs_matrix_view = Kokkos::View<real_type[2][p]>;

} // namespace nalu
} // namespace Sierra

#endif
