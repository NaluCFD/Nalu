/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef Hex27NGP_h
#define Hex27NGP_h

#include <AlgTraits.h>

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <array>
#include <type_traits>

#include <master_element/MasterElement.h>
#include <SimdInterface.h>
#include <Kokkos_Core.hpp>
#include <element_promotion/MasterElementUtils.h>
#include <master_element/MasterElementFunctions.h>
#include <stk_util/environment/ReportHandler.hpp>

namespace sierra{
namespace nalu{

class Hex27SCVNGP : public Hex27SCV
{
public:
  using AlgTraits = AlgTraitsHex27;
  using InterpWeightType = Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]>;
  using GradWeightType = Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_][AlgTraits::nDim_]>;

  Hex27SCVNGP()
  : Hex27SCV(),
    interpWeights_(copy_interpolation_weights_to_view(shapeFunctions_)),
    referenceGradWeights_(copy_deriv_weights_to_view(shapeDerivs_)),
    shiftedInterpWeights_(copy_interpolation_weights_to_view(shapeFunctionsShift_)),
    shiftedReferenceGradWeights_(copy_deriv_weights_to_view(shapeDerivsShift_))
  {};

  template <typename CoordViewType, typename OutputViewType>
  void determinant(CoordViewType coords, OutputViewType detj)
  {
    generic_determinant_3d<AlgTraits>(referenceGradWeights_, coords, detj);
    for (int ip = 0 ; ip < AlgTraits::numScsIp_; ++ip) {
      detj(ip) *= ipWeight_[ip];
    }
  }

private:
  InterpWeightType copy_interpolation_weights_to_view(std::vector<double>& interps)
  {
    InterpWeightType interpWeights{"interpolation_weights"};

    int shape_index = 0;
    for (int ip = 0; ip < AlgTraits::numScsIp_; ++ip) {
      for (int n = 0; n < AlgTraits::nodesPerElement_; ++n) {
        interpWeights(ip, n) = interps[shape_index];
        ++shape_index;
      }
    }
    return interpWeights;
  }

  GradWeightType copy_deriv_weights_to_view(std::vector<double>& derivs)
  {
    GradWeightType referenceGradWeights{"reference_gradient_weights"};

    int deriv_index = 0;
    for (int ip = 0; ip < AlgTraits::numScsIp_; ++ip) {
      for (int n = 0; n < AlgTraits::nodesPerElement_; ++n) {
        for (int d = 0; d < AlgTraits::nDim_; ++d) {
          referenceGradWeights(ip,n,d) = derivs[deriv_index];
          ++deriv_index;
        }
      }
    }
    return referenceGradWeights;
  }

  const InterpWeightType interpWeights_;
  const GradWeightType referenceGradWeights_;
  const InterpWeightType shiftedInterpWeights_;
  const GradWeightType shiftedReferenceGradWeights_;
};

class Hex27SCSNGP : public Hex27SCS
{
  using AlgTraits = AlgTraitsHex27;

  using InterpWeightType = Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]>;
  using GradWeightType = Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_][AlgTraits::nDim_]>;

public:
  Hex27SCSNGP()
: Hex27SCS(),
  interpWeights_(copy_interpolation_weights_to_view(shapeFunctions_)),
  referenceGradWeights_(copy_deriv_weights_to_view(shapeDerivs_)),
  shiftedInterpWeights_(copy_interpolation_weights_to_view(shapeFunctionsShift_)),
  shiftedReferenceGradWeights_(copy_deriv_weights_to_view(shapeDerivsShift_))
{};

  const InterpWeightType shape_functions()
  { return interpWeights_; }

  template <typename OutputViewType> void shape_functions(OutputViewType out)
  { Kokkos::deep_copy(out, interpWeights_); }

  const InterpWeightType shifted_functions()
  { return shiftedInterpWeights_; }

  template <typename OutputViewType> void shifted_shape_functions(OutputViewType out)
  { Kokkos::deep_copy(out, shiftedInterpWeights_); }

  template <typename CoordViewType, typename OutputViewType>
  void grad_op(CoordViewType coords, OutputViewType weights)
  { generic_grad_op_3d<AlgTraits>(referenceGradWeights_, coords, weights); }

  template <typename CoordViewType, typename OutputViewType>
  void shifted_grad_op(CoordViewType coords, OutputViewType weights)
  { generic_grad_op_3d<AlgTraits>(shiftedReferenceGradWeights_, coords, weights); }

  template <typename CoordViewType, typename OutputViewType>
  void area_vectors(CoordViewType coords, OutputViewType areav)
  { area_vectors(referenceGradWeights_, coords, areav); }

  template <typename CoordViewType, typename OutputViewType>
  void shifted_area_vectors(CoordViewType coords, OutputViewType areav)
  { area_vectors(shiftedReferenceGradWeights_, coords, areav); }

  template <typename CoordViewType, typename OutputViewType>
  void gij(CoordViewType coords, OutputViewType gupper, OutputViewType glower)
  { generic_gij_3d<AlgTraits>(referenceGradWeights_, coords, gupper, glower); }

  template <typename CoordViewType, typename OutputViewType>
  void shifted_gij(CoordViewType coords, OutputViewType gupper, OutputViewType glower)
  { generic_gij_3d<AlgTraits>(shiftedReferenceGradWeights_, coords, gupper, glower); }

private:
  template <int direction, typename CoordViewType, typename OutputViewType>
  void area_vector(int ip,  GradWeightType referenceGradWeights, CoordViewType coords, OutputViewType areav)
  {
    constexpr int s1Component = (direction == Jacobian::T_DIRECTION) ? Jacobian::S_DIRECTION : Jacobian::T_DIRECTION;
    constexpr int s2Component = (direction == Jacobian::U_DIRECTION) ? Jacobian::S_DIRECTION : Jacobian::U_DIRECTION;

    using ftype = typename CoordViewType::value_type;

    static_assert(std::is_same<ftype, typename OutputViewType::value_type>::value,
      "Incompatiable value type for views");

    static_assert(CoordViewType::Rank == 2, "Coordinate view assumed to be 2D");
    static_assert(OutputViewType::Rank == 2, "areav view assumed to be 2D");

    ftype jac[AlgTraits::nDim_][AlgTraits::nDim_-1] = { ftype(0) };

    for (int n = 0; n < AlgTraits::nodesPerElement_; ++n) {
      const ftype dn_ds1 = referenceGradWeights(ip, n, s1Component);
      const ftype dn_ds2 = referenceGradWeights(ip, n, s2Component);

      for (int d = 0; d < AlgTraits::nDim_; ++d) {
        jac[d][0] += dn_ds1 * coords(n,d);
        jac[d][1] += dn_ds2 * coords(n,d);
      }
    }

    //cross product
    areav(ip, 0) = jac[1][0] * jac[2][1] - jac[2][0] * jac[1][2];
    areav(ip, 1) = jac[2][0] * jac[0][1] - jac[0][0] * jac[2][1];
    areav(ip, 2) = jac[0][0] * jac[1][1] - jac[1][0] * jac[0][1];
  }

  template <typename CoordViewType, typename OutputViewType>
  void area_vectors(GradWeightType referenceGradWeights, CoordViewType coords, OutputViewType areav)
  {
    using ftype = typename CoordViewType::value_type;

    static_assert(std::is_same<ftype, typename OutputViewType::value_type>::value, "Incompatiable value type for views");
    static_assert(CoordViewType::Rank == 2, "Coordinate view assumed to be 2D");
    static_assert(OutputViewType::Rank == 2, "area_vector view assumed to be 2D");

    static_assert (AlgTraits::numScsIp_ % AlgTraits::nDim_ == 0, "Number of ips incorrect");
    constexpr int ipsPerDirection = AlgTraits::numScsIp_ / AlgTraits::nDim_;

    // this relies on the ips being laid out direction-by-direction,
    // specifically in the U->T->S order
    for (int ip = 0; ip < ipsPerDirection; ++ip) {
      ThrowAssert(ipInfo_[ip].direction == Jacobian::U_DIRECTION);
      area_vector<Jacobian::U_DIRECTION>(ip, referenceGradWeights, coords, areav);
    }

    for (int ip = ipsPerDirection; ip < 2 * ipsPerDirection; ++ip) {
      ThrowAssert(ipInfo_[ip].direction == Jacobian::T_DIRECTION);
      area_vector<Jacobian::T_DIRECTION>(ip, referenceGradWeights, coords, areav);
    }

    for (int ip = 2 * ipsPerDirection; ip < 3 * ipsPerDirection; ++ip) {
      ThrowAssert(ipInfo_[ip].direction == Jacobian::S_DIRECTION);
      area_vector<Jacobian::S_DIRECTION>(ip, referenceGradWeights, coords, areav);
    }

    for (int ip = 0; ip < AlgTraits::numScsIp_; ++ip) {
      ftype weight = ipInfo_[ip].weight;
      areav(ip, 0) *= weight;
      areav(ip, 1) *= weight;
      areav(ip, 2) *= weight;
    }
  }

  InterpWeightType copy_interpolation_weights_to_view(std::vector<double>& interps)
  {
    InterpWeightType interpWeights{"interpolation_weights"};

    int shape_index = 0;
    for (int ip = 0; ip < AlgTraits::numScsIp_; ++ip) {
      for (int n = 0; n < AlgTraits::nodesPerElement_; ++n) {
        interpWeights(ip, n) = interps[shape_index];
        ++shape_index;
      }
    }
    return interpWeights;
  }

  GradWeightType copy_deriv_weights_to_view(std::vector<double>& derivs)
  {
    GradWeightType referenceGradWeights{"reference_gradient_weights"};

    int deriv_index = 0;
    for (int ip = 0; ip < AlgTraits::numScsIp_; ++ip) {
      for (int n = 0; n < AlgTraits::nodesPerElement_; ++n) {
        for (int d = 0; d < AlgTraits::nDim_; ++d) {
          referenceGradWeights(ip,n,d) = derivs[deriv_index];
          ++deriv_index;
        }
      }
    }
    return referenceGradWeights;
  }

  const InterpWeightType interpWeights_;
  const GradWeightType referenceGradWeights_;
  const InterpWeightType shiftedInterpWeights_;
  const GradWeightType shiftedReferenceGradWeights_;
};

} // namespace nalu
} // namespace Sierra

#endif
