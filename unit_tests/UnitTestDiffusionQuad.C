#include <gtest/gtest.h>
#include <limits>

#include <element_promotion/operators/HighOrderOperatorsQuad.h>
#include <element_promotion/operators/HighOrderDiffusionQuad.h>
#include <element_promotion/operators/HighOrderGeometryQuadDiffusion.h>
#include <element_promotion/operators/HighOrderGeometryQuadVolume.h>
#include <element_promotion/operators/CoefficientMatrices.h>
#include <element_promotion/QuadratureRule.h>
#include <Kokkos_Core.hpp>
#include <Teuchos_BLAS.hpp>

#include <memory>
#include <tuple>
#include <random>

#include <element_promotion/ElementDescription.h>
#include "UnitTestUtils.h"

namespace {
  //--------------------------------------------------------------
  template <typename T, typename Y>
  void kron(const T& A,const T& B,const Y& C) {
    for (unsigned n = 0; n <  A.dimension_1(); ++n) {
      for (unsigned m = 0; m <  A.dimension_0(); ++m) {
        for (unsigned j = 0; j < A.dimension_1(); ++j) {
          for (unsigned i = 0; i < A.dimension_0(); ++i) {
            C(j*A.dimension_0()+i,n*A.dimension_0()+m) = A(n,j) * B(m,i);
          }
        }
      }
    }
  }
  //--------------------------------------------------------------
  template<typename T>
  T transpose(const T& A)
  {
    T At(A.label()+"_t");

    EXPECT_EQ(A.dimension_0(), A.dimension_1());
    for (unsigned j = 0; j < A.dimension_0(); ++j) {
      for (unsigned i = 0; i < A.dimension_0(); ++i) {
        At(i,j) = A(j,i);
      }
    }
    return At;
  }
  //--------------------------------------------------------------
  template <int p> void check_volume_metric()
  {
    Kokkos::View<double[p+1][p+1]> exactScvVolume{"vol"};
    std::vector<double> scsLocations1D = sierra::nalu::gauss_legendre_rule(p).first;
    std::vector<double> paddedScsLocations1D = sierra::nalu::pad_end_points(scsLocations1D,-1,+1); // add the element ends
    for (int j = 0; j < p+1; ++j) {
      double y_scsL = paddedScsLocations1D[j+0];
      double y_scsR = paddedScsLocations1D[j+1];
      for (int i = 0; i < p+1;++i) {
        double x_scsL = paddedScsLocations1D[i+0];
        double x_scsR = paddedScsLocations1D[i+1];
        exactScvVolume(j, i) = (x_scsR - x_scsL) * (y_scsR - y_scsL);
      }
    }

    Kokkos::View<double[2][p+1][p+1]> coords{"coords"};

    std::vector<double> coords1D = sierra::nalu::gauss_lobatto_legendre_rule(p+1).first;
    for (int j = 0; j < p+1; ++j) {
      double y = coords1D[j];
      for (int i = 0; i < p+1;++i) {
        double x = coords1D[i];
        coords(0,j,i) = x;
        coords(1,j,i) = y;
      }
    }

    auto ops = sierra::nalu::CVFEMQuadOperators<p>();
    Kokkos::View<double[p+1][p+1]> detj{""};
    sierra::nalu::high_order_metrics::compute_volume_metric_linear(ops, coords, detj);

    Kokkos::View<double[p+1][p+1]> computedScvVolume{""};
    ops.volume_2D(detj, computedScvVolume);
    EXPECT_VIEW_NEAR_2D(exactScvVolume, computedScvVolume, 1.0e-10);
  }
  //--------------------------------------------------------------
  template <int p> Kokkos::View<double[(p+1)*(p+1)][(p+1)*(p+1)]>
  single_square_element_laplacian()
  {
    auto mat = sierra::nalu::CoefficientMatrices<p>();

    auto diff = sierra::nalu::coefficients::difference_matrix<p>();
    auto deriv = sierra::nalu::coefficients::scs_derivative_weights<p>();
    auto integ = sierra::nalu::coefficients::nodal_integration_weights<p>();
    auto identity = sierra::nalu::coefficients::identity_matrix<p>();

    Kokkos::View<double[p + 1][p + 1]> laplacian1D{ "l1d" };
    Teuchos::BLAS<int, double>().GEMM(Teuchos::TRANS, Teuchos::TRANS,
      p + 1, p + 1, p + 1,
      1.0,
      diff.data(), p + 1, deriv.data(), p + 1,
      0.0,
      laplacian1D.data(),  p + 1);

    int ipiv[(p+1)*(p+1)] = {0};
    int info = 0;

    auto integ_t = transpose(integ);
    Teuchos::LAPACK<int,double>().GESV(p+1, p+1, integ_t.data(), p+1, ipiv, laplacian1D.data(), p+1, &info);
    EXPECT_EQ(info,0);

    constexpr int npe = (p+1)*(p+1);
    Kokkos::View<double[npe][npe]> lhsx{"m"};
    kron(identity, laplacian1D, lhsx ); //dxx

    Kokkos::View<double[npe][npe]> lhs{"laplacian"};
    kron(laplacian1D, identity, lhs); //dyy

    for (unsigned j = 0; j < lhs.size(); ++j) {
      lhs.data()[j] += lhsx.data()[j]; // dxx + dyy
    }

    return lhs;
  }

  template <int p>
  void check_diffusion()
  {
    struct MMSFunction {
      MMSFunction() : k(1.0), pi(std::acos(-1.0)) {};

      double val(double x, double y) {return (0.25*(std::cos(1.0*k*pi*x) + std::cos(2.0*k*pi*y))); }

      double laplacian(double x, double y) const {
        return ( -(k*pi)*(k*pi) * (std::cos(1.0*k*pi*x)/4. + std::cos(2.0*k*pi*y)) );
      };

      double k;
      double pi;
    };
    MMSFunction mms;

    Kokkos::View<double[2][p+1][p+1]> coords{"x"};
    Kokkos::View<double[p+1][p+1]> diffusivity{"d"};
    Kokkos::View<double[p+1][p+1]> scalar{"q"};
    Kokkos::View<double[p+1][p+1]> laplacian_of_scalar{"d2q"};

    std::vector<double> coords1D = sierra::nalu::gauss_lobatto_legendre_rule(p+1).first;
    for (int j = 0; j < p+1; ++j) {
      for (int i = 0; i < p+1;++i) {
        double x = coords1D[i];
        double y = coords1D[j];
        coords(0, j,i) = x;
        coords(1, j,i) = y;
        diffusivity(j,i) = 1.0;
        scalar(j,i) = mms.val(x,y);
        laplacian_of_scalar(j,i) = mms.laplacian(x,y);
      }
    }
    constexpr int npe = (p+1)*(p+1);
    Kokkos::View<double[npe][npe]> mass{"m"};
    Kokkos::View<double[p+1][p+1]> detj{""};
    auto ops = sierra::nalu::CVFEMQuadOperators<p>();
    sierra::nalu::high_order_metrics::compute_volume_metric_linear(ops, coords, detj);
    kron(ops.mat_.nodalWeights, ops.mat_.nodalWeights, mass);

    Kokkos::View<double[p+1][p+1]> m_times_laplacian_of_scalar{"rhs_e"};
    Teuchos::BLAS<int,double>().GEMV(
      Teuchos::NO_TRANS,
      npe, npe,
      +1.0,
      mass.ptr_on_device(), npe,
      laplacian_of_scalar.ptr_on_device(), 1,
      +0.0,
      m_times_laplacian_of_scalar.ptr_on_device(), 1
    );

    Kokkos::View<double[2][2][p][p+1]> AGJ{"AGJ"};
    sierra::nalu::high_order_metrics::compute_diffusion_metric_linear(ops, coords, diffusivity, AGJ);

    Kokkos::View<double[npe][npe]> lhs{"lhs"};
    sierra::nalu::tensor_assembly::elemental_diffusion_jacobian(ops, AGJ, lhs);

    Kokkos::View<double[p+1][p+1]> rhs{"rhs_l"};
    Teuchos::BLAS<int,double>().GEMV(
      Teuchos::TRANS, // row v column
      npe, npe,
      -1.0,
      lhs.ptr_on_device(), npe,
      scalar.ptr_on_device(), 1,
      +1.0,
      rhs.ptr_on_device(), 1
    );
    EXPECT_VIEW_NEAR_2D(m_times_laplacian_of_scalar, rhs, 1.0e-8);

    Kokkos::View<double[p+1][p+1]> rhs_jf{"rhs_j"};
    sierra::nalu::tensor_assembly::elemental_diffusion_action(ops, AGJ, scalar, rhs_jf);
    EXPECT_VIEW_NEAR_2D(rhs, rhs_jf, 1.0e-8);

    auto laplacian_operator = single_square_element_laplacian<p>();
    Kokkos::View<double[p+1][p+1]> numerical_laplacian{"num_l"};
    Teuchos::BLAS<int,double>().GEMV(
      Teuchos::TRANS, // row v column
      npe, npe,
      -1.0,
      laplacian_operator.ptr_on_device(), npe,
      scalar.ptr_on_device(), 1,
      +0.0,
      numerical_laplacian.ptr_on_device(), 1
    );
    EXPECT_VIEW_NEAR_2D(laplacian_of_scalar, numerical_laplacian, 1.0e-8);
  }
}

//--------------------------------------------------------------
#define TEST_POLY(x,y,z)  TEST(x, y##_##order_##z) { y<z>(); }
TEST_POLY(QuadDiffusion, check_volume_metric, 42);
TEST_POLY(QuadDiffusion, check_diffusion, 42);

