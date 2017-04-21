#ifndef _UnitTestViewUtils_h_
#define _UnitTestViewUtils_h_

#include <gtest/gtest.h>
#include <string>
#include <ostream>
#include <Kokkos_Core.hpp>

namespace unit_test_utils {

template <typename T>
void dump_2d_view(const T& v, bool clipzero = true)
{
  for (unsigned j = 0; j < v.dimension_1(); ++j) {
    for (unsigned i = 0; i < v.dimension_0(); ++i) {
      double vout = (std::abs(v(j,i)) < 1.0e-15  && clipzero) ? 0.0 : v(j,i);
      std::cout << "(" << i << ", " << j << ")" << ", " << std::setw(5) << v.label() << ": " << std::setw(12) << vout;
      if (i != v.dimension_0()-1) { std::cout << ", "; }
    }
    std::cout << std::endl;
  }
  std::cout << "--------------" << std::endl;
}

}

#define EXPECT_VIEW_NEAR_1D(x,y,tol)                 \
{                                                    \
  EXPECT_EQ(x.dimension_0(), y.dimension_0());       \
  for (unsigned i = 0; i < x.dimension_0(); ++i) {   \
    EXPECT_NEAR(x(i), y(i), tol);                    \
  }                                                  \
}

#define EXPECT_VIEW_NEAR_2D(x,y,tol)                  \
{                                                     \
  EXPECT_EQ(x.dimension_0(), y.dimension_0());        \
  for (unsigned j = 0; j < x.dimension_0();++j) {     \
    for (unsigned i = 0; i < x.dimension_1();++i) {   \
      EXPECT_NEAR(x(j,i), y(j,i), tol);               \
    }                                                 \
  }                                                   \
}

#define EXPECT_VIEW_NEAR_3D(x,y,tol)                   \
{                                                      \
  EXPECT_EQ(x.dimension_0(), y.dimension_0());         \
  EXPECT_EQ(x.dimension_1(), y.dimension_1());         \
  EXPECT_EQ(x.dimension_2(), y.dimension_2());         \
  for (unsigned k = 0; k < x.dimension_0();++k) {      \
    for (unsigned j = 0; j < x.dimension_1();++j) {    \
      for (unsigned i = 0; i < x.dimension_2();++i) {  \
        EXPECT_NEAR(x(k,j,i), y(k,j,i), tol);          \
      }                                                \
    }                                                  \
  }                                                    \
}

#define EXPECT_VIEW_NEAR_4D(x,y,tol)                     \
{                                                        \
  EXPECT_EQ(x.dimension_0(), y.dimension_0());           \
  EXPECT_EQ(x.dimension_1(), y.dimension_1());           \
  EXPECT_EQ(x.dimension_2(), y.dimension_2());           \
  EXPECT_EQ(x.dimension_3(), y.dimension_3());           \
  for (unsigned l = 0; l < x.dimension_0(); ++l) {       \
    for (unsigned k = 0; k < x.dimension_1(); ++k) {     \
      for (unsigned j = 0; j < x.dimension_2(); ++j) {   \
        for (unsigned i = 0; i < x.dimension_3(); ++i) { \
          EXPECT_NEAR(x(l,k,j,i), y(l,k,j,i), tol);      \
        }                                                \
      }                                                  \
    }                                                    \
  }                                                      \



#endif

