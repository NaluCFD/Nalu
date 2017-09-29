/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "UnitTestKokkosME.h"

void check_that_values_match(const sierra::nalu::SharedMemView<DoubleType*>& values,
                             const std::vector<double>& oldValues)
{
  for(size_t i=0; i<values.dimension(0); ++i) {
      EXPECT_NEAR(stk::simd::get_data(values(i),0), oldValues[i], tol)<<"i:"<<i;
  }
}

void check_that_values_match(const sierra::nalu::SharedMemView<DoubleType**>& values,
                             const std::vector<double>& oldValues)
{
  int counter = 0;
  for(size_t i=0; i<values.dimension(0); ++i) {
    for(size_t j=0; j<values.dimension(1); ++j) {
      EXPECT_NEAR(stk::simd::get_data(values(i,j),0), oldValues[counter++], tol)<<"i:"<<i<<", j:"<<j;
    }
  }
}

void check_that_values_match(const sierra::nalu::SharedMemView<DoubleType***>& values,
                             const std::vector<double>& oldValues)
{
  int counter = 0;
  for(size_t i=0; i<values.dimension(0); ++i) {
    for(size_t j=0; j<values.dimension(1); ++j) {
      for(size_t k=0; k<values.dimension(2); ++k) {
        EXPECT_NEAR(stk::simd::get_data(values(i,j,k),0), oldValues[counter++], tol)<<"i:"<<i<<", j:"<<j<<", k:"<<k;
      }
    }
  }
}

void copy_DoubleType0_to_double(const sierra::nalu::SharedMemView<DoubleType**>& view,
                                std::vector<double>& vec)
{
  const DoubleType* viewValues = view.data();
  int len = view.size();
  vec.resize(len);
  for(int i=0; i<len; ++i) {
    vec[i] = stk::simd::get_data(viewValues[i], 0);
  }
}

void compare_old_scv_volume( const sierra::nalu::SharedMemView<DoubleType**>& v_coords,
                            const sierra::nalu::SharedMemView<DoubleType*>& scv_volume,
                            sierra::nalu::MasterElement* meSCV)
{
  int len = scv_volume.dimension(0);
  std::vector<double> coords;
  copy_DoubleType0_to_double(v_coords, coords);
  std::vector<double> volume(len, 0.0);
  double error = 0;
  meSCV->determinant(1, coords.data(), volume.data(), &error);
  EXPECT_NEAR(error, 0.0, tol);
  check_that_values_match(scv_volume, volume);
}

void compare_old_scs_areav( const sierra::nalu::SharedMemView<DoubleType**>& v_coords,
                            const sierra::nalu::SharedMemView<DoubleType**>& scs_areav,
                            sierra::nalu::MasterElement* meSCS)
{
  int len = scs_areav.dimension(0)*scs_areav.dimension(1);
  std::vector<double> coords;
  copy_DoubleType0_to_double(v_coords, coords);
  std::vector<double> areav(len, 0.0);
  double error = 0;
  meSCS->determinant(1, coords.data(), areav.data(), &error);
  EXPECT_NEAR(error, 0.0, tol);
  check_that_values_match(scs_areav, areav);
}

void compare_old_scs_grad_op( const sierra::nalu::SharedMemView<DoubleType**>& v_coords,
                            const sierra::nalu::SharedMemView<DoubleType***>& scs_dndx,
                            const sierra::nalu::SharedMemView<DoubleType***>& scs_deriv,
                            sierra::nalu::MasterElement* meSCS)
{
  int len = scs_dndx.dimension(0)*scs_dndx.dimension(1)*scs_dndx.dimension(2);
  std::vector<double> coords;
  copy_DoubleType0_to_double(v_coords, coords);
  std::vector<double> grad_op(len, 0.0);
  std::vector<double> deriv(len, 0.0);
  std::vector<double> det_j(len, 0.0);
  double error = 0;
  meSCS->grad_op(1, coords.data(), grad_op.data(), deriv.data(), det_j.data(), &error);
  EXPECT_NEAR(error, 0.0, tol);
  check_that_values_match(scs_dndx, grad_op);
  check_that_values_match(scs_deriv, deriv);
}

void compare_old_scs_shifted_grad_op( const sierra::nalu::SharedMemView<DoubleType**>& v_coords,
                            const sierra::nalu::SharedMemView<DoubleType***>& scs_dndx,
                            const sierra::nalu::SharedMemView<DoubleType***>& scs_deriv,
                            sierra::nalu::MasterElement* meSCS)
{
  int len = scs_dndx.dimension(0)*scs_dndx.dimension(1)*scs_dndx.dimension(2);
  std::vector<double> coords;
  copy_DoubleType0_to_double(v_coords, coords);
  std::vector<double> grad_op(len, 0.0);
  std::vector<double> deriv(len, 0.0);
  std::vector<double> det_j(len, 0.0);
  double error = 0;
  meSCS->shifted_grad_op(1, coords.data(), grad_op.data(), deriv.data(), det_j.data(), &error);
  EXPECT_NEAR(error, 0.0, tol);
  check_that_values_match(scs_dndx, grad_op);
  check_that_values_match(scs_deriv, deriv);
}

template<typename AlgTraits>
void test_ME_views(const std::vector<sierra::nalu::ELEM_DATA_NEEDED>& requests)
{
  unit_test_utils::KokkosMEViews<AlgTraits> driver(true);

  // Passing `true` to constructor has already initialized everything
  // driver.fill_mesh_and_init_data(/* doPerturb = */ false);

  // Register ME data requests
  for(sierra::nalu::ELEM_DATA_NEEDED request : requests) {
    driver.dataNeeded_.add_master_element_call(request, sierra::nalu::CURRENT_COORDINATES);
  }

  // Execute the loop and perform all tests
  driver.execute([&](sierra::nalu::ScratchViews<DoubleType>& scratchViews,
                     sierra::nalu::MasterElement* meSCS,
                     sierra::nalu::MasterElement* meSCV) {
      // Extract data from scratchViews
      sierra::nalu::SharedMemView<DoubleType**>& v_coords = scratchViews.get_scratch_view_2D(
        *driver.coordinates_);
      auto& meViews = scratchViews.get_me_views(sierra::nalu::CURRENT_COORDINATES);

      if (meSCS != nullptr) {
        for(sierra::nalu::ELEM_DATA_NEEDED request : requests) {
          if (request == sierra::nalu::SCS_AREAV) {
            compare_old_scs_areav(v_coords, meViews.scs_areav, meSCS);
          }
          if (request == sierra::nalu::SCS_GRAD_OP) {
            compare_old_scs_grad_op(v_coords, meViews.dndx, meViews.deriv, meSCS);
          }
          if (request == sierra::nalu::SCS_SHIFTED_GRAD_OP) {
            compare_old_scs_shifted_grad_op(v_coords, meViews.dndx_shifted, meViews.deriv, meSCS);
          }
        }
      }
      if (meSCV != nullptr) {
        for(sierra::nalu::ELEM_DATA_NEEDED request : requests) {
          if (request == sierra::nalu::SCV_VOLUME) {
            compare_old_scv_volume(v_coords, meViews.scv_volume, meSCV);
          }
        }
      }
    });
}

TEST(KokkosME, test_hex8_views)
{
  test_ME_views<sierra::nalu::AlgTraitsHex8>(
    {sierra::nalu::SCS_AREAV,
     sierra::nalu::SCS_GRAD_OP,
//     sierra::nalu::SCS_SHIFTED_GRAD_OP,
     sierra::nalu::SCS_GIJ,
     sierra::nalu::SCV_VOLUME,
    }
  );
}

TEST(KokkosME, test_tet4_views)
{
  test_ME_views<sierra::nalu::AlgTraitsTet4>(
    {sierra::nalu::SCS_AREAV,
     sierra::nalu::SCS_GRAD_OP,
     sierra::nalu::SCS_SHIFTED_GRAD_OP,
     sierra::nalu::SCV_VOLUME
    }
  );
}

