/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "UnitTestKokkosME.h"

TEST(KokkosME_Hex8, test_hex8_views)
{
  unit_test_utils::KokkosMEViews<sierra::nalu::AlgTraitsHex8> driver(true);

  // Passing `true` to constructor has already initialized everything
  // driver.fill_mesh_and_init_data(/* doPerturb = */ false);

  // Register ME data requests
  driver.dataNeeded_.add_master_element_call(
    sierra::nalu::SCS_AREAV, sierra::nalu::CURRENT_COORDINATES);

  // Execute the loop and perform all tests
  driver.execute([&](sierra::nalu::ScratchViews<DoubleType>& scratchViews) {
      // Extract data from scratchViews
      sierra::nalu::SharedMemView<DoubleType**>& v_coords = scratchViews.get_scratch_view_2D(
        *driver.coordinates_);
      auto& meViews = scratchViews.get_me_views(sierra::nalu::CURRENT_COORDINATES);
      auto& v_scs_areav = meViews.scs_areav;

      // Perform tests here...
      EXPECT_EQ(v_coords.dimension(0), 8);
      EXPECT_EQ(v_coords.dimension(1), 3);
      EXPECT_EQ(v_scs_areav.dimension(0), 12);
    });
}
