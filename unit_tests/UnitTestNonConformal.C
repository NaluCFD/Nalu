#include <gtest/gtest.h>
#include <NaluEnv.h>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_topology/topology.hpp>

#include <element_promotion/PromotedPartHelper.h>
#include <element_promotion/ElementDescription.h>
#include <element_promotion/PromoteElement.h>
#include <element_promotion/PromotedElementIO.h>
#include <nalu_make_unique.h>

#include "UnitTestUtils.h"
#include "UnitTestKokkosUtils.h"


#include <stk_unit_tests/stk_mesh_fixtures/HexFixture.hpp>
#include <stk_unit_tests/stk_mesh_fixtures/TetFixture.hpp>
#include <stk_unit_tests/stk_mesh_fixtures/WedgeFixture.hpp>
#include <stk_unit_tests/stk_mesh_fixtures/PyramidFixture.hpp>

#include <algorithm>
#include <string>
#include <array>
#include <random>
#include "UnitTestNonConformalMeshHelperUtil.h"

TEST(NonConformal, TetTet) {
  using TetFix = stk::mesh::fixtures::TetFixture;

  int nA = 5;
  int nB = 8;
  auto meshHelper = sierra::nalu::make_unique<unit_test_utils::NonConformalMeshHelper<TetFix,TetFix>>(nA, nB);
  const auto& dgInterface = meshHelper->dgInterface;
  stk::mesh::BulkData& bulk = meshHelper->bulk_;
  const stk::mesh::MetaData& meta = meshHelper->meta_;

  ASSERT_NE(meta.get_part("block_1"), nullptr);
  stk::topology topoFromBlock1 = meta.get_part("block_1")->topology();
  stk::topology topoA = dgInterface.first.elemTopo;

  ASSERT_NE(meta.get_part("block_2"), nullptr);
  stk::topology topoFromBlock2 = meta.get_part("block_2")->topology();
  stk::topology topoB = dgInterface.second.elemTopo;

  stk::topology tetTopo = stk::topology::TETRAHEDRON_4;

  EXPECT_EQ(topoA.value(), tetTopo.value());
  EXPECT_EQ(topoA.value(), topoFromBlock1.value());
  EXPECT_EQ(topoB.value(), tetTopo.value());
  EXPECT_EQ(topoB.value(), topoFromBlock2.value());

  if (bulk.parallel_size() == 1) {
    const auto& elem_buckets = bulk.get_buckets(stk::topology::ELEM_RANK, meta.universal_part());
    size_t elem_count = 0u;
    for(const auto* ib : elem_buckets) {
      elem_count += ib->size();
    }
    const size_t numTetsPerHex = 6;
    const size_t expectedElemCount = numTetsPerHex*(nA*nA*nA) + numTetsPerHex*(nB*nB*nB);
    EXPECT_EQ(elem_count, expectedElemCount);
  }

  std::string fileName = "nonconformal_" + topoA.name() + "_" + topoB.name() + ".g";
  EXPECT_NO_THROW(unit_test_utils::dump_mesh(bulk, {}, fileName));
}


