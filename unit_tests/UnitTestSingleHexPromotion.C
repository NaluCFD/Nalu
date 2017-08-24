#include <gtest/gtest.h>
#include <limits>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/SkinMesh.hpp>

#include <master_element/MasterElementHO.h>
#include <element_promotion/PromotedPartHelper.h>
#include <element_promotion/PromoteElement.h>
#include <element_promotion/PromotedElementIO.h>

#include <nalu_make_unique.h>
#include <NaluEnv.h>
#include <BucketLoop.h>

#include <memory>
#include <tuple>
#include <random>

#include <element_promotion/ElementDescription.h>
#include "UnitTestUtils.h"

TEST(SingleHexPromotion, coords_p2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }

  // Hex 27 standard node locations for a [0,1]^3 element, with the center node moved from index 20 to index 26.
  std::vector<std::vector<double>> expectedCoords =
  {
      {+0.0,+0.0,+0.0}, {+1.0,+0.0,+0.0}, {+1.0,+1.0,+0.0}, {+0.0,+1.0,+0.0},
      {+0.0,+0.0,+1.0}, {+1.0,+0.0,+1.0}, {+1.0,+1.0,+1.0}, {+0.0,+1.0,+1.0},
      {+0.5,+0.0,+0.0}, {+1.0,+0.5,+0.0}, {+0.5,+1.0,+0.0}, {+0.0,+0.5,+0.0},
      {+0.0,+0.0,+0.5}, {+1.0,+0.0,+0.5}, {+1.0,+1.0,+0.5}, {+0.0,+1.0,+0.5},
      {+0.5,+0.0,+1.0}, {+1.0,+0.5,+1.0}, {+0.5,+1.0,+1.0}, {+0.0,+0.5,+1.0},
      {+0.5,+0.5,+0.0}, {+0.5,+0.5,+1.0},
      {+0.0,+0.5,+0.5}, {+1.0,+0.5,+0.5},
      {+0.5,+0.0,+0.5}, {+0.5,+1.0,+0.5},
      {+0.5,+0.5,+0.5}
  };

  int dim = 3;
  int polynomialOrder = 2;

  stk::mesh::MetaData meta(dim);
  stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);

  std::string singleElemMeshSpec = "generated:1x1x1";
  unit_test_utils::fill_and_promote_hex_mesh(singleElemMeshSpec, bulk, polynomialOrder);
  const stk::mesh::PartVector promotedElemParts = sierra::nalu::only_super_elem_parts(meta.get_parts());
  const stk::mesh::Selector promotedElemSelector = stk::mesh::selectUnion(promotedElemParts);
  const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::ELEM_RANK, promotedElemSelector);

  stk::mesh::EntityVector elems;
  stk::mesh::get_selected_entities(promotedElemSelector, buckets, elems);
  ASSERT_EQ(elems.size(), 1u);

  VectorFieldType* coordField = meta.get_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates");
  for (stk::mesh::Entity elem : elems) {
    const stk::mesh::Entity* elemNodeRelations = bulk.begin_nodes(elem);
    for (unsigned k = 0; k < bulk.num_nodes(elem); ++k) {
      const stk::mesh::Entity node = elemNodeRelations[k];
      const double* stkCoordsForNodeK = stk::mesh::field_data(*coordField, node);
      for (int d = 0; d < dim; ++d) {
        EXPECT_NEAR(stkCoordsForNodeK[d],  expectedCoords[k][d], 1.0e-14);
      }
    }
  }

  bool doOutput = false;
  if (doOutput) {
    unit_test_utils::dump_promoted_mesh_file(bulk, polynomialOrder);
  }
}
