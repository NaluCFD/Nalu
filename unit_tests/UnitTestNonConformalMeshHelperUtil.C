#include "UnitTestNonConformalMeshHelperUtil.h"

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

namespace unit_test_utils {

  template<class FixtureA, class FixtureB>
  NonConformalMeshHelper<FixtureA, FixtureB>::NonConformalMeshHelper(int nA, int nB)
: NA(nA), NB(nB), meta_(3), bulk_(meta_, MPI_COMM_WORLD),
  fixture_A(meta_, bulk_, NA, NA, NA, 1, 1),
  fixture_B(meta_, bulk_, NB, NB, NB, fixture_A.num_nodes() + 1, fixture_A.num_elements() + 1)
  {
    auto& block_1 = meta_.declare_part_with_topology("block_1", fixture_A.m_elem_topology);
    fixture_A.m_elem_parts.push_back(&block_1);
    stk::io::put_io_part_attribute(block_1);

    auto& block_2 = meta_.declare_part_with_topology("block_2", fixture_B.m_elem_topology);
    fixture_B.m_elem_parts.push_back(&block_2);
    stk::io::put_io_part_attribute(block_2);

    meta_.declare_part("all_block_boundary_surface");

    dgInterface.first.elemTopo = fixture_A.m_elem_topology;
    dgInterface.first.faceTopo = fixture_A.m_face_topology;
    dgInterface.first.part = &meta_.declare_part_with_topology("surface_1", fixture_A.m_face_topology);

    dgInterface.second.elemTopo = fixture_B.m_elem_topology;
    dgInterface.second.faceTopo = fixture_B.m_face_topology;
    dgInterface.second.part = &meta_.declare_part_with_topology("surface_2", fixture_B.m_face_topology);

    stk::io::put_io_part_attribute(*dgInterface.first.part);
    stk::io::put_io_part_attribute(*dgInterface.second.part);

    declare_fields_and_commit();
    create_mesh();
  }

  template <class FixtureA, class FixtureB>
  void NonConformalMeshHelper<FixtureA, FixtureB>::declare_fields_and_commit()
  {
    auto& velField = meta_.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
    stk::mesh::put_field(velField, meta_.universal_part(), 3);

    auto& pressField = meta_.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure");
    stk::mesh::put_field(pressField, meta_.universal_part(), 1);

    meta_.commit();
  }

  template<class FixtureA, class FixtureB>
  void NonConformalMeshHelper<FixtureA, FixtureB>::create_mesh()
  {
    using BlockMapping = stk::mesh::fixtures::FixedCartesianCoordinateMapping;
    fixture_A.generate_mesh(BlockMapping(NA, NA, NA, 1.0, 1.0,1.0));
    fixture_B.generate_mesh(BlockMapping(NB, NB, NB, 1.0, 1.0, 1.0));
    const VectorFieldType& coordField = *static_cast<const VectorFieldType*>(meta_.coordinate_field());

    stk::mesh::Selector block1_selector = *meta_.get_part("block_1");
    stk::mesh::Selector block2_selector = *meta_.get_part("block_2");
    stk::mesh::BucketVector block1_node_buckets = bulk_.get_buckets(stk::topology::NODE_RANK, block1_selector);
    stk::mesh::BucketVector block2_node_buckets = bulk_.get_buckets(stk::topology::NODE_RANK, block2_selector);

    // Move fitxture B over in the x-direction by 1.0, then rotate 90 degrees so processor boundaries are not aligned
    for(const auto* ib : block2_node_buckets)  {
      const stk::mesh::Bucket& b = *ib;
      auto length = b.size();

      for(size_t k = 0u; k < length; ++k) {
        double* coords = stk::mesh::field_data(coordField , b[k]);

        // x' = 1 + x, y '= z, z' = 1-y
        coords[0] += 1.0;
        std::swap(coords[1], coords[2]);
        coords[2] = 1 - coords[2];
      }
    }

    // create sides of both blocks
    stk::mesh::create_exposed_block_boundary_sides(bulk_,meta_.universal_part(), {meta_.get_part("all_block_boundary_surface")});

    // add faces to the left/right (surface_1/surface_2) of the interface, knowing that the interface is at x = 0
    stk::mesh::EntityVector entities;
    std::vector<stk::mesh::PartVector> add_parts;
    std::vector<stk::mesh::PartVector> remove_parts;
    stk::mesh::BucketVector side_buckets = bulk_.get_buckets(meta_.side_rank(), meta_.locally_owned_part());
    for (const auto* ib : side_buckets) {
      const stk::mesh::Bucket& b = *ib;
      auto length = b.size();

      for (size_t k = 0u; k < length; ++k) {
        stk::mesh::Entity side = b[k];
        const auto* sideNodes = bulk_.begin_nodes(side);
        unsigned nSideNodes = bulk_.num_nodes(side);
        double x_mean = 0.0;
        for (unsigned n = 0; n < nSideNodes; ++n) {
          double* cPtr = stk::mesh::field_data(coordField, sideNodes[n]);
          x_mean += cPtr[0];
        }
        x_mean /= nSideNodes;

        stk::mesh::Entity elem = b.begin_elements(k)[0];
        const auto* elemNodes = bulk_.begin_nodes(elem);
        unsigned nElemNodes = bulk_.num_nodes(elem);

        double xe_mean = 0.0;
        for (unsigned n = 0; n < nElemNodes; ++n) {
          double* cPtr = stk::mesh::field_data(coordField, elemNodes[n]);
          xe_mean += cPtr[0];
        }
        xe_mean /= nElemNodes;

        if (std::abs(x_mean - 1.0) < 1e-8) {
          entities.push_back(side);
          remove_parts.push_back({ });
          if (xe_mean < 1.0) {
            add_parts.push_back({dgInterface.first.part});
          }
          else {
            add_parts.push_back({dgInterface.second.part});
          }
        }
      }
    }
    bulk_.batch_change_entity_parts(entities, add_parts, remove_parts);
  }

} // namespace unit_test_utils

// ETI all 3d combos
template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::HexFixture, stk::mesh::fixtures::HexFixture> ;
template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::HexFixture, stk::mesh::fixtures::TetFixture> ;
template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::HexFixture, stk::mesh::fixtures::PyramidFixture> ;
template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::HexFixture, stk::mesh::fixtures::WedgeFixture> ;

template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::PyramidFixture, stk::mesh::fixtures::HexFixture> ;
template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::PyramidFixture, stk::mesh::fixtures::TetFixture> ;
template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::PyramidFixture, stk::mesh::fixtures::WedgeFixture> ;
template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::PyramidFixture, stk::mesh::fixtures::PyramidFixture> ;

template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::TetFixture, stk::mesh::fixtures::HexFixture> ;
template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::TetFixture, stk::mesh::fixtures::PyramidFixture> ;
template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::TetFixture, stk::mesh::fixtures::WedgeFixture> ;
template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::TetFixture, stk::mesh::fixtures::TetFixture> ;

template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::WedgeFixture, stk::mesh::fixtures::HexFixture> ;
template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::WedgeFixture, stk::mesh::fixtures::PyramidFixture> ;
template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::WedgeFixture, stk::mesh::fixtures::TetFixture> ;
template class unit_test_utils::NonConformalMeshHelper<stk::mesh::fixtures::WedgeFixture, stk::mesh::fixtures::WedgeFixture> ;

