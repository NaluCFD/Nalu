#include <gtest/gtest.h>
#include <NaluEnv.h>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_topology/topology.hpp>

#include "UnitTestUtils.h"

#include <algorithm>
#include <string>
#include <array>

namespace unit_test_utils {

void fill_mesh_1_elem_per_proc_hex8(stk::mesh::BulkData& bulk)
{
    int nprocs = bulk.parallel_size();
    std::string meshSpec = "generated:1x1x"+std::to_string(nprocs);
    fill_hex8_mesh(meshSpec, bulk);
}

void fill_hex8_mesh(const std::string& meshSpec, stk::mesh::BulkData& bulk)
{
    stk::io::StkMeshIoBroker io(bulk.parallel());
    io.set_bulk_data(bulk);
    io.add_mesh_database(meshSpec, stk::io::READ_MESH);
    io.create_input_mesh();
    io.populate_bulk_data();
}

std::ostream& nalu_out()
{
  return sierra::nalu::NaluEnv::self().naluOutputP0();
}

void create_one_element(
  stk::mesh::BulkData& bulk,
  stk::topology topo, const
  std::vector<std::vector<double>>& nodeLocations)
{
  // create just one element

   auto& meta = bulk.mesh_meta_data();
   stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", topo);
   stk::mesh::PartVector allSurfaces = { &meta.declare_part("all_surfaces", meta.side_rank()) };

   // set a coordinate field
   using vector_field_type = stk::mesh::Field<double, stk::mesh::Cartesian3d>;
   auto& coordField = meta.declare_field<vector_field_type>(stk::topology::NODE_RANK, "coordinates");
   stk::mesh::put_field(coordField, block_1);
   stk::mesh::put_field(coordField, stk::mesh::selectUnion(allSurfaces));
   meta.set_coordinate_field(&coordField);
   meta.commit();

   stk::mesh::EntityIdVector nodeIds(topo.num_nodes());
   std::iota(nodeIds.begin(), nodeIds.end(), 1);

   bulk.modification_begin();

   for (auto id : nodeIds) {
     bulk.declare_entity(stk::topology::NODE_RANK, id);
   }
   auto elem = stk::mesh::declare_element (bulk, block_1, 1, nodeIds);
   stk::mesh::create_all_sides(bulk, block_1, allSurfaces, false);

   bulk.modification_end();

   const auto* nodes = bulk.begin_nodes(elem);
   for (unsigned j = 0; j  < bulk.num_nodes(elem); ++j) {
     for (unsigned i = 0; i < topo.dimension(); ++i) {
       stk::mesh::field_data(coordField, nodes[j])[i] = nodeLocations.at(j).at(i);
     }
   }
}

void create_one_reference_quad4_element(stk::mesh::BulkData& bulk)
{
  std::vector<std::vector<double>> nodeLocations =
  {
      {-0.5,-0.5}, {+0.5,-0.5},
      {+0.5,+0.5}, {-0.5,+0.5}
  };
  create_one_element(bulk, stk::topology::QUADRILATERAL_4_2D, nodeLocations);
}

void create_one_reference_quad9_element(stk::mesh::BulkData& bulk)
{
  std::vector<std::vector<double>> nodeLocations =
  {
      {-1.0,-1.0}, {+1.0,-1.0},
      {+1.0,+1.0}, {-1.0,+1.0},
      {0.0, -1.0}, {+1.0, 0.0}, {0.0, +1.0}, {-1.0, 0.0},
      {0.0, 0.0}
  };
  create_one_element(bulk, stk::topology::QUADRILATERAL_9_2D, nodeLocations);
}

void create_one_reference_tri3_element(stk::mesh::BulkData& bulk)
{
  std::vector<std::vector<double>> nodeLocations =
  {
      {-0.5,-0.5}, {+0.5,-0.5}, {-0.5,+0.5}
  };
  create_one_element(bulk, stk::topology::TRIANGLE_3_2D, nodeLocations);
}

void create_one_reference_tet4_element(stk::mesh::BulkData& bulk)
{
   std::vector<std::vector<double>> nodeLocations =
   {
       {-0.5,-0.5,-0.5}, {+0.5,-0.5,-0.5}, {-0.5,+0.5,-0.5}, {-0.5,-0.5,+0.5}
   };
   create_one_element(bulk, stk::topology::TET_4, nodeLocations);
}

void create_one_reference_hex8_element(stk::mesh::BulkData& bulk)
{
   std::vector<std::vector<double>> nodeLocations =
   {
       {-0.5,-0.5,-0.5}, {+0.5,-0.5,-0.5}, {+0.5,+0.5,-0.5}, {-0.5,+0.5,-0.5},
       {-0.5,-0.5,+0.5}, {+0.5,-0.5,+0.5}, {+0.5,+0.5,+0.5}, {-0.5,+0.5,+0.5}
   };
   create_one_element(bulk, stk::topology::HEX_8, nodeLocations);
}

void create_one_reference_hex27_element(stk::mesh::BulkData& bulk)
{
   std::vector<std::vector<double>> nodeLocations =
   {
       {-1.0,-1.0,-1.0}, {+1.0,-1.0,-1.0}, {+1.0,+1.0,-1.0}, {-1.0,+1.0,-1.0},
       {-1.0,-1.0,+1.0}, {+1.0,-1.0,+1.0}, {+1.0,+1.0,+1.0}, {-1.0,+1.0,+1.0},
       {+0.0,-1.0,-1.0}, {+1.0,+0.0,-1.0}, {+0.0,+1.0,-1.0}, {-1.0,+0.0,-1.0},
       {-1.0,-1.0,+0.0}, {+1.0,-1.0,+0.0}, {+1.0,+1.0, 0.0}, {-1.0,+1.0,+0.0},
       {+0.0,-1.0,+1.0}, {+1.0,+0.0,+1.0}, {+0.0,+1.0,+1.0}, {-1.0,+0.0,+1.0},
       {+0.0,+0.0,+0.0},
       {+0.0,+0.0,-1.0}, {+0.0,+0.0,+1.0},
       {-1.0,+0.0,+0.0}, {+1.0,+0.0,+0.0},
       {+0.0,-1.0,+0.0}, {+0.0,+1.0,+0.0}
   };
   create_one_element(bulk, stk::topology::HEX_27, nodeLocations);
}

void create_one_reference_pyramid5_element(stk::mesh::BulkData& bulk)
{
  std::vector<std::vector<double>> nodeLocations =
  {
      {-1.0, -1.0, +0.0}, {+1.0, -1.0, +0.0}, {+1.0, +1.0, +0.0}, {-1.0, +1.0, +0.0},
      {0.0, 0.0, +1.0}
  };
   create_one_element(bulk, stk::topology::PYRAMID_5, nodeLocations);
}

void create_one_reference_wedge6_element(stk::mesh::BulkData& bulk)
{
   std::vector<std::vector<double>> nodeLocations =
   {
       {-0.5,-0.5,-1.0}, {+0.5,-0.5,-1.0}, {-0.5,+0.5,-1.0},
       {-0.5,-0.5,+1.0}, {+0.5,-0.5,+1.0}, {-0.5,+0.5,+1.0}
   };
   create_one_element(bulk, stk::topology::WEDGE_6, nodeLocations);
}

void create_one_reference_element(stk::mesh::BulkData& bulk, stk::topology topo)
{
  switch (topo.value())
  {
    case stk::topology::TRIANGLE_3_2D:
      create_one_reference_tri3_element(bulk);
      break;

    case stk::topology::QUADRILATERAL_4_2D:
      create_one_reference_quad4_element(bulk);
      break;

    case stk::topology::QUADRILATERAL_9_2D:
      create_one_reference_quad9_element(bulk);
      break;

    case stk::topology::TETRAHEDRON_4:
      create_one_reference_tet4_element(bulk);
      break;

    case stk::topology::PYRAMID_5:
      create_one_reference_pyramid5_element(bulk);
      break;

    case stk::topology::WEDGE_6:
      create_one_reference_wedge6_element(bulk);
      break;

    case stk::topology::HEXAHEDRON_8:
      create_one_reference_hex8_element(bulk);
      break;

    case stk::topology::HEXAHEDRON_27:
      create_one_reference_hex27_element(bulk);
      break;

    default: FAIL() << " Element type " + topo.name() + " not implemented";
  }
}



}

