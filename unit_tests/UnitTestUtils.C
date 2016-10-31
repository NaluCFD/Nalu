

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

#include <algorithm>
#include <string>
#include <array>

namespace unit_test_utils {

void fill_mesh_1_elem_per_proc_hex8(stk::mesh::BulkData& bulk)
{
    int nprocs = bulk.parallel_size();
    std::string meshSpec = "generated:1x1x"+std::to_string(nprocs);

    stk::io::StkMeshIoBroker io(bulk.parallel());
    io.set_bulk_data(bulk);
    io.add_mesh_database(meshSpec, stk::io::READ_MESH);
    io.create_input_mesh();
    io.populate_bulk_data();
}
void create_one_reference_hex8_element(stk::mesh::BulkData& bulk)
{
  // Create one -1/2:1/2 reference element for the hex8 topology

   auto& meta = bulk.mesh_meta_data();

   stk::topology hex8 = stk::topology::HEX_8;
   stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", hex8);

   stk::topology quad4 = stk::topology::QUAD_4;
   stk::mesh::PartVector allSurfaces = { &meta.declare_part_with_topology("all_surfaces", quad4) };

   // set a coordinate field
   using vector_field_type = stk::mesh::Field<double, stk::mesh::Cartesian3d>;
   auto& coordField = meta.declare_field<vector_field_type>(stk::topology::NODE_RANK, "coordinates");
   stk::mesh::put_field(coordField, block_1);
   stk::mesh::put_field(coordField, stk::mesh::selectUnion(allSurfaces));
   meta.set_coordinate_field(&coordField);
   meta.commit();

   stk::mesh::EntityIdVector nodeIds(hex8.num_nodes());
   std::iota(nodeIds.begin(), nodeIds.end(), 1);

   bulk.modification_begin();

   for (auto id : nodeIds) {
     bulk.declare_entity(stk::topology::NODE_RANK, id);
   }
   auto elem = stk::mesh::declare_element (bulk, block_1, 1, nodeIds);
   stk::mesh::create_all_sides(bulk, block_1, allSurfaces, false);

   bulk.modification_end();

   // coordinate map between ordinal and position for a -1/2:1/2 hex8 element
   // from exodus standard for P=1 hex element
   std::array<std::array<double, 3>, 8> locMap =
   {{
       {{-0.5,-0.5,-0.5}}, {{+0.5,-0.5,-0.5}}, {{+0.5,+0.5,-0.5}}, {{-0.5,+0.5,-0.5}},
       {{-0.5,-0.5,+0.5}}, {{+0.5,-0.5,+0.5}}, {{+0.5,+0.5,+0.5}}, {{-0.5,+0.5,+0.5}}
   }};

   const auto* nodes = bulk.begin_nodes(elem);
   for (unsigned j = 0; j  < bulk.num_nodes(elem); ++j) {
     for (unsigned i = 0; i < 3; ++i) {
       stk::mesh::field_data(coordField, nodes[j])[i] = locMap.at(j).at(i);
     }
   }
}


void create_one_reference_hex27_element(stk::mesh::BulkData& bulk)
{
  // Create one -1:1 reference element for the hex27 topology

   auto& meta = bulk.mesh_meta_data();

   stk::topology hex27 = stk::topology::HEX_27;
   stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", hex27);

   stk::topology quad9 = stk::topology::QUAD_9;
   stk::mesh::PartVector allSurfaces = { &meta.declare_part_with_topology("all_surfaces", quad9) };

   // set a coordinate field
   using vector_field_type = stk::mesh::Field<double, stk::mesh::Cartesian3d>;
   auto& coordField = meta.declare_field<vector_field_type>(stk::topology::NODE_RANK, "coordinates");
   stk::mesh::put_field(coordField, block_1);
   stk::mesh::put_field(coordField, stk::mesh::selectUnion(allSurfaces));
   meta.set_coordinate_field(&coordField);
   meta.commit();

   stk::mesh::EntityIdVector nodeIds(hex27.num_nodes());
   std::iota(nodeIds.begin(), nodeIds.end(), 1);

   bulk.modification_begin();

   for (auto id : nodeIds) {
     bulk.declare_entity(stk::topology::NODE_RANK, id);
   }
   auto elem = stk::mesh::declare_element (bulk, block_1, 1, nodeIds);
   stk::mesh::create_all_sides(bulk, block_1, allSurfaces, false);

   bulk.modification_end();

   // coordinate map between ordinal and position for a -1:1 hex 27 element
   // from exodus standard for P=2 hex element
   std::array<std::array<double, 3>, 27> locMap =
   {{
       {{-1.0,-1.0,-1.0}}, {{+1.0,-1.0,-1.0}}, {{+1.0,+1.0,-1.0}}, {{-1.0,+1.0,-1.0}},
       {{-1.0,-1.0,+1.0}}, {{+1.0,-1.0,+1.0}}, {{+1.0,+1.0,+1.0}}, {{-1.0,+1.0,+1.0}},
       {{+0.0,-1.0,-1.0}}, {{+1.0,+0.0,-1.0}}, {{+0.0,+1.0,-1.0}}, {{-1.0,+0.0,-1.0}},
       {{-1.0,-1.0,+0.0}}, {{+1.0,-1.0,+0.0}}, {{+1.0,+1.0, 0.0}}, {{-1.0,+1.0,+0.0}},
       {{+0.0,-1.0,+1.0}}, {{+1.0,+0.0,+1.0}}, {{+0.0,+1.0,+1.0}}, {{-1.0,+0.0,+1.0}},
       {{+0.0,+0.0,+0.0}},
       {{+0.0,+0.0,-1.0}}, {{+0.0,+0.0,+1.0}},
       {{-1.0,+0.0,+0.0}}, {{+1.0,+0.0,+0.0}},
       {{+0.0,-1.0,+0.0}}, {{+0.0,+1.0,+0.0}}
   }};

   const auto* nodes = bulk.begin_nodes(elem);
   for (unsigned j = 0; j  < bulk.num_nodes(elem); ++j) {
     for (unsigned i = 0; i < 3; ++i) {
       stk::mesh::field_data(coordField, nodes[j])[i] = locMap.at(j).at(i);
     }
   }
}

}

