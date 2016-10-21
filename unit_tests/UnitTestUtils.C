#include <string>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

#include <algorithm>

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

void create_one_hex27_element(stk::mesh::BulkData& bulk)
{
  // no coords yet
   auto& meta = bulk.mesh_meta_data();

   stk::topology hex27 = stk::topology::HEX_27;
   stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", hex27);

   stk::mesh::PartVector surfaces(6);
   stk::topology quad9 = stk::topology::QUAD_9;
   surfaces[0] = &meta.declare_part_with_topology("surface_1", quad9);
   surfaces[1] = &meta.declare_part_with_topology("surface_2", quad9);
   surfaces[2] = &meta.declare_part_with_topology("surface_3", quad9);
   surfaces[3] = &meta.declare_part_with_topology("surface_4", quad9);
   surfaces[4] = &meta.declare_part_with_topology("surface_5", quad9);
   surfaces[5] = &meta.declare_part_with_topology("surface_6", quad9);
   meta.commit();

   stk::mesh::EntityIdVector nodeIds(27);
   std::iota(nodeIds.begin(), nodeIds.end(), 1);

   std::vector<unsigned> sideNodeOrdinals(quad9.num_nodes());

   bulk.modification_begin();
   {
     for (auto id : nodeIds) {
       bulk.declare_entity(stk::topology::NODE_RANK, id);
     }
     stk::mesh::Entity elem = stk::mesh::declare_element(bulk, block_1, 1, nodeIds);
     const stk::mesh::Entity* elem_node_rels = bulk.begin_nodes(elem);

     for (unsigned j = 0; j < hex27.num_faces(); ++j) {
       stk::mesh::Entity face = bulk.declare_solo_side(j+1, {surfaces[j]});
       hex27.side_node_ordinals(j, sideNodeOrdinals.begin());
       for (unsigned i = 0; i < quad9.num_nodes(); ++i) {
         const stk::mesh::Entity node = elem_node_rels[sideNodeOrdinals[i]];
         bulk.declare_relation(face, node, i);
       }
       bulk.declare_relation(elem, face, j);
     }
   }
   bulk.modification_end();
}

}

