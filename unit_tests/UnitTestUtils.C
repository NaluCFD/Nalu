#include <string>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

#include "UnitTestUtils.h"

namespace unit_test_utils {

void fill_mesh_1_elem_per_proc_hex8(stk::mesh::BulkData& bulk)
{
    int nprocs = bulk.parallel_size();
    std::string meshSpec = "generated:1x1x"+std::to_string(nprocs);
    fill_mesh_hex8(bulk, meshSpec);
}

void fill_mesh_hex8(stk::mesh::BulkData& bulk, const std::string& meshSpec)
{
    stk::io::StkMeshIoBroker io(bulk.parallel());
    io.set_bulk_data(bulk);
    io.add_mesh_database(meshSpec, stk::io::READ_MESH);
    io.create_input_mesh();
    io.populate_bulk_data();
}

}// namespace unit_test_utils

