#include <string>

#include <stk_mesh/base/BulkData.hpp>

namespace unit_test_utils {

void fill_mesh_1_elem_per_proc_hex8(stk::mesh::BulkData& bulk);
void fill_hex8_mesh(const std::string& meshSpec, stk::mesh::BulkData& bulk);

void create_one_reference_hex8_element(stk::mesh::BulkData& bulk);
void create_one_reference_hex27_element(stk::mesh::BulkData& bulk);

}

