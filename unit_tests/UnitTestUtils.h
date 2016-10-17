
#include <stk_mesh/base/BulkData.hpp>

namespace unit_test_utils {

void fill_mesh_1_elem_per_proc_hex8(stk::mesh::BulkData& bulk);
void fill_mesh_hex8(stk::mesh::BulkData& bulk, const std::string& meshSpec);

}

