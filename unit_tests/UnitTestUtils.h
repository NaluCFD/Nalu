#include <string>
#include <ostream>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_topology/topology.hpp>

namespace unit_test_utils {

void fill_mesh_1_elem_per_proc_hex8(stk::mesh::BulkData& bulk);
void fill_hex8_mesh(const std::string& meshSpec, stk::mesh::BulkData& bulk);
std::ostream& nalu_out();
void create_one_reference_element(stk::mesh::BulkData& bulk, stk::topology topo);


}

