#include <gtest/gtest.h>
#include <limits>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <master_element/MasterElement.h>

#include "UnitTestUtils.h"

namespace {

void check_HexSCV_determinant(const stk::mesh::BulkData& bulk)
{
    typedef stk::mesh::Field<double,stk::mesh::Cartesian> CoordFieldType;
    CoordFieldType* coordField = bulk.mesh_meta_data().get_field<CoordFieldType>(stk::topology::NODE_RANK, "coordinates");
    EXPECT_TRUE(coordField != nullptr);

    stk::mesh::EntityVector elems;
    stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, elems);

    const double tolerance = std::numeric_limits<double>::epsilon();

    stk::topology hex8 = stk::topology::HEX_8;

    const unsigned spatialDim = bulk.mesh_meta_data().spatial_dimension();
    std::vector<double> hex8_node_coords(hex8.num_nodes()*spatialDim, 0.0);
    std::vector<double> hex8_scvolumes(hex8.num_nodes(), 0.0);

    sierra::nalu::HexSCV hexSCV;
    double error[1] = {0};
    for(stk::mesh::Entity elem : elems) {
        EXPECT_EQ(hex8, bulk.bucket(elem).topology());

        const stk::mesh::Entity* nodes = bulk.begin_nodes(elem);
        unsigned numNodes = bulk.num_nodes(elem);
        unsigned counter = 0;
        for(unsigned i=0; i<numNodes; ++i) {
            double* nodeCoords = stk::mesh::field_data(*coordField, nodes[i]);
            
            for(unsigned d=0; d<spatialDim; ++d) {
                hex8_node_coords[counter++] = nodeCoords[d];
            }
        }

        hexSCV.determinant(1, hex8_node_coords.data(), hex8_scvolumes.data(), error);

        EXPECT_EQ(0, error[0]);

        for(unsigned i=0; i<hex8.num_nodes(); ++i) {
           EXPECT_NEAR(0.125, hex8_scvolumes[i], tolerance);
        }
    }
}

}//namespace

TEST(HexSCV, determinant)
{
    stk::ParallelMachine comm = MPI_COMM_WORLD;

    unsigned spatialDimension = 3;
    stk::mesh::MetaData meta(spatialDimension);
    stk::mesh::BulkData bulk(meta, comm);

    unit_test_utils::fill_mesh_1_elem_per_proc_hex8(bulk);

    check_HexSCV_determinant(bulk);
}

