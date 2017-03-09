#ifndef _UnitTestUtils_h_
#define _UnitTestUtils_h_

#include <string>
#include <ostream>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <Kokkos_Core.hpp>

#include <master_element/MasterElement.h>

typedef stk::mesh::Field<double> ScalarFieldType;
typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType;
typedef stk::mesh::Field<double,stk::mesh::Cartesian,stk::mesh::Cartesian> TensorFieldType;

namespace unit_test_utils {

void fill_mesh_1_elem_per_proc_hex8(stk::mesh::BulkData& bulk);
void fill_hex8_mesh(const std::string& meshSpec, stk::mesh::BulkData& bulk);
std::ostream& nalu_out();

stk::mesh::Entity create_one_reference_element(stk::mesh::BulkData& bulk, stk::topology topo);

double quadratic(double a, const double* b, const double* H, const double* x);

double initialize_quadratic_scalar_field(const stk::mesh::BulkData& bulk,
                                      const VectorFieldType& coordField,
                                      const ScalarFieldType& qField);

sierra::nalu::MasterElement *
get_surface_master_element(const stk::topology & theTopo);

sierra::nalu::MasterElement *
get_volume_master_element(const stk::topology & theTopo);


}

const double tol = 1.e-10;

#ifndef KOKKOS_HAVE_CUDA
class Hex8Mesh : public ::testing::Test
{
protected:
    Hex8Mesh()
    : comm(MPI_COMM_WORLD), spatialDimension(3),
      meta(spatialDimension), bulk(meta, comm),
      topo(stk::topology::HEX_8),
      elemCentroidField(&meta.declare_field<VectorFieldType>(stk::topology::ELEM_RANK, "elemCentroid")),   
      nodalPressureField(&meta.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "nodalPressure")), 
      discreteLaplacianOfPressure(&meta.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "discreteLaplacian")),
      partVec(),
      coordField(nullptr),
      exactLaplacian(0.0)
    { 
      stk::mesh::put_field(*elemCentroidField, meta.universal_part(), spatialDimension, (double*)nullptr); 
      double one = 1.0;
      stk::mesh::put_field(*nodalPressureField, meta.universal_part(), 1, &one);
      stk::mesh::put_field(*discreteLaplacianOfPressure, meta.universal_part(), 1, 0.0);
    }

    ~Hex8Mesh() {}

    void fill_mesh(const std::string& meshSpec = "generated:20x20x20")
    { 
      unit_test_utils::fill_hex8_mesh(meshSpec, bulk);
    }

    void fill_mesh_and_initialize_test_fields(const std::string& meshSpec = "generated:20x20x20")
    {   
        fill_mesh(meshSpec);

        partVec = {meta.get_part("block_1")};

        coordField = static_cast<const VectorFieldType*>(meta.coordinate_field());
        EXPECT_TRUE(coordField != nullptr);

        exactLaplacian = unit_test_utils::initialize_quadratic_scalar_field(bulk, *coordField, *nodalPressureField);
        stk::mesh::field_fill(0.0, *discreteLaplacianOfPressure);
    }

    void check_discrete_laplacian(double exactLaplacian);

    stk::ParallelMachine comm;
    unsigned spatialDimension;
    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk;
    stk::topology topo; 
    VectorFieldType* elemCentroidField;
    ScalarFieldType* nodalPressureField;
    ScalarFieldType* discreteLaplacianOfPressure;
    stk::mesh::PartVector partVec;
    const VectorFieldType* coordField;
    double exactLaplacian; 
};

#endif

#endif

