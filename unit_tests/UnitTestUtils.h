#ifndef _UnitTestUtils_h_
#define _UnitTestUtils_h_

#include <string>
#include <ostream>

#include <SimdInterface.h>
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
typedef stk::mesh::Field<double, stk::mesh::SimpleArrayTag>  GenericFieldType;
typedef sierra::nalu::DoubleType DoubleType;

namespace stk { namespace mesh { class FieldBase; } }

namespace unit_test_utils {

void fill_mesh_1_elem_per_proc_hex8(stk::mesh::BulkData& bulk);
void fill_hex8_mesh(const std::string& meshSpec, stk::mesh::BulkData& bulk);
void perturb_coord_hex_8(stk::mesh::BulkData& bulk, double perturbationSize = 0.125);


void fill_and_promote_hex_mesh(const std::string& meshSpec, stk::mesh::BulkData& bulk, int polyOrder);

void dump_mesh(stk::mesh::BulkData& bulk, std::vector<stk::mesh::FieldBase*> fields);

void dump_promoted_mesh_file(stk::mesh::BulkData& bulk, int polyOrder);

std::ostream& nalu_out();

stk::mesh::Entity create_one_reference_element(stk::mesh::BulkData& bulk, stk::topology topo);
stk::mesh::Entity create_one_perturbed_element(stk::mesh::BulkData& bulk, stk::topology topo);

double quadratic(double a, const double* b, const double* H, const double* x);

double vector_norm(const std::vector<double> & vec, const stk::ParallelMachine& comm = MPI_COMM_WORLD);

double initialize_quadratic_scalar_field(const stk::mesh::BulkData& bulk,
                                      const VectorFieldType& coordField,
                                      const ScalarFieldType& qField);

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
      scalarQ(&meta.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "scalarQ")),
      diffFluxCoeff(&meta.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "diffFluxCoeff")),
      partVec(),
      coordField(nullptr),
      exactLaplacian(0.0)
    { 
      stk::mesh::put_field(*elemCentroidField, meta.universal_part(), spatialDimension, (double*)nullptr); 
      double one = 1.0;
      double zero = 0.0;
      stk::mesh::put_field(*nodalPressureField, meta.universal_part(), 1, &one);
      stk::mesh::put_field(*discreteLaplacianOfPressure, meta.universal_part(), 1, 0.0);
      stk::mesh::put_field(*discreteLaplacianOfPressure, meta.universal_part(), 1, &zero);
      stk::mesh::put_field(*scalarQ, meta.universal_part(), 1, &zero);
      stk::mesh::put_field(*diffFluxCoeff, meta.universal_part(), 1, &zero);
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
        stk::mesh::field_fill(0.1, *scalarQ);
        stk::mesh::field_fill(0.2, *diffFluxCoeff);
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
    ScalarFieldType* scalarQ;
    ScalarFieldType* diffFluxCoeff;
    stk::mesh::PartVector partVec;
    const VectorFieldType* coordField;
    double exactLaplacian; 
};

class Hex8MeshWithNSOFields : public Hex8Mesh
{
protected:
    Hex8MeshWithNSOFields()
    : Hex8Mesh(),
      massFlowRate(&meta.declare_field<VectorFieldType>(stk::topology::ELEM_RANK, "massFlowRate")),   
      Gju(&meta.declare_field<GenericFieldType>(stk::topology::NODE_RANK, "Gju", 1/*num-states*/)), 
      velocity(&meta.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity", 3/*num-states*/)), 
      density(&meta.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "density", 3/*num-states*/)), 
      viscosity(&meta.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "viscosity")),
      pressure(&meta.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "pressure"))
    {
      double one = 1.0;
      sierra::nalu::HexSCS hex8SCS;
      stk::mesh::put_field(*massFlowRate, meta.universal_part(), hex8SCS.numIntPoints_, &one);
      stk::mesh::put_field(*Gju, meta.universal_part(), 3, &one);
      stk::mesh::put_field(*velocity, meta.universal_part(), 3, &one);
      stk::mesh::put_field(*density, meta.universal_part(), 1, &one);
      stk::mesh::put_field(*viscosity, meta.universal_part(), 1, &one);
      stk::mesh::put_field(*pressure, meta.universal_part(), 1, &one);
    }

    VectorFieldType* massFlowRate;
    GenericFieldType* Gju;
    VectorFieldType* velocity;
    ScalarFieldType* density;
    ScalarFieldType* viscosity;
    ScalarFieldType* pressure;
};

#endif

#endif

