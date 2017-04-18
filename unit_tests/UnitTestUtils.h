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
void fill_and_promote_hex_mesh(const std::string& meshSpec, stk::mesh::BulkData& bulk, int polyOrder);
void dump_promoted_mesh_file(stk::mesh::BulkData& bulk, int polyOrder);

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


template <typename T>
void dump_2d_view(const T& v)
{
  for (unsigned j = 0; j < v.dimension_1(); ++j) {
    for (unsigned i = 0; i < v.dimension_0(); ++i) {
      double vout = (std::abs(v(j,i)) < 1.0e-15) ? 0.0 : v(j,i);
      std::cout << "(" << i << ", " << j << ")" << ", " << std::setw(5) << v.label() << ": " << std::setw(12) << vout;
      if (i != v.dimension_0()-1) { std::cout << ", "; }
    }
    std::cout << std::endl;
  }
  std::cout << "--------------" << std::endl;
}
}

#define EXPECT_VIEW_NEAR_1D(x,y,tol)                 \
{                                                    \
  EXPECT_EQ(x.dimension_0(), y.dimension_0());       \
  for (unsigned i = 0; i < x.dimension_0(); ++i) {   \
    EXPECT_NEAR(x(i), y(i), tol);                    \
  }                                                  \
}

#define EXPECT_VIEW_NEAR_2D(x,y,tol)                  \
{                                                     \
  EXPECT_EQ(x.dimension_0(), y.dimension_0());        \
  for (unsigned j = 0; j < x.dimension_0();++j) {     \
    for (unsigned i = 0; i < x.dimension_1();++i) {   \
      EXPECT_NEAR(x(j,i), y(j,i), tol);               \
    }                                                 \
  }                                                   \
}

#define EXPECT_VIEW_NEAR_3D(x,y,tol)                   \
{                                                      \
  EXPECT_EQ(x.dimension_0(), y.dimension_0());         \
  EXPECT_EQ(x.dimension_1(), y.dimension_1());         \
  EXPECT_EQ(x.dimension_2(), y.dimension_2());         \
  for (unsigned k = 0; k < x.dimension_0();++k) {      \
    for (unsigned j = 0; j < x.dimension_1();++j) {    \
      for (unsigned i = 0; i < x.dimension_2();++i) {  \
        EXPECT_NEAR(x(k,j,i), y(k,j,i), tol);          \
      }                                                \
    }                                                  \
  }                                                    \
}

#define EXPECT_VIEW_NEAR_4D(x,y,tol)                     \
{                                                        \
  EXPECT_EQ(x.dimension_0(), y.dimension_0());           \
  EXPECT_EQ(x.dimension_1(), y.dimension_1());           \
  EXPECT_EQ(x.dimension_2(), y.dimension_2());           \
  EXPECT_EQ(x.dimension_3(), y.dimension_3());           \
  for (unsigned l = 0; l < x.dimension_0(); ++l) {       \
    for (unsigned k = 0; k < x.dimension_1(); ++k) {     \
      for (unsigned j = 0; j < x.dimension_2(); ++j) {   \
        for (unsigned i = 0; i < x.dimension_3(); ++i) { \
          EXPECT_NEAR(x(l,k,j,i), y(l,k,j,i), tol);      \
        }                                                \
      }                                                  \
    }                                                    \
  }                                                      \
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

