#ifndef _UnitTestNonConformalMeshUtil_h_
#define _UnitTestNonConformalMeshUtil_h_

#include <string>
#include <ostream>

#include <SimdInterface.h>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_topology/topology.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <Kokkos_Core.hpp>

#include <master_element/Hex8CVFEM.h>

typedef stk::mesh::Field<double> ScalarFieldType;
typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType;
typedef stk::mesh::Field<double,stk::mesh::Cartesian,stk::mesh::Cartesian> TensorFieldType;
typedef stk::mesh::Field<double, stk::mesh::SimpleArrayTag>  GenericFieldType;
typedef sierra::nalu::DoubleType DoubleType;

namespace stk { namespace mesh { class FieldBase; } }

namespace unit_test_utils {

  struct DGSurface {
    stk::topology::topology_t elemTopo;
    stk::topology::topology_t faceTopo;
    stk::mesh::Part* part;
  };

  // makes a mesh with two different blocks generated through stk::mesh::fixtures
  template<class FixtureA, class FixtureB>
  class NonConformalMeshHelper
  {
  public:
    NonConformalMeshHelper(int nA, int nB);
    void declare_fields_and_commit();
    void create_mesh();

    int NA;
    int NB;
    stk::mesh::MetaData meta_;
    stk::mesh::BulkData bulk_;
    FixtureA fixture_A;
    FixtureB fixture_B;
    std::pair<DGSurface, DGSurface> dgInterface;
  };
}

#endif
