#include <gtest/gtest.h>
#include <NaluEnv.h>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_topology/topology.hpp>

#include "UnitTestUtils.h"
#include "UnitTestKokkosUtils.h"

#include <algorithm>
#include <string>
#include <array>
#include <random>

namespace unit_test_utils {

void fill_mesh_1_elem_per_proc_hex8(stk::mesh::BulkData& bulk)
{
    int nprocs = bulk.parallel_size();
    std::string meshSpec = "generated:1x1x"+std::to_string(nprocs);
    fill_hex8_mesh(meshSpec, bulk);
}

void fill_hex8_mesh(const std::string& meshSpec, stk::mesh::BulkData& bulk)
{
    stk::io::StkMeshIoBroker io(bulk.parallel());
    io.set_bulk_data(bulk);
    io.add_mesh_database(meshSpec, stk::io::READ_MESH);
    io.create_input_mesh();
    io.populate_bulk_data();
}

std::ostream& nalu_out()
{
  return sierra::nalu::NaluEnv::self().naluOutputP0();
}

stk::mesh::Entity create_one_element(
  stk::mesh::BulkData& bulk,
  stk::topology topo, const
  std::vector<std::vector<double>>& nodeLocations)
{
  // create just one element

   auto& meta = bulk.mesh_meta_data();
   stk::mesh::Part& block_1 = meta.declare_part_with_topology("block_1", topo);
   stk::mesh::PartVector allSurfaces = { &meta.declare_part("all_surfaces", meta.side_rank()) };

   // set a coordinate field
   using vector_field_type = stk::mesh::Field<double, stk::mesh::Cartesian3d>;
   auto& coordField = meta.declare_field<vector_field_type>(stk::topology::NODE_RANK, "coordinates");
   stk::mesh::put_field(coordField, block_1);
   stk::mesh::put_field(coordField, stk::mesh::selectUnion(allSurfaces));
   meta.set_coordinate_field(&coordField);
   meta.commit();

   stk::mesh::EntityIdVector nodeIds(topo.num_nodes());
   std::iota(nodeIds.begin(), nodeIds.end(), 1);

   bulk.modification_begin();

   for (auto id : nodeIds) {
     bulk.declare_entity(stk::topology::NODE_RANK, id, {});
   }
   auto elem = stk::mesh::declare_element (bulk, block_1, 1, nodeIds);
   stk::mesh::create_all_sides(bulk, block_1, allSurfaces, false);

   bulk.modification_end();

   const auto* nodes = bulk.begin_nodes(elem);
   for (unsigned j = 0; j  < bulk.num_nodes(elem); ++j) {
     for (unsigned i = 0; i < topo.dimension(); ++i) {
       stk::mesh::field_data(coordField, nodes[j])[i] = nodeLocations.at(j).at(i);
     }
   }

   return elem;
}

stk::mesh::Entity create_one_reference_quad4_element(stk::mesh::BulkData& bulk)
{
  std::vector<std::vector<double>> nodeLocations =
  {
      {-0.5,-0.5}, {+0.5,-0.5},
      {+0.5,+0.5}, {-0.5,+0.5}
  };
  return create_one_element(bulk, stk::topology::QUADRILATERAL_4_2D, nodeLocations);
}

stk::mesh::Entity create_one_reference_quad9_element(stk::mesh::BulkData& bulk)
{
  std::vector<std::vector<double>> nodeLocations =
  {
      {-1.0,-1.0}, {+1.0,-1.0},
      {+1.0,+1.0}, {-1.0,+1.0},
      {0.0, -1.0}, {+1.0, 0.0}, {0.0, +1.0}, {-1.0, 0.0},
      {0.0, 0.0}
  };
  return create_one_element(bulk, stk::topology::QUADRILATERAL_9_2D, nodeLocations);
}

stk::mesh::Entity create_one_reference_tri3_element(stk::mesh::BulkData& bulk)
{
  std::vector<std::vector<double>> nodeLocations =
  {
      {0.0,0.0}, {1.0,0}, {0.0,1.0}
  };
  return create_one_element(bulk, stk::topology::TRIANGLE_3_2D, nodeLocations);
}

stk::mesh::Entity create_one_reference_tet4_element(stk::mesh::BulkData& bulk)
{
  std::vector<std::vector<double>> nodeLocations =
  {
      {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}
  };
   return create_one_element(bulk, stk::topology::TET_4, nodeLocations);
}

stk::mesh::Entity create_one_reference_hex8_element(stk::mesh::BulkData& bulk)
{
   std::vector<std::vector<double>> nodeLocations =
   {
       {-0.5,-0.5,-0.5}, {+0.5,-0.5,-0.5}, {+0.5,+0.5,-0.5}, {-0.5,+0.5,-0.5},
       {-0.5,-0.5,+0.5}, {+0.5,-0.5,+0.5}, {+0.5,+0.5,+0.5}, {-0.5,+0.5,+0.5}
   };
   return create_one_element(bulk, stk::topology::HEX_8, nodeLocations);
}

stk::mesh::Entity create_one_reference_hex27_element(stk::mesh::BulkData& bulk)
{
   std::vector<std::vector<double>> nodeLocations =
   {
       {-1.0,-1.0,-1.0}, {+1.0,-1.0,-1.0}, {+1.0,+1.0,-1.0}, {-1.0,+1.0,-1.0},
       {-1.0,-1.0,+1.0}, {+1.0,-1.0,+1.0}, {+1.0,+1.0,+1.0}, {-1.0,+1.0,+1.0},
       {+0.0,-1.0,-1.0}, {+1.0,+0.0,-1.0}, {+0.0,+1.0,-1.0}, {-1.0,+0.0,-1.0},
       {-1.0,-1.0,+0.0}, {+1.0,-1.0,+0.0}, {+1.0,+1.0, 0.0}, {-1.0,+1.0,+0.0},
       {+0.0,-1.0,+1.0}, {+1.0,+0.0,+1.0}, {+0.0,+1.0,+1.0}, {-1.0,+0.0,+1.0},
       {+0.0,+0.0,+0.0},
       {+0.0,+0.0,-1.0}, {+0.0,+0.0,+1.0},
       {-1.0,+0.0,+0.0}, {+1.0,+0.0,+0.0},
       {+0.0,-1.0,+0.0}, {+0.0,+1.0,+0.0}
   };
   return create_one_element(bulk, stk::topology::HEX_27, nodeLocations);
}

stk::mesh::Entity create_one_reference_pyramid5_element(stk::mesh::BulkData& bulk)
{
  std::vector<std::vector<double>> nodeLocations =
  {
      {-1.0, -1.0, +0.0}, {+1.0, -1.0, +0.0}, {+1.0, +1.0, +0.0}, {-1.0, +1.0, +0.0},
      {0.0, 0.0, +1.0}
  };
  return create_one_element(bulk, stk::topology::PYRAMID_5, nodeLocations);
}

stk::mesh::Entity create_one_reference_wedge6_element(stk::mesh::BulkData& bulk)
{
  std::vector<std::vector<double>> nodeLocations =
  {
      {0.0,0.0,-1.0}, {+1.0,0.0,-1.0}, {0.0,+1.0,-1.0},
      {0.0,0.0,+1.0}, {+1.0,0.0,+1.0}, {0.0,+1.0,+1.0}
  };
  return create_one_element(bulk, stk::topology::WEDGE_6, nodeLocations);
}

stk::mesh::Entity create_one_reference_element(stk::mesh::BulkData& bulk, stk::topology topo)
{
  switch (topo.value())
  {
    case stk::topology::TRIANGLE_3_2D:
      return create_one_reference_tri3_element(bulk);

    case stk::topology::QUADRILATERAL_4_2D:
      return create_one_reference_quad4_element(bulk);

    case stk::topology::QUADRILATERAL_9_2D:
      return create_one_reference_quad9_element(bulk);

    case stk::topology::TETRAHEDRON_4:
      return create_one_reference_tet4_element(bulk);

    case stk::topology::PYRAMID_5:
      return create_one_reference_pyramid5_element(bulk);

    case stk::topology::WEDGE_6:
      return create_one_reference_wedge6_element(bulk);

    case stk::topology::HEXAHEDRON_8:
      return create_one_reference_hex8_element(bulk);

    case stk::topology::HEXAHEDRON_27:
      return create_one_reference_hex27_element(bulk);

    default:
      EXPECT_TRUE(false) << " Element type " + topo.name() + " not implemented";
      return stk::mesh::Entity(0u);
  }
}

double linear(double a, const double* b, const double* x) {
  return (a + b[0]*x[0] + b[1]*x[1] + b[2]*x[2]);
}

double quadratic(double a, const double* b, const double* H, const double* x)
{
  double quad = x[0] * (H[0]*x[0] + H[1]*x[1] + H[2]*x[2])
              + x[1] * (H[3]*x[0] + H[4]*x[1] + H[5]*x[2])
              + x[2] * (H[6]*x[0] + H[7]*x[1] + H[8]*x[2]);

  return (linear(a,b,x) + 0.5*quad);
}

sierra::nalu::MasterElement *
get_surface_master_element(const stk::topology & theTopo)
{
  sierra::nalu::MasterElement *theElem = NULL;

  static std::map<stk::topology, sierra::nalu::MasterElement*> s_topo_masterelem_map;

  std::map<stk::topology, sierra::nalu::MasterElement *>::iterator it =
    s_topo_masterelem_map.find(theTopo);
  if ( it == s_topo_masterelem_map.end() ) {
    theElem = sierra::nalu::MasterElement::create_surface_master_element(theTopo);
    ThrowRequire(theElem != nullptr);

    s_topo_masterelem_map[theTopo] = theElem;
  }
  else {
    theElem = it->second;
  }

  return theElem;
}

sierra::nalu::MasterElement *
get_volume_master_element(const stk::topology & theTopo)
{
  sierra::nalu::MasterElement *theElem = NULL;

  static std::map<stk::topology, sierra::nalu::MasterElement*> s_topov_masterelem_map;

  std::map<stk::topology, sierra::nalu::MasterElement *>::iterator it =
    s_topov_masterelem_map.find(theTopo);
  if ( it == s_topov_masterelem_map.end() ) {
    theElem = sierra::nalu::MasterElement::create_volume_master_element(theTopo);
    ThrowRequire(theElem != nullptr);

    s_topov_masterelem_map[theTopo] = theElem;
  }
  else {
    theElem = it->second;
  }

  return theElem;
}

#ifndef KOKKOS_HAVE_CUDA

double initialize_linear_scalar_field(
  const stk::mesh::BulkData& bulk,
  const VectorFieldType& coordField,
  const ScalarFieldType& qField)
{
  // q = a + b^T x + 1/2 x^T H x
  std::mt19937 rng;
  rng.seed(0); // fixed seed
  std::uniform_real_distribution<double> coeff(-1.0, 1.0);

  double a  = coeff(rng);

  double b[3];
  for (unsigned j = 0; j < 3; ++j) {
    b[j] = coeff(rng);
  }

  const auto& meta = bulk.mesh_meta_data();
  EXPECT_EQ(meta.spatial_dimension(), 3u);

  const stk::mesh::Selector selector = meta.locally_owned_part() | meta.globally_shared_part();
  const auto& buckets = bulk.get_buckets(stk::topology::NODE_RANK, selector);
  kokkos_thread_team_bucket_loop(buckets, [&](stk::mesh::Entity node)
  {
    const double* coords = stk::mesh::field_data(coordField, node);
    *stk::mesh::field_data(qField, node) = linear(a, b, coords);
  });

  return std::sqrt(b[0] * b[0] + b[1]* b[1] + b[2] * b[2]);
}


double initialize_quadratic_scalar_field(
  const stk::mesh::BulkData& bulk,
  const VectorFieldType& coordField,
  const ScalarFieldType& qField)
{
  // q = a + b^T x + 1/2 x^T H x
  std::mt19937 rng;
  rng.seed(0); // fixed seed
  std::uniform_real_distribution<double> coeff(-1.0, 1.0);

  double a  = coeff(rng);

  double b[3];
  for (unsigned j = 0; j < 3; ++j) {
    b[j] = coeff(rng);
  }

  double H[9];
  for (unsigned j = 0; j < 9; ++j) {
    H[j] = coeff(rng);
  }
  H[3] = H[1];
  H[6] = H[2];
  H[7] = H[5];

  const auto& meta = bulk.mesh_meta_data();
  EXPECT_EQ(meta.spatial_dimension(), 3u);

  const stk::mesh::Selector selector = meta.locally_owned_part() | meta.globally_shared_part();
  const auto& buckets = bulk.get_buckets(stk::topology::NODE_RANK, selector);
  kokkos_thread_team_bucket_loop(buckets, [&](stk::mesh::Entity node)
  {
    const double* coords = stk::mesh::field_data(coordField, node);
    *stk::mesh::field_data(qField, node) = quadratic(a,b,H, coords);
  });

  double traceOfHessian = H[0] + H[4] + H[8];

  return (traceOfHessian);
}

}//namespace unit_test_utils

void Hex8Mesh::check_discrete_laplacian(double exactLaplacian)
{
   const stk::mesh::Selector selector = meta.locally_owned_part() & !meta.globally_shared_part();
   const stk::mesh::BucketVector& nodeBuckets = bulk.get_buckets(stk::topology::NODE_RANK, selector);
   kokkos_thread_team_bucket_loop(nodeBuckets, [&](stk::mesh::Entity node)
   {
     if (bulk.num_elements(node) == 8) {
         EXPECT_NEAR(*stk::mesh::field_data(*discreteLaplacianOfPressure, node), exactLaplacian    , tol);
     }
   });
}

#endif

