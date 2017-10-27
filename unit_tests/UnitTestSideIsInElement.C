#include <gtest/gtest.h>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/SkinMesh.hpp>

#include <nalu_make_unique.h>
#include <NaluEnv.h>

#include <memory>
#include <tuple>
#include <random>
#include <limits>

#include "UnitTestUtils.h"



namespace
{
  using VectorFieldType = stk::mesh::Field<double, stk::mesh::Cartesian>;

  std::array<std::array<double,3>,3> random_rotation_matrix(int dim, std::mt19937& rng)
  {
    // always 3D, but set rotation around Z so only the relevant 2x2 subset matters
    std::uniform_real_distribution<double> angle(0, 2*M_PI);
    std::uniform_real_distribution<double> axis(0, 1);

    const double theta = angle(rng);
    const double cosTheta = cos(theta);
    const double sinTheta = sin(theta);
    std::array<double,3> n = {{(dim==3) ? axis(rng) : 0, (dim==3) ? axis(rng) : 0, (dim==3) ? axis(rng) : 1}};
    std::array<std::array<double,3>,3> rot = {{ {{cosTheta,0,0}},{{0,cosTheta,0}},{{0,0,cosTheta}} }};
    std::array<std::array<double,3>,3> nX = {{ {{0,0,0}},{{0,0,0}},{{0,0,0}} }};

    if( dim == 3 ) {
      double nMag = 0.0;
      for (int j = 0; j < 3; ++j) {
        nMag += n[j]*n[j];
      }
      nMag = 1.0/std::sqrt(nMag);
      for (int j = 0; j < dim; ++j) {
        n[j] *= nMag;
      }
    }

    nX[0][1] = -n[2];
    nX[0][2] =  n[1];
    nX[1][0] =  n[2];
    nX[1][2] = -n[0];
    nX[2][0] = -n[1];
    nX[2][1] =  n[0];

    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        rot[i][j] += (1-cosTheta)*n[i]*n[j] + sinTheta*nX[i][j];
      }
    }

    return rot;
  }

  void randomly_perturb_element_coords(
    std::mt19937& rng,
    int npe,
    const stk::mesh::Entity* node_rels,
    const VectorFieldType& coordField)
  {
    int dim = coordField.max_size(stk::topology::NODE_RANK);
    auto rot = random_rotation_matrix(dim, rng);

    std::uniform_real_distribution<double> random_perturb(0.125, 0.25);

    std::vector<double> coords_rt(npe * dim, 0.0);
    for (int n = 0; n < npe; ++n) {
      const double* coords =  stk::mesh::field_data(coordField, node_rels[n]);
      for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
          coords_rt[n*dim + i] += rot[i][j]*coords[j];
        }
      }
    }

    for (int n = 0; n < npe; ++n) {
      double* coords =  stk::mesh::field_data(coordField, node_rels[n]);
      for (int d = 0; d < dim; ++d) {
        coords[d] = coords_rt[n*dim+d] + random_perturb(rng);
      }
    }
  }


  void check_side_is_in_element(stk::topology topo)
  {
    if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) { return; }

    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::uniform_real_distribution<double> random_pt(-0.1,1.1); // can be outside of element side

    for (int ntrials = 0; ntrials < 25; ++ntrials) {

      std::vector<double> sideIsoPoint = { random_pt(rng), random_pt(rng) };

      int dim = topo.dimension();
      stk::mesh::MetaData meta(topo.dimension());
      stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD);

      auto elem = unit_test_utils::create_one_reference_element(bulk, topo);
      const VectorFieldType& coordField = *static_cast<const VectorFieldType*>(meta.coordinate_field());
      randomly_perturb_element_coords(rng, topo.num_nodes(), bulk.begin_nodes(elem), coordField);

      const auto& face_buckets = bulk.get_buckets(meta.side_rank(), meta.universal_part());
      for (const auto* ib : face_buckets) {
        const auto& b = *ib;

        auto* meSide = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());
        std::vector<double> sideInterpWeights(meSide->numIntPoints_ * meSide->nodesPerElement_, 0.0);
        meSide->shape_fcn(sideInterpWeights.data());

        for (size_t k = 0; k < b.size(); ++k) {
          auto face = b[k];
          const auto* face_node_rels = bulk.begin_nodes(face);

          std::vector<double> faceCoords(dim * meSide->nodesPerElement_);
          for (int d = 0; d < dim; ++d) {
            for (int n = 0; n < meSide->nodesPerElement_; ++n) {
              stk::mesh::Entity face_node = face_node_rels[n];
              double* coords =  stk::mesh::field_data(coordField, face_node);
              faceCoords[d*meSide->nodesPerElement_ + n] = coords[d];
            }
          }

          // interpolat to the parametric coordinates
          std::vector<double> faceInterpCoords(dim, 0.0);
          meSide->interpolatePoint(dim, sideIsoPoint.data(), faceCoords.data(), faceInterpCoords.data());

          // check if isInElement can reproduce the original parametric coords
          std::vector<double> isoFacePoint(dim-1 , 0.0);
          meSide->isInElement(faceCoords.data(), faceInterpCoords.data(), isoFacePoint.data());


          for (int d = 0; d < dim-1; ++d) {
            EXPECT_NEAR(isoFacePoint[d], sideIsoPoint[d], 1.0e-8);
          }
        }
      }
    }
  }

} // namespace

// Pyramids won't work.  Edge32D SCS has no isInElement implementation yet

#define TEST_ALL_VALID_TOPOS(x,y) \
    TEST(x, tri3##_##y)   { y(stk::topology::TRI_3_2D); } \
    TEST(x, quad4##_##y)  { y(stk::topology::QUAD_4_2D); } \
    TEST(x, tet4##_##y)   { y(stk::topology::TET_4); } \
    TEST(x, wedge6##_##y) { y(stk::topology::WEDGE_6); } \
    TEST(x, hex8##_##y)   { y(stk::topology::HEX_8); } \
    TEST(x, hex27##_##y)  { y(stk::topology::HEX_27); }

TEST_ALL_VALID_TOPOS(SideIsInElement, check_side_is_in_element);

