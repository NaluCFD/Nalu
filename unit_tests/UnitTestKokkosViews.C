#include <gtest/gtest.h>
#include <limits>
#include <random>

#include "UnitTestUtils.h"

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <Kokkos_Core.hpp>

#include <master_element/MasterElement.h>

namespace {

const double tol = 1.0e-10;

typedef stk::mesh::Field<double> ScalarFieldType;
typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType;

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

//========= kokkos helper stuff ==============================

typedef Kokkos::Schedule<Kokkos::Dynamic> DynamicScheduleType;
typedef typename Kokkos::TeamPolicy<typename Kokkos::DefaultExecutionSpace, DynamicScheduleType>::member_type TeamHandleType;

using DeviceShmem = Kokkos::DefaultExecutionSpace::scratch_memory_space;
template<typename T>
using SharedMemView = Kokkos::View<T, Kokkos::LayoutRight, DeviceShmem, Kokkos::MemoryUnmanaged>;
using DeviceTeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
using DeviceTeam = DeviceTeamPolicy::member_type;

inline DeviceTeamPolicy get_team_policy(const size_t sz, const size_t bytes_per_team,
    const size_t bytes_per_thread)
{
  DeviceTeamPolicy policy(sz, Kokkos::AUTO);
  return policy.set_scratch_size(0, Kokkos::PerTeam(bytes_per_team), Kokkos::PerThread(bytes_per_thread));
}

inline
SharedMemView<double*> get_shmem_view_1D(const TeamHandleType& team, size_t len)
{
  return Kokkos::subview(SharedMemView<double**>(team.team_shmem(), team.team_size(), len), team.team_rank(), Kokkos::ALL());
}

inline
SharedMemView<double**> get_shmem_view_2D(const TeamHandleType& team, size_t len1, size_t len2)
{
  return Kokkos::subview(SharedMemView<double***>(team.team_shmem(), team.team_size(), len1, len2), team.team_rank(), Kokkos::ALL(), Kokkos::ALL());
}

//========= end of kokkos helper stuff ==============================

//============= loop abstraction functions ==========================

template<class OUTER_LOOP_BODY, class INNER_LOOP_BODY>
void bucket_loop_serial_only(const stk::mesh::BucketVector& buckets, const OUTER_LOOP_BODY& outer_loop_body, const INNER_LOOP_BODY& inner_loop_body)
{
    for(const stk::mesh::Bucket* bptr : buckets)
    {
        const stk::mesh::Bucket& bkt = *bptr;
        stk::topology topo = bkt.topology();
        sierra::nalu::MasterElement* meSCS = get_surface_master_element(topo);

        outer_loop_body(topo,*meSCS);

        for(size_t j=0; j<bkt.size(); ++j)
        {
            inner_loop_body(bkt[j], topo, *meSCS);
        }
    }
}

template<class LOOP_BODY>
void kokkos_bucket_loop(const stk::mesh::BucketVector& buckets, LOOP_BODY inner_loop_body)
{
    Kokkos::parallel_for(buckets.size(), [&](const size_t& i)
    {
        const stk::mesh::Bucket& bkt = *buckets[i];
        for(size_t j=0; j<bkt.size(); ++j)
        {
            inner_loop_body(bkt[j]);
        }
    });
}

template<class LOOP_BODY>
void kokkos_thread_team_bucket_loop(const stk::mesh::BucketVector& buckets, LOOP_BODY inner_loop_body)
{
    Kokkos::parallel_for(Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>(buckets.size(), Kokkos::AUTO), KOKKOS_LAMBDA(const TeamHandleType& team)
    {
        const stk::mesh::Bucket& bkt = *buckets[team.league_rank()];
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, bkt.size()), [&](const size_t& j)
        {
            inner_loop_body(bkt[j]);
        });
    });
}

template<class LOOP_BODY>
void kokkos_thread_team_bucket_loop_with_topo(const stk::mesh::BucketVector& buckets,
                                    const LOOP_BODY& inner_loop_body)
{
    Kokkos::parallel_for(Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>(buckets.size(), Kokkos::AUTO), KOKKOS_LAMBDA(const TeamHandleType& team)
    {
        const stk::mesh::Bucket& bkt = *buckets[team.league_rank()];
        stk::topology topo = bkt.topology();
        sierra::nalu::MasterElement* meSCS = get_surface_master_element(topo);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, bkt.size()), [&](const size_t& j)
        {
            inner_loop_body(bkt[j], topo, *meSCS);
        });
    });
}

//============= end of loop abstraction functions ==========================

double quadratic(double a, const double* b, const double* H, const double* x)
{
  double lin = b[0]*x[0] + b[1]*x[1] + b[2]*x[2];
  double quad = x[0] * (H[0]*x[0] + H[1]*x[1] + H[2]*x[2])
              + x[1] * (H[3]*x[0] + H[4]*x[1] + H[5]*x[2])
              + x[2] * (H[6]*x[0] + H[7]*x[1] + H[8]*x[2]);

  return (a + lin + 0.5*quad);
}

#ifndef KOKKOS_HAVE_CUDA
//following tests can't run on cuda due to variety of reasons, including
//use of std::vectors, use of MasterElement functions (defined for host), etc.

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

  double H[9];
  for (unsigned j = 0; j < 9; ++j) {
    H[j] = coeff(rng);
  }

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

    void fill_mesh(const std::string& meshSpec = "generated:70x70x70")
    {
      unit_test_utils::fill_hex8_mesh(meshSpec, bulk);
    }

    void fill_mesh_and_initialize_test_fields(const std::string& meshSpec = "generated:70x70x70")
    {
        fill_mesh(meshSpec);

        partVec = {&meta.locally_owned_part()};

        coordField = static_cast<const VectorFieldType*>(meta.coordinate_field());
        EXPECT_TRUE(coordField != nullptr);

        exactLaplacian = initialize_linear_scalar_field(bulk, *coordField, *nodalPressureField);
        stk::mesh::field_fill(0.0, *discreteLaplacianOfPressure);
    }

    void check_discrete_laplacian(double exactLaplacian)
    {
       const stk::mesh::Selector selector = meta.locally_owned_part() & !meta.globally_shared_part();
       const stk::mesh::BucketVector& nodeBuckets = bulk.get_buckets(stk::topology::NODE_RANK, selector);
       kokkos_thread_team_bucket_loop(nodeBuckets, [&](stk::mesh::Entity node)
       {
         // we didn't include parallel communication or boundary stencil modification, so
         // only expect that the Laplacian calculation works at regularly connected nodes
         // in the domain
         if (bulk.num_elements(node) == 8) {
           EXPECT_NEAR(*stk::mesh::field_data(*discreteLaplacianOfPressure, node), exactLaplacian, tol);
         }
       });
    }

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

void find_max_nodes_and_ips(const stk::mesh::BucketVector& buckets,
                            int& maxNodesPerElement, int& maxScsIp)
{
  maxNodesPerElement = 0;
  maxScsIp = 0;
  size_t numEntities = 0;
  for(const stk::mesh::Bucket* bptr : buckets) {
    stk::topology topo = bptr->topology();
    maxNodesPerElement = std::max(maxNodesPerElement, (int)topo.num_nodes());
    sierra::nalu::MasterElement *meSCS = get_surface_master_element(topo);
    maxScsIp = std::max(maxScsIp, meSCS->numIntPoints_);
    numEntities += bptr->size();
  }
  std::cout<<"num entities: "<<numEntities<<std::endl;
}

//=========== Test class that mimics an element algorithm ==============
//=========== and uses std::vectors for scratch arrays    ==============
//
class TestElemAlgorithmWithVectors
{
public:
  TestElemAlgorithmWithVectors(stk::mesh::BulkData& bulk, const stk::mesh::PartVector& partVec,
                    const VectorFieldType* coord, ScalarFieldType* discreteLaplacian,
                    ScalarFieldType* nodalPressure)
  : bulkData_(bulk), partVec_(partVec),
    discreteLaplacianOfPressure(discreteLaplacian),
    nodalPressureField(nodalPressure), coordField(coord)
  {}

  void execute()
  {
      double scs_error = 0.0;
      const stk::mesh::MetaData& meta = bulkData_.mesh_meta_data();
  
      const int nDim = meta.spatial_dimension();
  
      std::vector<double> elemNodeCoords;
      std::vector<double> elemNodePressures;
  
      std::vector<double> scs_areav;
      std::vector<double> dndx;
      std::vector<double> deriv;
      std::vector<double> det_j;
  
      auto resizer = [&](int nodesPerElem, int numScsIp)
      {
          elemNodeCoords.resize(nodesPerElem*nDim);
          elemNodePressures.resize(nodesPerElem);
          scs_areav.resize(numScsIp*nDim);
          dndx.resize(nDim*numScsIp*nodesPerElem);
          deriv.resize(nDim*numScsIp*nodesPerElem);
          det_j.resize(numScsIp);
      };

      const stk::mesh::BucketVector& elemBuckets = bulkData_.get_buckets(stk::topology::ELEM_RANK, meta.locally_owned_part());
  
      int maxNodesPerElement = 0, maxNumScsIp = 0;
      find_max_nodes_and_ips(elemBuckets, maxNodesPerElement, maxNumScsIp);
      resizer(maxNodesPerElement, maxNumScsIp);

      bucket_loop_serial_only(elemBuckets,
          [&](stk::topology topo, sierra::nalu::MasterElement& meSCS)
          {
              const int nodesPerElem = topo.num_nodes();
              resizer(nodesPerElem, meSCS.numIntPoints_);
          }
          ,
          [&](stk::mesh::Entity elem, stk::topology topo, sierra::nalu::MasterElement& meSCS)
          {
              const stk::mesh::Entity* elemNodes = bulkData_.begin_nodes(elem);
      
              double* p_elemNodeCoords = elemNodeCoords.data();
              double* p_elemNodePressures = elemNodePressures.data();
      
              double* p_scs_areav = scs_areav.data();
              double* p_dndx = dndx.data();
              double* p_deriv = deriv.data();
              double* p_det_j = det_j.data();
              const int* lrscv = meSCS.adjacentNodes();
      
              const int numScsIp = meSCS.numIntPoints_;
              const int nodesPerElem = topo.num_nodes();
              for(int n=0; n<nodesPerElem; ++n) {
                  const double* nodeCoords = stk::mesh::field_data(*coordField, elemNodes[n]);
      
                  const int nodeOffset = n*nDim;
                  for(int d=0; d<nDim; ++d) {
                      p_elemNodeCoords[nodeOffset+d] = nodeCoords[d];
                  }
                  const double* nodePressure = stk::mesh::field_data(*nodalPressureField, elemNodes[n]);
                  p_elemNodePressures[n] = nodePressure[0];
              }
      
              meSCS.determinant(1, p_elemNodeCoords, p_scs_areav, &scs_error);
              meSCS.grad_op(1, p_elemNodeCoords, p_dndx, p_deriv, p_det_j, &scs_error);
      
              for (int ip = 0; ip < numScsIp; ++ip ) {
      
                double dpdxIp = 0.0;
                const int ipOffset = nDim*nodesPerElem*ip;
                for ( int ic = 0; ic < nodesPerElem; ++ic) {
                  const int offSetDnDx = ipOffset + ic*nDim;
                  for ( int j = 0; j < nDim; ++j ) {
                    dpdxIp += p_dndx[offSetDnDx+j]*p_elemNodePressures[ic]*p_scs_areav[ip*nDim+j];
                  }
                }
                EXPECT_TRUE(std::abs(dpdxIp) > tol);
      
                const stk::mesh::Entity lNode = elemNodes[lrscv[2*ip+0]];
                const stk::mesh::Entity rNode = elemNodes[lrscv[2*ip+1]];
      
                Kokkos::atomic_add(stk::mesh::field_data(*discreteLaplacianOfPressure, lNode), dpdxIp);
                Kokkos::atomic_add(stk::mesh::field_data(*discreteLaplacianOfPressure, rNode), -dpdxIp);
              }
          }
      );
  }

private:
  stk::mesh::BulkData& bulkData_;
  const stk::mesh::PartVector& partVec_;
  ScalarFieldType* discreteLaplacianOfPressure;
  ScalarFieldType* nodalPressureField;
  const VectorFieldType* coordField;
};


//======= templated element kernel function ==================

template<int nodesPerElem, int numScsIp>
void element_discrete_laplacian_kernel_3d(stk::mesh::BulkData& bulkData, stk::mesh::Entity elem,
                       sierra::nalu::MasterElement& meSCS,
                       ScalarFieldType* discreteLaplacianOfPressure,
                       ScalarFieldType* nodalPressureField,
                       const VectorFieldType* coordField)
{
    const int nDim = 3;
    const stk::mesh::Entity* elemNodes = bulkData.begin_nodes(elem);

    double p_elemNodeCoords[nodesPerElem*nDim];
    double p_elemNodePressures[nodesPerElem];

    double p_scs_areav[numScsIp*nDim];
    double p_dndx[nDim*numScsIp*nodesPerElem];
    double p_deriv[nDim*numScsIp*nodesPerElem];
    double p_det_j[numScsIp];
    const int* lrscv = meSCS.adjacentNodes();

    for(int n=0; n<nodesPerElem; ++n) {
        const double* nodeCoords = stk::mesh::field_data(*coordField, elemNodes[n]);

        const int nodeOffset = n*nDim;
        for(int d=0; d<nDim; ++d) {
            p_elemNodeCoords[nodeOffset+d] = nodeCoords[d];
        }
        const double* nodePressure = stk::mesh::field_data(*nodalPressureField, elemNodes[n]);
        p_elemNodePressures[n] = nodePressure[0];
    }

    double scs_error = 0;
    meSCS.determinant(1, p_elemNodeCoords, p_scs_areav, &scs_error);
    meSCS.grad_op(1, p_elemNodeCoords, p_dndx, p_deriv, p_det_j, &scs_error);

    for (int ip = 0; ip < numScsIp; ++ip ) {

      double dpdxIp = 0.0;
      const int ipOffset = nDim*nodesPerElem*ip;
      for ( int ic = 0; ic < nodesPerElem; ++ic) {
        const int offSetDnDx = ipOffset + ic*nDim;
        for ( int j = 0; j < nDim; ++j ) {
          dpdxIp += p_dndx[offSetDnDx+j]*p_elemNodePressures[ic]*p_scs_areav[ip*nDim+j];
        }
      }
      EXPECT_TRUE(std::abs(dpdxIp) > tol);

      const stk::mesh::Entity lNode = elemNodes[lrscv[2*ip+0]];
      const stk::mesh::Entity rNode = elemNodes[lrscv[2*ip+1]];

      Kokkos::atomic_add(stk::mesh::field_data(*discreteLaplacianOfPressure, lNode), dpdxIp);
      Kokkos::atomic_add(stk::mesh::field_data(*discreteLaplacianOfPressure, rNode), -dpdxIp);
    }
}

//=========== Test class that mimics an element algorithm ==============
//=========== and uses a templated kernel (with compile-time scratch arrays)  ==============
//
class TestElemAlgorithmWithTemplate
{
public:
  TestElemAlgorithmWithTemplate(stk::mesh::BulkData& bulk, const stk::mesh::PartVector& partVec,
                    const VectorFieldType* coord, ScalarFieldType* discreteLaplacian,
                    ScalarFieldType* nodalPressure)
  : bulkData_(bulk), partVec_(partVec),
    discreteLaplacianOfPressure(discreteLaplacian),
    nodalPressureField(nodalPressure), coordField(coord)
  {}

  void execute()
  {
      const stk::mesh::MetaData& meta = bulkData_.mesh_meta_data();
  
      const stk::mesh::BucketVector& elemBuckets = bulkData_.get_buckets(stk::topology::ELEM_RANK, meta.locally_owned_part());
  
      kokkos_thread_team_bucket_loop_with_topo(elemBuckets,
          [&](stk::mesh::Entity elem, stk::topology topo, sierra::nalu::MasterElement& meSCS)
          {
             //this is an incomplete switch, doesn't handle all possible topologies...
             //just an illustration for this test.
              switch(topo) {
              case stk::topology::HEX_8:
                  element_discrete_laplacian_kernel_3d<8,12>(bulkData_, elem, meSCS,
                       discreteLaplacianOfPressure, nodalPressureField, coordField);
                  break;
              case stk::topology::HEX_27:
                  element_discrete_laplacian_kernel_3d<27,216>(bulkData_, elem, meSCS,
                       discreteLaplacianOfPressure, nodalPressureField, coordField);
                  break;
              case stk::topology::TET_4:
                  element_discrete_laplacian_kernel_3d<4,6>(bulkData_, elem, meSCS,
                       discreteLaplacianOfPressure, nodalPressureField, coordField);
                  break;
              case stk::topology::PYRAMID_5:
                  element_discrete_laplacian_kernel_3d<5,8>(bulkData_, elem, meSCS,
                       discreteLaplacianOfPressure, nodalPressureField, coordField);
                  break;
              case stk::topology::WEDGE_6:
                  element_discrete_laplacian_kernel_3d<6,9>(bulkData_, elem, meSCS,
                       discreteLaplacianOfPressure, nodalPressureField, coordField);
                  break;
              default:
                  std::cerr<<"ERROR! Unhandled topology: "<<topo<<std::endl;
                  break;
              }
          }
      );
  }

private:
  stk::mesh::BulkData& bulkData_;
  const stk::mesh::PartVector& partVec_;
  ScalarFieldType* discreteLaplacianOfPressure;
  ScalarFieldType* nodalPressureField;
  const VectorFieldType* coordField;
};

//=========== Test class that mimics an element algorithm ==============
//=========== and uses Kokkos::views for scratch arrays    ==============
//
class TestElemAlgorithmWithViews
{
public:
  TestElemAlgorithmWithViews(stk::mesh::BulkData& bulk, const stk::mesh::PartVector& partVec,
                    const VectorFieldType* coord, ScalarFieldType* discreteLaplacian,
                    ScalarFieldType* nodalPressure)
  : bulkData_(bulk), partVec_(partVec),
    discreteLaplacianOfPressure(discreteLaplacian),
    nodalPressureField(nodalPressure), coordField(coord)
  {}

  void execute()
  {
    double scs_error = 0.0;
    const stk::mesh::MetaData& meta = bulkData_.mesh_meta_data();

    const int nDim = meta.spatial_dimension();

    const stk::mesh::BucketVector& elemBuckets = bulkData_.get_buckets(stk::topology::ELEM_RANK, meta.locally_owned_part());

    int maxNodesPerElement = 0, maxNumScsIp = 0;
    find_max_nodes_and_ips(elemBuckets, maxNodesPerElement, maxNumScsIp);

    const int bytes_per_team = 0;
    const int bytes_per_thread =
       SharedMemView<double**>::shmem_size(maxNodesPerElement, nDim) +
       SharedMemView<double*>::shmem_size(maxNodesPerElement) +
       SharedMemView<double**>::shmem_size(maxNumScsIp, nDim) +
       SharedMemView<double**>::shmem_size(maxNumScsIp, maxNodesPerElement*nDim) +
       SharedMemView<double**>::shmem_size(maxNumScsIp, maxNodesPerElement*nDim) +
       SharedMemView<double*>::shmem_size(maxNumScsIp);

    auto team_exec = get_team_policy(elemBuckets.size(), bytes_per_team, bytes_per_thread);
    Kokkos::parallel_for(team_exec, [&](const TeamHandleType& team)
    {
        const stk::mesh::Bucket& bkt = *elemBuckets[team.league_rank()];
        stk::topology topo = bkt.topology();
        sierra::nalu::MasterElement& meSCS = *get_surface_master_element(topo);

        const int nodesPerElem = topo.num_nodes();
        const int numScsIp = meSCS.numIntPoints_;

        SharedMemView<double**> elemNodeCoords;
        SharedMemView<double*> elemNodePressures;
     
        SharedMemView<double**> scs_areav;
        SharedMemView<double**> dndx;
        SharedMemView<double**> deriv;
        SharedMemView<double*> det_j;

        elemNodeCoords = get_shmem_view_2D(team, nodesPerElem, nDim);
        elemNodePressures = get_shmem_view_1D(team, nodesPerElem);
        scs_areav = get_shmem_view_2D(team, numScsIp, nDim);
        dndx = get_shmem_view_2D(team, numScsIp, nodesPerElem*nDim);
        deriv = get_shmem_view_2D(team, numScsIp, nodesPerElem*nDim);
        det_j = get_shmem_view_1D(team, numScsIp);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, bkt.size()), [&](const size_t& jj)
        {
             stk::mesh::Entity elem = bkt[jj];
              const stk::mesh::Entity* elemNodes = bulkData_.begin_nodes(elem);
              for(int n=0; n<nodesPerElem; ++n) {
                  const double* nodeCoords = stk::mesh::field_data(*coordField, elemNodes[n]);
      
                  for(int d=0; d<nDim; ++d) {
                      elemNodeCoords(n,d) = nodeCoords[d];
                  }
                  const double* nodePressure = stk::mesh::field_data(*nodalPressureField, elemNodes[n]);
                  elemNodePressures[n] = nodePressure[0];
              }
      
              meSCS.determinant(1, &elemNodeCoords(0,0), &scs_areav(0,0), &scs_error);
              meSCS.grad_op(1, &elemNodeCoords(0,0), &dndx(0,0), &deriv(0,0), &det_j(0), &scs_error);
              const int* lrscv = meSCS.adjacentNodes();
      
              for (int ip = 0; ip < numScsIp; ++ip ) {
      
                double dpdxIp = 0.0;
                for ( int ic = 0; ic < nodesPerElem; ++ic) {
                  for ( int j = 0; j < nDim; ++j ) {
                    dpdxIp += dndx(ip,ic*nDim+j)*elemNodePressures(ic)*scs_areav(ip,j);
                  }
                }
                EXPECT_TRUE(std::abs(dpdxIp) > tol);
      
                const stk::mesh::Entity lNode = elemNodes[lrscv[2*ip+0]];
                const stk::mesh::Entity rNode = elemNodes[lrscv[2*ip+1]];
      
                Kokkos::atomic_add(stk::mesh::field_data(*discreteLaplacianOfPressure, lNode), dpdxIp);
                Kokkos::atomic_add(stk::mesh::field_data(*discreteLaplacianOfPressure, rNode), -dpdxIp);
              }
        });
    });
  }

private:
  stk::mesh::BulkData& bulkData_;
  const stk::mesh::PartVector& partVec_;
  ScalarFieldType* discreteLaplacianOfPressure;
  ScalarFieldType* nodalPressureField;
  const VectorFieldType* coordField;
};

//========= below are the test 'main's... ===============

TEST_F(Hex8Mesh, indexing_vectors)
{
    fill_mesh_and_initialize_test_fields();

    TestElemAlgorithmWithVectors testAlgorithm(bulk, partVec, coordField,
                          discreteLaplacianOfPressure, nodalPressureField);

    testAlgorithm.execute();

    check_discrete_laplacian(exactLaplacian);
}

TEST_F(Hex8Mesh, indexing_template_raw_arrays)
{
    fill_mesh_and_initialize_test_fields();

    TestElemAlgorithmWithTemplate testAlgorithm(bulk, partVec, coordField,
                          discreteLaplacianOfPressure, nodalPressureField);

    testAlgorithm.execute();

    check_discrete_laplacian(exactLaplacian);
}

TEST_F(Hex8Mesh, indexing_views)
{
    fill_mesh_and_initialize_test_fields();

    TestElemAlgorithmWithViews testAlgorithm(bulk, partVec, coordField,
                          discreteLaplacianOfPressure, nodalPressureField);

    testAlgorithm.execute();

    check_discrete_laplacian(exactLaplacian);
}

//end of stuff that's ifndef'd for KOKKOS_HAVE_CUDA
#endif

}

