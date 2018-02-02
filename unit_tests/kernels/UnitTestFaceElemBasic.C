#include "UnitTestUtils.h"
#include "UnitTestHelperObjects.h"
#include <stk_util/parallel/Parallel.hpp>


#include "BcAlgTraits.h"
#include "ElemDataRequests.h"
#include "ScratchViews.h"
#include "SolutionOptions.h"
#include "CopyAndInterleave.h"

#include "MomentumSymmetryElemKernel.h"

#include <gtest/gtest.h>

namespace {
void fill_with_node_ids(stk::mesh::BulkData& bulk, IdFieldType* idField)
{
    const stk::mesh::BucketVector& nodeBuckets = bulk.buckets(stk::topology::NODE_RANK);
    for(const stk::mesh::Bucket* bptr : nodeBuckets) {
        for(stk::mesh::Entity node : *bptr) {
            double* idptr = stk::mesh::field_data(*idField, node);
            *idptr = bulk.identifier(node);
        }
    }
}

void verify_faces_exist(const stk::mesh::BulkData& bulk)
{
    EXPECT_TRUE(bulk.buckets(stk::topology::FACE_RANK).size() > 0);
}

class TestFaceElemKernel {
public:
    TestFaceElemKernel(stk::topology faceTopo, stk::topology elemTopo,
                       IdFieldType* idField,
                       sierra::nalu::ElemDataRequests& faceDataNeeded,
                       sierra::nalu::ElemDataRequests& elemDataNeeded)
    : numTimesExecuted_(0), faceTopo_(faceTopo), elemTopo_(elemTopo), idField_(idField)
    {
        faceDataNeeded.add_gathered_nodal_field(*idField, 1);
        elemDataNeeded.add_gathered_nodal_field(*idField, 1);
    }

    void execute(/*perhaps bundle args into something like SharedMemData... */
                 sierra::nalu::ScratchViews<DoubleType>& faceViews,
                 sierra::nalu::ScratchViews<DoubleType>& elemViews,
                 int numSimdFaces,
                 const int* elemFaceOrdinals)
    {
        sierra::nalu::SharedMemView<DoubleType*>& faceNodeIds = faceViews.get_scratch_view_1D(*idField_);
        sierra::nalu::SharedMemView<DoubleType*>& elemNodeIds = elemViews.get_scratch_view_1D(*idField_);
        EXPECT_EQ(faceTopo_.num_nodes(), faceNodeIds.size());
        EXPECT_EQ(elemTopo_.num_nodes(), elemNodeIds.size());

        std::vector<int> faceNodeOrdinals(faceTopo_.num_nodes());

        for(int simdIndex=0; simdIndex<numSimdFaces; ++simdIndex) {
           elemTopo_.face_node_ordinals(elemFaceOrdinals[simdIndex], faceNodeOrdinals.begin());
           for(unsigned i=0; i<faceTopo_.num_nodes(); ++i) {
               DoubleType faceNodeId = faceNodeIds(i);
               DoubleType elemNodeId = elemNodeIds(faceNodeOrdinals[i]);
               EXPECT_NEAR(stk::simd::get_data(faceNodeId,simdIndex), stk::simd::get_data(elemNodeId,simdIndex), 1.e-9);
           }
        }

        ++numTimesExecuted_;
    }

    unsigned numTimesExecuted_;
    stk::topology faceTopo_;
    stk::topology elemTopo_;
    IdFieldType* idField_;
};

constexpr int simdLen = stk::simd::ndoubles;

int
calculate_shared_mem_bytes_per_thread(int lhsSize, int rhsSize, int scratchIdsSize, int nDim,
                                      sierra::nalu::ElemDataRequests& faceDataNeeded,
                                      sierra::nalu::ElemDataRequests& elemDataNeeded,
                                      const sierra::nalu::ScratchMeInfo &meInfo)
{
    int bytes_per_thread = (rhsSize + lhsSize)*sizeof(double) + (2*scratchIdsSize)*sizeof(int)
                         + sierra::nalu::get_num_bytes_pre_req_data<double>(faceDataNeeded, nDim)
                         + sierra::nalu::get_num_bytes_pre_req_data<double>(elemDataNeeded, nDim, meInfo);
    bytes_per_thread *= 2*simdLen;
    return bytes_per_thread;
}

size_t get_simd_bucket_length(size_t bktLength)
{
    size_t simdBucketLen = bktLength/simdLen;
    const size_t remainder = bktLength%simdLen;
    if (remainder > 0) {
      simdBucketLen += 1;
    }
    return simdBucketLen;
}
int get_next_num_elems_simd(int bktIndex, int bktLength)
{
  int numElems = simdLen;
  if (bktLength - bktIndex*simdLen < simdLen) {
    numElems = bktLength - bktIndex*simdLen;
  }
  if (numElems < 0 || numElems > simdLen) {
    std::cout<<"ERROR, simdElems="<<numElems<<" shouldn't happen!!"<<std::endl;
    numElems = 0;
  }
  return numElems;
}

class TestFaceElemAlgorithm {
public:
    TestFaceElemAlgorithm(stk::mesh::Part* part, stk::mesh::EntityRank rank, unsigned nodesPerEntity, unsigned numDof=1)
    : kernels_(), part_(part), entityRank_(rank), nodesPerEntity_(nodesPerEntity), numDof_(numDof)
    {
    }

    template<typename LambdaFunction>
    void run_face_elem_algorithm(stk::mesh::BulkData& bulk, LambdaFunction func)
    {
        int nDim = bulk.mesh_meta_data().spatial_dimension();

        // hard code for quad4/hex8 element pairs
        sierra::nalu::ScratchMeInfo meElemInfo;
        meElemInfo.nodalGatherSize_ = 8;
        meElemInfo.nodesPerFace_ = 4;
        meElemInfo.nodesPerElement_ = 8;
        meElemInfo.numFaceIp_ = 4;
        meElemInfo.numScsIp_ = 12;
        meElemInfo.numScvIp_ = 8;

        int rhsSize = meElemInfo.nodalGatherSize_*numDof_, lhsSize = rhsSize*rhsSize, scratchIdsSize = rhsSize;

        const int bytes_per_team = 0;
        const int bytes_per_thread = calculate_shared_mem_bytes_per_thread(lhsSize, rhsSize, scratchIdsSize,
                                                                         nDim, faceDataNeeded_, elemDataNeeded_, meElemInfo);

        const bool interleaveMeViews = false;

        stk::mesh::Selector s_locally_owned_union = bulk.mesh_meta_data().locally_owned_part() & *part_;

        stk::mesh::BucketVector const& buckets = bulk.get_buckets(entityRank_, s_locally_owned_union );

        auto team_exec = sierra::nalu::get_team_policy(buckets.size(), bytes_per_team, bytes_per_thread);
        Kokkos::parallel_for(team_exec, [&](const sierra::nalu::TeamHandleType& team)
        {
          stk::mesh::Bucket & b = *buckets[team.league_rank()];

          ThrowAssertMsg(b.topology().num_nodes() == (unsigned)nodesPerEntity_,
                         "TestFaceElemAlgorithm expected nodesPerEntity_ = "
                         <<nodesPerEntity_<<", but b.topology().num_nodes() = "<<b.topology().num_nodes());

          std::unique_ptr<sierra::nalu::ScratchViews<double>> faceViews[simdLen];
          std::unique_ptr<sierra::nalu::ScratchViews<double>> elemViews[simdLen];

          for(int i=0; i<simdLen; ++i) {
              faceViews[i] = std::unique_ptr<sierra::nalu::ScratchViews<double>>(
                      new sierra::nalu::ScratchViews<double>(team, bulk, nodesPerEntity_, faceDataNeeded_));
              elemViews[i] = std::unique_ptr<sierra::nalu::ScratchViews<double>>(
                      new sierra::nalu::ScratchViews<double>(team, bulk, meElemInfo, elemDataNeeded_));
          }

          sierra::nalu::ScratchViews<DoubleType> simdFaceViews(team, bulk, nodesPerEntity_, faceDataNeeded_);
          sierra::nalu::ScratchViews<DoubleType> simdElemViews(team, bulk, meElemInfo, elemDataNeeded_);

          sierra::nalu::SharedMemView<DoubleType*> simdrhs = sierra::nalu::get_shmem_view_1D<DoubleType>(team, rhsSize);
          sierra::nalu::SharedMemView<DoubleType**> simdlhs = sierra::nalu::get_shmem_view_2D<DoubleType>(team, rhsSize, rhsSize);
          sierra::nalu::SharedMemView<double*> rhs = sierra::nalu::get_shmem_view_1D<double>(team, rhsSize);
          sierra::nalu::SharedMemView<double**> lhs = sierra::nalu::get_shmem_view_2D<double>(team, rhsSize, rhsSize);

          int elemFaceOrdinals[simdLen] = {-1};
          const size_t bucketLen   = b.size();
          const size_t simdBucketLen = get_simd_bucket_length(bucketLen);

          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, simdBucketLen), [&](const size_t& bktIndex)
          {
            int numSimdFaces = get_next_num_elems_simd(bktIndex, bucketLen);

            for(int simdFaceIndex=0; simdFaceIndex<numSimdFaces; ++simdFaceIndex) {
              stk::mesh::Entity face = b[bktIndex*simdLen + simdFaceIndex];
              elemFaceOrdinals[simdFaceIndex] = bulk.begin_element_ordinals(face)[0];
              sierra::nalu::fill_pre_req_data(faceDataNeeded_, bulk, face, *faceViews[simdFaceIndex], interleaveMeViews);

              const stk::mesh::Entity* elems = bulk.begin_elements(face);
              ThrowAssertMsg(bulk.num_elements(face)==1, "Expecting just 1 element attaced to face!");
              sierra::nalu::fill_pre_req_data(elemDataNeeded_, bulk, elems[0], *elemViews[simdFaceIndex], interleaveMeViews);
            }

            copy_and_interleave(faceViews, numSimdFaces, simdFaceViews, interleaveMeViews);
            copy_and_interleave(elemViews, numSimdFaces, simdElemViews, interleaveMeViews);
            fill_master_element_views(faceDataNeeded_, bulk, simdFaceViews, elemFaceOrdinals);
            fill_master_element_views(elemDataNeeded_, bulk, simdElemViews);

            func(/*other args here?*/ simdlhs, simdrhs, simdFaceViews, simdElemViews, numSimdFaces, elemFaceOrdinals);
          });
        });
    }

    sierra::nalu::ElemDataRequests faceDataNeeded_;
    sierra::nalu::ElemDataRequests elemDataNeeded_;
    std::vector<TestFaceElemKernel*> kernels_;

private:
    stk::mesh::Part* part_;
    stk::mesh::EntityRank entityRank_;
    unsigned nodesPerEntity_;
    unsigned numDof_;
};

} //anonymous namespace

TEST_F(Hex8Mesh, faceElemBasic)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }
  fill_mesh("generated:1x1x1|sideset:xXyYzZ");
  verify_faces_exist(bulk);
  fill_with_node_ids(bulk, idField);

  stk::topology faceTopo = stk::topology::QUAD_4;
  stk::topology elemTopo = stk::topology::HEX_8;
  sierra::nalu::MasterElement* meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(faceTopo);
  sierra::nalu::MasterElement* meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(elemTopo);
  sierra::nalu::MasterElement* meSCV = sierra::nalu::MasterElementRepo::get_volume_master_element(elemTopo);

  stk::mesh::Part* surface1 = meta.get_part("surface_1");
  TestFaceElemAlgorithm faceElemAlg(surface1, stk::topology::FACE_RANK, faceTopo.num_nodes());
  faceElemAlg.faceDataNeeded_.add_cvfem_face_me(meFC);
  faceElemAlg.elemDataNeeded_.add_cvfem_surface_me(meSCS);
  faceElemAlg.elemDataNeeded_.add_cvfem_volume_me(meSCV);

  TestFaceElemKernel faceElemKernel(faceTopo, elemTopo, idField,
                                    faceElemAlg.faceDataNeeded_, faceElemAlg.elemDataNeeded_);

  faceElemAlg.run_face_elem_algorithm(bulk,
          [&](sierra::nalu::SharedMemView<DoubleType**> &lhs,
              sierra::nalu::SharedMemView<DoubleType *> &rhs,
              sierra::nalu::ScratchViews<DoubleType>& faceViews,
              sierra::nalu::ScratchViews<DoubleType>& elemViews,
              int numSimdFaces,
              const int* elemFaceOrdinals)
      {
          faceElemKernel.execute(faceViews, elemViews, numSimdFaces, elemFaceOrdinals);
      });

  unsigned expectedNumFaces = 6;
  EXPECT_EQ(expectedNumFaces, faceElemKernel.numTimesExecuted_);
}

TEST_F(Hex8ElementWithBCFields, faceElemMomentumSymmetry)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 1) {
    return;
  }
  verify_faces_exist(bulk);

  sierra::nalu::SolutionOptions solnOptions;

  stk::topology faceTopo = stk::topology::QUAD_4;
  stk::mesh::Part* surface1 = meta.get_part("all_surfaces");
  TestFaceElemAlgorithm faceElemAlg(surface1, stk::topology::FACE_RANK, faceTopo.num_nodes(), sierra::nalu::BcAlgTraitsHex8Quad4::nDim_);

  sierra::nalu::MomentumSymmetryElemKernel<sierra::nalu::BcAlgTraitsHex8Quad4>  momentumSymmetryElemKernel(meta, solnOptions, &velocity, &viscosity,
                                    faceElemAlg.faceDataNeeded_, faceElemAlg.elemDataNeeded_);

  faceElemAlg.run_face_elem_algorithm(bulk,
          [&](
              sierra::nalu::SharedMemView<DoubleType**> &lhs,
              sierra::nalu::SharedMemView<DoubleType *> &rhs,
              sierra::nalu::ScratchViews<DoubleType>& faceViews,
              sierra::nalu::ScratchViews<DoubleType>& elemViews,
              int numSimdFaces,
              const int* elemFaceOrdinals)
      {
      momentumSymmetryElemKernel.execute(lhs, rhs, faceViews, elemViews);
      });
}
