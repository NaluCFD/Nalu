#include "UnitTestUtils.h"
#include "UnitTestHelperObjects.h"
#include <stk_util/parallel/Parallel.hpp>


#include "BcAlgTraits.h"
#include "ElemDataRequests.h"
#include "ScratchViews.h"
#include "SolutionOptions.h"
#include "CopyAndInterleave.h"

#include "AssembleFaceElemSolverAlgorithm.h"
#include "kernel/MomentumSymmetryElemKernel.h"

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
  unsigned numDof = 3;

  unit_test_utils::HelperObjects helperObjs(bulk, elemTopo, numDof, surface1);

  sierra::nalu::AssembleFaceElemSolverAlgorithm faceElemAlg(helperObjs.realm, surface1, &helperObjs.eqSystem,
                                                          faceTopo.num_nodes(), elemTopo.num_nodes());
  faceElemAlg.faceDataNeeded_.add_cvfem_face_me(meFC);
  faceElemAlg.elemDataNeeded_.add_cvfem_surface_me(meSCS);
  faceElemAlg.elemDataNeeded_.add_cvfem_volume_me(meSCV);

  TestFaceElemKernel faceElemKernel(faceTopo, elemTopo, idField,
                                    faceElemAlg.faceDataNeeded_, faceElemAlg.elemDataNeeded_);

  faceElemAlg.run_face_elem_algorithm(bulk,
          [&](sierra::nalu::SharedMemData_FaceElem &smdata)
      {
          faceElemKernel.execute(smdata.simdFaceViews, smdata.simdElemViews, smdata.numSimdFaces, smdata.elemFaceOrdinals);
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
  stk::topology elemTopo = stk::topology::HEX_8;
  stk::mesh::Part* surface1 = meta.get_part("all_surfaces");
  unit_test_utils::HelperObjects helperObjs(bulk, elemTopo, sierra::nalu::BcAlgTraitsHex8Quad4::nDim_, surface1);

  sierra::nalu::AssembleFaceElemSolverAlgorithm faceElemAlg(helperObjs.realm, surface1, &helperObjs.eqSystem,
                                                          faceTopo.num_nodes(), elemTopo.num_nodes());

  auto  momentumSymmetryElemKernel =
      new sierra::nalu::MomentumSymmetryElemKernel<sierra::nalu::BcAlgTraitsHex8Quad4>(meta, solnOptions, &velocity, &viscosity,
                                    faceElemAlg.faceDataNeeded_, faceElemAlg.elemDataNeeded_);

  faceElemAlg.activeKernels_.push_back(momentumSymmetryElemKernel);

  faceElemAlg.execute();
}
